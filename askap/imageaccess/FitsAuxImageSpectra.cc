/// @file FitsAuxImageSpectra.cc
/// @brief Access FITS image
/// @details This class implements IImageAccess interface for FITS image
///
///
/// @copyright (c) 2016 CSIRO
/// Australia Telescope National Facility (ATNF)
/// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
/// PO Box 76, Epping NSW 1710, Australia
/// atnf-enquiries@csiro.au
///
/// This file is part of the ASKAP software distribution.
///
/// The ASKAP software distribution is free software: you can redistribute it
/// and/or modify it under the terms of the GNU General Public License as
/// published by the Free Software Foundation; either version 2 of the License,
/// or (at your option) any later version.
///
/// This program is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with this program; if not, write to the Free Software
/// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
///
/// @author 
///

#include <askap/imageaccess/FitsAuxImageSpectra.h>
#include <askap/imageaccess/Utils.h>
#include <askap/askap/AskapLogging.h>

using namespace askap;
using namespace askap::accessors;

const static int FITS_BIN_TABLE_COMPONENT_NAME_MAX_LEN = 128;

ASKAP_LOGGER(logger, ".FitsAuxImageSpectra");

/// class static function
void FitsAuxImageSpectra::PrintError(int status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/

    char status_str[FLEN_STATUS];
    fits_get_errstatus(status, status_str);

    if (status) {
        ASKAPLOG_ERROR_STR(logger, "FitsIO error: " << status_str); /* print error report */

        exit(status);      /* terminate the program, returning error status */
    }
}
FitsAuxImageSpectra::FitsAuxImageSpectra(const std::string& fitsFileName,
                                         const int nChannels, const int nrows,
                                         const casacore::CoordinateSystem& coord,
                                         const casacore::RecordInterface &tableInfo)
    : itsStatus(0), itsName(fitsFileName), itsNChannels(nChannels),
      itsCurrentRow(1), itsIA(new FitsImageAccess())
{
    // we want to use image accessor to create the fits file and populate
    // the keywords for us.
    casacore::IPosition dummyShape(4,1,1,1,1);
    itsIA->create(itsName,dummyShape,coord);

    // image accessor creates the file with the .fits extension even if the filename
    // does not have the .fits extension
    if ( fitsFileName.rfind(".fits") == std::string::npos ) {
        itsName.append(".fits");
    }

    auto index = itsName.rfind(".fits");
    itsBasename = std::string(itsName,0,index);

    // open the file to create the binary table 
    if ( fits_open_file(&itsFitsPtr,itsName.c_str(), READWRITE, &itsStatus) )
         PrintError(itsStatus);
    
    // create the fits binary table
    create(tableInfo,nChannels,nrows);

    if ( fits_close_file(itsFitsPtr,&itsStatus) )       /* close the FITS file */
         PrintError(itsStatus);
}

void 
FitsAuxImageSpectra::create(const casacore::RecordInterface &tableInfo,
                            const int nChannels, const int nrows)
{
    const char extname[] = "PolSpec";

    // check to see if the user specifies the column name for the stoke
    casacore::String stokes("spectrum"); // spectrum for I,Q,U and V
    if ( tableInfo.isDefined("Stokes") ) {
        casacore::RecordFieldId fieldId(tableInfo.fieldNumber("Stokes"));
        tableInfo.get(fieldId,stokes);
    }
    unsigned int nCol = 2;
    CPointerWrapper cPointerWrapper { nCol };
    std::string col1 = "Id";
    cPointerWrapper.itsTType[0] = new char[sizeof(char) * col1.length() + 1];
    std::memset(cPointerWrapper.itsTType[0],'\0',col1.length() + 1);
    std::memcpy(cPointerWrapper.itsTType[0],col1.data(),col1.length());
    std::string compMaxLenTForm = std::to_string(FITS_BIN_TABLE_COMPONENT_NAME_MAX_LEN) + 'A';
    cPointerWrapper.itsTForm[0] = new char[sizeof(char)*compMaxLenTForm.length()+1];
    std::memset(cPointerWrapper.itsTForm[0],'\0',compMaxLenTForm.length()+1);
    std::memcpy(cPointerWrapper.itsTForm[0],compMaxLenTForm.data(),compMaxLenTForm.length());
    cPointerWrapper.itsUnits[0] = new char[2];
    std::memset(cPointerWrapper.itsUnits[0],'\0',2);
    std::memcpy(cPointerWrapper.itsUnits[0],"",1);

    std::string col2(stokes.data());
    cPointerWrapper.itsTType[1] = new char[sizeof(char) * col2.length() + 1];
    std::memset(cPointerWrapper.itsTType[1],'\0',col2.length() + 1);
    std::memcpy(cPointerWrapper.itsTType[1],col2.data(),col2.length());
    // Dont know why QE does not work on galaxy. It always shows QE(1)
    // instead of QE(288) in fv
    //std::string s = std::string("QE(") + std::to_string(nChannels) + ")";
    std::string s = std::to_string(nChannels) +  std::string("E");
    cPointerWrapper.itsTForm[1] = new char[s.length()+1];
    std::memset(cPointerWrapper.itsTForm[1],'\0',s.length()+1);
    std::memcpy(cPointerWrapper.itsTForm[1],s.data(),s.length());
    cPointerWrapper.itsUnits[1] = new char[2];
    std::memset(cPointerWrapper.itsUnits[1],'\0',2);
    std::memcpy(cPointerWrapper.itsUnits[1],"",1);

    
    itsStatus = 0;

    int tfields = 2;
    auto t = cPointerWrapper.itsTType;
    auto f = cPointerWrapper.itsTForm;
    auto u = cPointerWrapper.itsUnits;

//char *tform[] = { "a8",     "288E"};
    if ( fits_create_tbl(itsFitsPtr, BINARY_TBL, nrows, tfields,t,f,u,extname,&itsStatus) )
         PrintError(itsStatus);    


    // Get the values of CRVAL4 and CDELT4 which are the reference and delta frequency respectively
    std::pair<std::string,std::string> crval4 = itsIA->getMetadataKeyword(itsBasename,"CRVAL4");
    std::pair<std::string,std::string> cdelt4 = itsIA->getMetadataKeyword(itsBasename,"CDELT4");
    // also get the CRPIX4
    std::pair<std::string,std::string> crpix4 = itsIA->getMetadataKeyword(itsBasename,"CRPIX4");
    // Add these values to the binary table header as 
    double refFreq = 0.0;
    double deltFreq = 0.0;
    double crpix = 0.0;
    if ( crval4.first != "" ) {
        refFreq = std::stod(crval4.first);
    } else {
        ASKAPLOG_WARN_STR(logger,"CRVAL4 keyword is empty in " << itsName);
    }
    if ( cdelt4.first != "" ) {
        deltFreq = std::stod(cdelt4.first);
    } else {
        ASKAPLOG_WARN_STR(logger,"CDELT4 keyword is empty in " << itsName);
    }
    if ( crpix4.first != "" ) {
        crpix = std::stod(crpix4.first);
    } else {
        ASKAPLOG_WARN_STR(logger,"CRPIX4 keyword is empty in " << itsName);
    }

    int nchan = nChannels;
    if ( fits_update_key(itsFitsPtr, TINT, "NCHAN", &nchan, "Number of Channels", &itsStatus) )
        PrintError(itsStatus);

    if ( fits_update_key(itsFitsPtr, TDOUBLE, "REF_FREQ", &refFreq, crval4.second.c_str(), &itsStatus) )
        PrintError(itsStatus);

    if ( fits_update_key(itsFitsPtr, TDOUBLE, "FREQ_INC", &deltFreq, cdelt4.second.c_str(), &itsStatus) )
        PrintError(itsStatus);
    
    if ( fits_update_key(itsFitsPtr, TDOUBLE, "REF_FPIX", &crpix, crpix4.second.c_str(), &itsStatus) )
        PrintError(itsStatus);
}

void 
FitsAuxImageSpectra::add(const std::string& id, const SpectrumT& spectrum)
{
    itsStatus = 0;

    // open the file again
    if ( fits_open_file(&itsFitsPtr,itsName.c_str(), READWRITE, &itsStatus) )
         PrintError(itsStatus);

    // move to the HDU 2 which is the spectrum binary table
    int hduType;
    if ( fits_movabs_hdu(itsFitsPtr,spectrumHDU(),&hduType,&itsStatus) )
         PrintError(itsStatus);

    long firstElem = 1;
    char* bptr[] = {const_cast<char *> (id.c_str())};
    ASKAPCHECK(id.length() <= FITS_BIN_TABLE_COMPONENT_NAME_MAX_LEN,"column name (ie component name) exceeded 128 characters");
    if (fits_write_col(itsFitsPtr,TSTRING,1,itsCurrentRow,firstElem,1,
                       bptr,&itsStatus) )
        PrintError(itsStatus);

    if (fits_write_col(itsFitsPtr,TFLOAT,2,itsCurrentRow,firstElem,spectrum.size(),
                       const_cast<float *> (spectrum.data()),&itsStatus) )
        PrintError(itsStatus);

    // close the FITS file
    if (fits_close_file(itsFitsPtr,&itsStatus))
         PrintError(itsStatus);

    itsId2RowMap.insert(std::make_pair(id,itsCurrentRow));
    itsCurrentRow += 1;
    
}

void FitsAuxImageSpectra::add(const std::vector<std::string>& ids, 
                              const ArrayOfSpectrumT& arrayOfSpectrums)
{
    ASKAPLOG_INFO_STR(logger,"---> FitsAuxImageSpectra::add");
    auto s1 = ids.size();
    auto s2 = arrayOfSpectrums.nrow();
    ASKAPCHECK(s1 == s2, "ids and arrayOfSpectrums are not the same. " 
                << " s1 = " << s1 << ", s2 = " << s2);

    if ( ! ids.empty() ) {
        // open the file again
        if ( fits_open_file(&itsFitsPtr,itsName.c_str(), READWRITE, &itsStatus) )
            PrintError(itsStatus);

        // move to the HDU 2 which is the spectrum binary table
        int hduType;
        if ( fits_movabs_hdu(itsFitsPtr,spectrumHDU(),&hduType,&itsStatus) )
            PrintError(itsStatus);

        long firstElem = 1;

        auto nrow = arrayOfSpectrums.nrow();
        for (auto r = 0; r < nrow; r++) {
            this->add(ids[r],arrayOfSpectrums.row(r));
        }
    }
}

void
FitsAuxImageSpectra::get(const long row, SpectrumT& spectrum)
{
    std::cout << "get spectrum"  << std::endl;
    itsStatus = 0;
    // open the file again
    if ( fits_open_file(&itsFitsPtr,itsName.c_str(), READWRITE, &itsStatus) )
         PrintError(itsStatus);

    // move to the HDU 2 which is the spectrum binary table
    int hduType;
    if ( fits_movabs_hdu(itsFitsPtr,spectrumHDU(),&hduType,&itsStatus) )
         PrintError(itsStatus);

    float fnull = 0.0;    
    int anynull = 0;
    const int spectrumCol = 2; // for now
    const long felem = 1;
    const long nelem = itsNChannels;
    
    spectrum.resize(itsNChannels);
    float* data = const_cast<float *> (spectrum.data());
    if (fits_read_col(itsFitsPtr,TFLOAT,spectrumCol,row,felem,nelem,&fnull,data,&anynull,&itsStatus))
        PrintError(itsStatus);

    // close the FITS file
    if (fits_close_file(itsFitsPtr,&itsStatus))
         PrintError(itsStatus);

    std::cout << "spectrum size: " << spectrum.size() << std::endl;
}

void
FitsAuxImageSpectra::get(const std::string& id, SpectrumT& spectrum)
{
    const auto& iter = itsId2RowMap.find(id);
    if ( iter != itsId2RowMap.end() ) {
        this->get(iter->second,spectrum);
    }
}
