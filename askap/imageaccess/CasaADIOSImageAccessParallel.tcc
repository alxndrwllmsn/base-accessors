/// @file CasaADIOSImageAccessParallel.cc
/// @brief Access casa image via ADIOS Storage Manager using parallel interface
/// @details This class implements IImageAccess interface for CASA image and specifies use of the ADIOS storage manager
///
/// @copyright (c) 2007 CSIRO
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
/// @author Alexander Williamson <alex.williamson@icrar.org>
/// @author Pascal Jahan Elahi <pascal.elahi@pawsey.org.au>

// have to remove this due to templating #include <askap_accessors.h>

#include <askap/imageaccess/CasaADIOSImageAccessParallel.h>

#include <askap/askap/AskapLogging.h>
#include <askap/imageaccess/ADIOSImage.h>
#include <casacore/images/Images/SubImage.h>
#include <casacore/images/Regions/ImageRegion.h>
#include <casacore/images/Regions/RegionHandler.h>

#include <mpi.h>
#include <askap/askapparallel/MPIComms.h>
#include <askap/askapparallel/AskapParallel.h>

ASKAP_LOGGER(casaADIOSImAccessParallelLogger, ".casaADIOSImageAccessorParallel");

using namespace askap;
using namespace askap::accessors;

/// @brief constructor
/// @param[in] comms, MPI communicator
/// @param[in] config, configuration file name 
template <class T>
CasaADIOSImageAccessParallel<T>::CasaADIOSImageAccessParallel(askapparallel::AskapParallel &comms, std::string config, size_t comm_index):
    itsComms(comms), itsParallel(-1), itsCommIndex(comm_index)
{
    configname = config;
    if (config != "") {
        ASKAPLOG_INFO_STR(casaADIOSImAccessParallelLogger, "Creating parallel ADIOS accessor with configuration file " << config);
    }
    else {
        ASKAPLOG_INFO_STR(casaADIOSImAccessParallelLogger, "Creating parallel ADIOS accessor with default configuration");
    }
}

// reading methods

/// @brief obtain the shape
/// @param[in] name image name
/// @return full shape of the given image
template <class T>
casacore::IPosition CasaADIOSImageAccessParallel<T>::shape(const std::string &name) const
{
    ADIOSImage<T> img(itsComms, name);
    return img.shape();
}

/// @brief read full image
/// @param[in] name image name
/// @return array with pixels
template <class T>
casacore::Array<T> CasaADIOSImageAccessParallel<T>::read(const std::string &name) const
{
    ASKAPLOG_INFO_STR(casaADIOSImAccessParallelLogger, "Reading parallel CASA ADIOS image " << name);
    imagePtr_p.reset();
    imagePtr_p.reset(new ADIOSImage<T>(itsComms, name));
    if (imagePtr_p->hasPixelMask()) {
        ASKAPLOG_INFO_STR(casaADIOSImAccessParallelLogger, " - setting unmasked pixels to zero");
        // generate an Array of zeros and copy the elements for which the mask is true
        casacore::Array<T> tempArray(imagePtr_p->get().shape(), static_cast<T>(0.0));
        tempArray = casacore::MaskedArray<T>(imagePtr_p->get(), imagePtr_p->getMask(), casacore::True);
        return tempArray;
    } else {
        return imagePtr_p->get();
    }
}

/// @brief read part of the image
/// @param[in] name image name
/// @param[in] blc bottom left corner of the selection
/// @param[in] trc top right corner of the selection
/// @return array with pixels for the selection only
template <class T>
casacore::Array<T> CasaADIOSImageAccessParallel<T>::read(const std::string &name, const casacore::IPosition &blc,
        const casacore::IPosition &trc) const
{
    ASKAPLOG_INFO_STR(casaADIOSImAccessParallelLogger, "Reading a slice of the parallel CASA ADIOS image " << name << " from " << blc << " to " << trc);
    imagePtr_p.reset();
    imagePtr_p.reset(new ADIOSImage<T>(itsComms, name));

    if (imagePtr_p->hasPixelMask()) {
        ASKAPLOG_INFO_STR(casaADIOSImAccessParallelLogger, " - setting unmasked pixels to zero");
        // generate an Array of zeros and copy the elements for which the mask is true
        const casacore::Slicer slicer(blc, trc, casacore::Slicer::endIsLast);
        casacore::Array<T> tempSlice(imagePtr_p->getSlice(slicer).shape(), static_cast<T>(0.0));
        tempSlice = casacore::MaskedArray<T>(imagePtr_p->getSlice(slicer), imagePtr_p->getMaskSlice(slicer), casacore::True);
        return tempSlice;
    } else {
        return imagePtr_p->getSlice(casacore::Slicer(blc, trc, casacore::Slicer::endIsLast));
    }
}

/// @brief Determine whether an image has a mask
/// @param[in] nam image name
/// @return True if image has a mask, False if not.
template <class T>
bool CasaADIOSImageAccessParallel<T>::isMasked(const std::string &name) const
{
    ADIOSImage<T> img(itsComms, name);
    return img.hasPixelMask();
}

/// @brief read the mask for the full image
/// @param[in] name image name
/// @return bool array with mask values - 1=good, 0=bad
template <class T>
casacore::LogicalArray CasaADIOSImageAccessParallel<T>::readMask(const std::string &name) const
{

    ADIOSImage<T> img(itsComms, name);
    if (img.hasPixelMask()) {
        return img.getMask();
    } else {
        casacore::LogicalArray mask(img.shape(),true);
        return mask;
    }

}

/// @brief read the mask for part of the image
/// @param[in] name image name
/// @param[in] blc bottom left corner of the selection
/// @param[in] trc top right corner of the selection
/// @return bool array with mask values - 1=good, 0=bad
template <class T>
casacore::LogicalArray CasaADIOSImageAccessParallel<T>::readMask(const std::string &name, const casacore::IPosition &blc,
                                                    const casacore::IPosition &trc) const
{
    ADIOSImage<T> img(itsComms, name);
    const casacore::Slicer slicer(blc, trc, casacore::Slicer::endIsLast);
    if (img.hasPixelMask()) {
        return img.getMaskSlice(slicer);
    } else {
        casacore::LogicalArray mask(slicer.length(),true);
        return mask;
    }

}

/// @brief obtain coordinate system info
/// @param[in] name image name
/// @return coordinate system object
template <class T>
casacore::CoordinateSystem CasaADIOSImageAccessParallel<T>::coordSys(const std::string &name) const
{
    ADIOSImage<T> img(itsComms, name);
    return img.coordinates();
}

template <class T>
casacore::CoordinateSystem CasaADIOSImageAccessParallel<T>::coordSysSlice(const std::string &name, const casacore::IPosition &blc,
        const casacore::IPosition &trc) const
{
    casacore::Slicer slc(blc, trc, casacore::Slicer::endIsLast);
    ASKAPLOG_INFO_STR(casaADIOSImAccessParallelLogger, " CasaADIOSImageAccessParallel - Slicer " << slc);
    casacore::PagedImage<T> img(name);
    casacore::SubImage<T> si = casacore::SubImage<T>(img, slc, casacore::AxesSpecifier(casacore::True));
    return si.coordinates();


}
/// @brief obtain beam info
/// @param[in] name image name
/// @return beam info vector
template <class T>
casacore::Vector<casacore::Quantum<double> > CasaADIOSImageAccessParallel<T>::beamInfo(const std::string &name) const
{
    ADIOSImage<T> img(itsComms, name);
    casacore::ImageInfo ii = img.imageInfo();
    if (!img.imageInfo().hasMultipleBeams()) {
      return ii.restoringBeam().toVector();
    } else {
      return casacore::Vector<casacore::Quantum<double>>();
    }
}

/// @brief get restoring beam info
/// @param[in] name image name
/// @return beamlist  list of beams, beamlist will be empty if image only has a single beam
template <class T>
BeamList CasaADIOSImageAccessParallel<T>::beamList(const std::string &name) const
{
    ADIOSImage<T> img(itsComms, name);
    casacore::ImageInfo ii = img.imageInfo();
    BeamList bl;
    if (img.imageInfo().hasMultipleBeams()) {
      for (int chan = 0; chan < ii.nChannels(); chan++) {
        casacore::GaussianBeam gb = ii.restoringBeam(chan,0);
        bl[chan] = gb.toVector();
      }
    }
    return bl;
}

template <class T>
std::string CasaADIOSImageAccessParallel<T>::getUnits(const std::string &name) const
{
    casacore::Table tmpTable(name);
    std::string units = tmpTable.keywordSet().asString("units");
    return units;
}

/// @brief Get a particular keyword from the image metadata (A.K.A header)
/// @details This reads a given keyword to the image metadata.
/// @param[in] name Image name
/// @param[in] keyword The name of the metadata keyword
/// @return pair of strings - keyword value and comment
template <class T>
std::pair<std::string, std::string> CasaADIOSImageAccessParallel<T>::getMetadataKeyword(const std::string &name, const std::string &keyword) const
{

    ADIOSImage<T> img(itsComms, name);
    casacore::TableRecord miscinfo = img.miscInfo();
    std::string value = "";
    std::string comment = "";
    if (miscinfo.isDefined(keyword)) {
        value = miscinfo.asString(keyword);
        comment = miscinfo.comment(keyword);
    } else {
        ASKAPLOG_DEBUG_STR(casaADIOSImAccessParallelLogger, "Keyword " << keyword << " is not defined in metadata for image " << name);
    }
    return std::pair<std::string,std::string>(value,comment);

}


// writing methods

/// @brief create a new image
/// @details A call to this method should preceed any write calls. The actual
/// image may be created only upon the first write call. Details depend on the
/// implementation.
/// @param[in] name image name
/// @param[in] shape full shape of the image
/// @param[in] csys coordinate system of the full image
template <class T>
void CasaADIOSImageAccessParallel<T>::create(const std::string &name, const casacore::IPosition &shape,
                             const casacore::CoordinateSystem &csys)
{
    ASKAPLOG_INFO_STR(casaADIOSImAccessParallelLogger, "Creating a new CASA ADIOS image " << name << " with the shape " << shape);
    imagePtr_p.reset(new ADIOSImage<T>(itsComms, itsCommIndex, 
        casacore::TiledShape(shape), csys, name, configname));
}

/// @brief write full image
/// @param[in] name image name
/// @param[in] arr array with pixels
template <class T>
void CasaADIOSImageAccessParallel<T>::write(const std::string &name, const casacore::Array<T> &arr)
{
    ASKAPLOG_INFO_STR(casaADIOSImAccessParallelLogger, "Writing an array with the shape " << arr.shape() << " into a CASA ADIOS image " << name);
    imagePtr_p->put(arr);
}

/// @brief write a slice of an image
/// @param[in] name image name
/// @param[in] arr array with pixels
/// @param[in] where bottom left corner where to put the slice to (trc is deduced from the array shape)
template <class T>
void CasaADIOSImageAccessParallel<T>::write(const std::string &name, const casacore::Array<T> &arr,
                               const casacore::IPosition &where)
{
    ASKAPLOG_INFO_STR(casaADIOSImAccessParallelLogger, "Writing a slice with the shape " << arr.shape() << " into a CASA ADIOS image " <<
                      name << " at " << where);
    imagePtr_p->putSlice(arr, where);
}

/// @brief write a slice of an image and mask
/// @param[in] name image name
/// @param[in] arr array with pixels
/// @param[in] mask array with mask
template <class T>
void CasaADIOSImageAccessParallel<T>::write(const std::string &name, const casacore::Array<T> &arr,
                               const casacore::Array<bool> &mask)
{
    ASKAPLOG_INFO_STR(casaADIOSImAccessParallelLogger, "Writing image & mask with the shape " << arr.shape() << " into a CASA ADIOS image " <<
                      name);
    imagePtr_p->put(arr);
    imagePtr_p->pixelMask().put(mask);
}


/// @brief write a slice of an image and mask
/// @param[in] name image name
/// @param[in] arr array with pixels
/// @param[in] mask array with mask
/// @param[in] where bottom left corner where to put the slice to (trc is deduced from the array shape)
template <class T>
void CasaADIOSImageAccessParallel<T>::write(const std::string &name, const casacore::Array<T> &arr,
                               const casacore::Array<bool> &mask, const casacore::IPosition &where)
{
    ASKAPLOG_INFO_STR(casaADIOSImAccessParallelLogger, "Writing a slice with the shape " << arr.shape() << " into a CASA ADIOS image " <<
                      name << " at " << where);
    imagePtr_p->putSlice(arr, where);
    imagePtr_p->pixelMask().putSlice(mask, where);
}


/// @brief write a slice of an image mask
/// @param[in] name image name
/// @param[in] arr array with pixels
/// @param[in] where bottom left corner where to put the slice to (trc is deduced from the array shape)
template <class T>
void CasaADIOSImageAccessParallel<T>::writeMask(const std::string &name, const casacore::Array<bool> &mask,
                                const casacore::IPosition &where)
{
    ASKAPLOG_INFO_STR(casaADIOSImAccessParallelLogger, "Writing a slice with the shape " << mask.shape() << " into a CASA ADIOS image " <<
                      name << " at " << where);
    imagePtr_p->pixelMask().putSlice(mask, where);
}

/// @brief write a slice of an image mask
/// @param[in] name image name
/// @param[in] arr array with pixels
template <class T>
void CasaADIOSImageAccessParallel<T>::writeMask(const std::string &name, const casacore::Array<bool> &mask)
{
    ASKAPLOG_INFO_STR(casaADIOSImAccessParallelLogger, "Writing a full mask with the shape " << mask.shape() << " into a CASA ADIOS image " <<
                      name);
    imagePtr_p->pixelMask().put(mask);
}
/// @brief set brightness units of the image
/// @details
/// @param[in] name image name
/// @param[in] units string describing brightness units of the image (e.g. "Jy/beam")
template <class T>
void CasaADIOSImageAccessParallel<T>::setUnits(const std::string &name, const std::string &units)
{
    imagePtr_p->setUnits(casacore::Unit(units));
}

/// @brief set restoring beam info
/// @details For the restored image we want to carry size and orientation of the restoring beam
/// with the image. This method allows to assign this info.
/// @param[in] name image name
/// @param[in] maj major axis in radians
/// @param[in] min minor axis in radians
/// @param[in] pa position angle in radians
template <class T>
void CasaADIOSImageAccessParallel<T>::setBeamInfo(const std::string &name, double maj, double min, double pa)
{

    casacore::ImageInfo ii = imagePtr_p->imageInfo();
    ii.setRestoringBeam(casacore::Quantity(maj, "rad"), casacore::Quantity(min, "rad"), casacore::Quantity(pa, "rad"));
    imagePtr_p->setImageInfo(ii);
}

/// @brief set restoring beam info
/// @details For the restored image we want to carry size and orientation of the restoring beam
/// with the image. This method allows to assign this info for all channels
/// @param[in] name image name
/// @param[in] beamlist  list of beams
template <class T>
void CasaADIOSImageAccessParallel<T>::setBeamInfo(const std::string &name, const BeamList & beamlist)
{

    casacore::ImageInfo ii = imagePtr_p->imageInfo();
    ii.setAllBeams(beamlist.size(),1,casacore::GaussianBeam());
    for (const auto& beam : beamlist) {
      ASKAPDEBUGASSERT(beam.second.size()==3);
      ii.setBeam(beam.first,0,beam.second[0],beam.second[1],beam.second[2]);
    }
    imagePtr_p->setImageInfo(ii);
}

/// @brief apply mask to image
/// @details Deteails depend upon the implemenation - CASA images will have the pixel mask assigned
/// but FITS images will have it applied to the pixels ... which is an irreversible process
/// @param[in] name image name
/// @param[in] the mask
template <class T>
void CasaADIOSImageAccessParallel<T>::makeDefaultMask(const std::string &name)
{
    // Create a mask and make it default region.
    // need to assert sizes etc ...
    imagePtr_p->makeMask("mask", casacore::True, casacore::True);
    casacore::Array<casacore::Bool> mask(imagePtr_p->shape(),casacore::True);
    imagePtr_p->pixelMask().put(mask);

}


/// @brief Set a particular keyword for the metadata (A.K.A header)
/// @details This adds a given keyword to the image metadata.
/// @param[in] name Image name
/// @param[in] keyword The name of the metadata keyword
/// @param[in] value The value for the keyword, in string format
/// @param[in] desc A description of the keyword
template <class T>
void CasaADIOSImageAccessParallel<T>::setMetadataKeyword(const std::string &name, const std::string &keyword,
        const std::string value, const std::string &desc)
{
    casacore::TableRecord miscinfo = imagePtr_p->miscInfo();
    miscinfo.define(keyword, value);
    miscinfo.setComment(keyword, desc);
    imagePtr_p->setMiscInfo(miscinfo);

}

template <class T>
void CasaADIOSImageAccessParallel<T>::setMetadataKeywords(const std::string &name, const LOFAR::ParameterSet &keywords)
{
    casacore::TableRecord miscinfo = imagePtr_p->miscInfo();
    // Note: we could sort through the keywords here and pick out ones that need to go in places
    // other than miscInfo to be more compatible with casacore
    for (auto &elem : keywords) {
      const string keyword = elem.first;
      const std::vector<string> valanddesc = elem.second.getStringVector();
      if (valanddesc.size() > 0) {
        const string value = valanddesc[0];
        const string desc = (valanddesc.size() > 1 ? valanddesc[1] : "");

        const string type = (valanddesc.size() > 2 ? toUpper(valanddesc[2]) : "STRING");
        if (type == "INT") {
          try {
            const int intVal = std::stoi(value);
            miscinfo.define(keyword, intVal);
            miscinfo.setComment(keyword, desc);
          } catch (const std::invalid_argument&) {
            ASKAPLOG_WARN_STR(casaADIOSImAccessParallelLogger, "Invalid int value for header keyword "<<keyword<<" : "<<value);
          } catch (const std::out_of_range&) {
            ASKAPLOG_WARN_STR(casaADIOSImAccessParallelLogger, "Out of range int value for header keyword "<<keyword<<" : "<<value);
          }
        } else if (type == "DOUBLE") {
          try {
            const double doubleVal = std::stod(value);
            miscinfo.define(keyword, doubleVal);
            miscinfo.setComment(keyword, desc);
          } catch (const std::invalid_argument&) {
            ASKAPLOG_WARN_STR(casaADIOSImAccessParallelLogger, "Invalid double value for header keyword "<<keyword<<" : "<<value);
          } catch (const std::out_of_range&) {
            ASKAPLOG_WARN_STR(casaADIOSImAccessParallelLogger, "Out of range double value for header keyword "<<keyword<<" : "<<value);
          }
        } else if (type == "STRING") {
          miscinfo.define(keyword, value);
          miscinfo.setComment(keyword, desc);
        } else {
          ASKAPLOG_WARN_STR(casaADIOSImAccessParallelLogger, "Invalid type for header keyword "<<keyword<<" : "<<type);
        }
      }
    }
    imagePtr_p->setMiscInfo(miscinfo);
}

/// @brief Add HISTORY messages to the image metadata
/// @details Adds a list of strings detailing the history of the image
/// @param[in] name Image name
/// @param[in] historyLines History comments to add
template <class T>
void CasaADIOSImageAccessParallel<T>::addHistory(const std::string &name, const std::vector<std::string> &historyLines)
{
    casacore::LogIO log = imagePtr_p->logSink();
    for (const auto& history : historyLines) {
        log << history << casacore::LogIO::POST;
    }
}


/// @brief set info for image that can vary by e.g., channel
/// @details Add arbitrary info to the image as either keywords or a table
/// @param[in] name image name
/// @param[in] info record with information
template <class T>
void CasaADIOSImageAccessParallel<T>::setInfo(const std::string &name, const casacore::RecordInterface & info)
{
    // make a copy of the table record
    casacore::TableRecord updateTableRecord = imagePtr_p->miscInfo();
    // find the name of the info table.  this is the name field of the info sub record
    std::string infoTableName = "notfound";
    casacore::uInt subRecordFieldId = 0;
    casacore::uInt nFields = info.nfields();
    for(casacore::uInt f = 0; f < nFields; f++) {
        std::string name = info.name(f);
        casacore::DataType type = info.dataType(f);
        if ( type == casacore::DataType::TpRecord ) {
           infoTableName = name;
            break;
        }
    }

    // add the info to it
    updateTableRecord.defineRecord(infoTableName,info);
    // now set the updated record back to the image
    imagePtr_p->setMiscInfo(updateTableRecord);
}

/// @brief this methods retrieves the table(s) in the image and stores them in the casacore::Record
/// @param[in] name - image name
/// @param[in] tblName - name of the table to retrieve the data. if tblName = "All" then retrieve all
///                      the tables in the image
/// @param[out] info - casacore::Record to contain the tables' data.
template <class T>
void CasaADIOSImageAccessParallel<T>::getInfo(const std::string &name, const std::string& tableName, casacore::Record &info)
{
    ADIOSImage<T> img(itsComms, name);
    //casacore::TableRecord tableRecord = imagePtr_p->miscInfo().toRecord();
    casacore::Record tableRecord = img.miscInfo().toRecord();
    casacore::uInt nFields = tableRecord.nfields();
    for(casacore::uInt f = 0; f < nFields; f++) {
        std::string name = tableRecord.name(f);
        casacore::DataType type = tableRecord.dataType(f);
        // we know tables are stored as sub records in the image
        if ( type == casacore::DataType::TpRecord ) {
            if ( (tableName != name) && (tableName != "All") ) {
                // tableName is not the same as current sub record's name
                // AND table is not "All", process the next field
                continue;
            }
            // add the sub record which is a table data in this case to
            // the info object
            info.defineRecord(name,tableRecord.asRecord(f));
        }
    }
}
