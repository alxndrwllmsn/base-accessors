/// @file CasaADIOSImageAccessParallel.h
/// @brief Access casa image via ADIOS Storage Manager using parallel io
/// @details This class implements IImageAccess interface for CASA image and specifies use of the ADIOS storage manager
///
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

#ifndef ASKAP_ACCESSORS_CASA_ADIOS_IMAGE_ACCESS_PARALLEL_H
#define ASKAP_ACCESSORS_CASA_ADIOS_IMAGE_ACCESS_PARALLEL_H

#include <boost/shared_ptr.hpp>

#include <askap/imageaccess/IImageAccess.h>
#include <askap/imageaccess/ADIOSImage.h>
#include <askap/imageaccess/CasaADIOSImageAccess.h>
#include <mpi.h>
#include <askap/askapparallel/MPIComms.h>
#include <askap/askapparallel/AskapParallel.h>

//PJE - still a work in progress 

namespace askap {
namespace accessors {

/// @brief Access casa image via ADIOS2
/// @details This class implements IImageAccess interface for CASA image via ADIOS2
/// @ingroup imageaccess

template <class T = casacore::Float>
struct CasaADIOSImageAccessParallel : public CasaADIOSImageAccess<T> {

public:
    /// @brief constructor
    /// @param[in] comms, MPI communicator
    /// @param[in] config, string pointing to adios configuration
    CasaADIOSImageAccessParallel(askapparallel::AskapParallel &comms, std::string config = "", size_t comm_index = 0);

    // PJE -- uncertain how many functions would really require overwriting so as to 
    // do something with the comms. 

    //////////////////
    // Reading methods
    //////////////////

    /// @brief obtain the shape
    /// @param[in] name image name
    /// @return full shape of the given image
    virtual casacore::IPosition shape(const std::string &name) const override;

    /// @brief read full image
    /// @param[in] name image name
    /// @return array with pixels
    virtual casacore::Array<T> read(const std::string &name) const override;

    /// @brief read part of the image
    /// @param[in] name image name
    /// @param[in] blc bottom left corner of the selection
    /// @param[in] trc top right corner of the selection
    /// @return array with pixels for the selection only
    virtual casacore::Array<T> read(const std::string &name, const casacore::IPosition &blc,
                                    const casacore::IPosition &trc) const override;

    /// @brief Determine whether an image has a mask
    /// @param[in] nam image name
    /// @return True if image has a mask, False if not.
    virtual bool isMasked(const std::string &name) const override;

    /// @brief read the mask for the full image
    /// @param[in] name image name
    /// @return bool array with mask values - 1=good, 0=bad
    virtual casacore::LogicalArray readMask(const std::string &name) const override;

    /// @brief read the mask for part of the image
    /// @param[in] name image name
    /// @param[in] blc bottom left corner of the selection
    /// @param[in] trc top right corner of the selection
    /// @return bool array with mask values - 1=good, 0=bad
    virtual casacore::LogicalArray readMask(const std::string &name, const casacore::IPosition &blc,
                                    const casacore::IPosition &trc) const override;

    /// @brief obtain coordinate system info
    /// @param[in] name image name
    /// @return coordinate system object
    virtual casacore::CoordinateSystem coordSys(const std::string &name) const override;

    /// @brief obtain coordinate system info for part of an image
    /// @param[in] name image name
    /// @return coordinate system object
    virtual casacore::CoordinateSystem coordSysSlice(const std::string &name, const casacore::IPosition &blc,
            const casacore::IPosition &trc) const override;

    /// @brief obtain beam info
    /// @param[in] name image name
    /// @return beam info vector
    virtual casacore::Vector<casacore::Quantum<double> > beamInfo(const std::string &name) const override;

    /// @brief obtain beam info
    /// @param[in] name image name
    /// @return beam info list
    virtual BeamList beamList(const std::string &name) const override;

    /// @brief obtain pixel units
    /// @param[in] name image name
    /// @return units string
    virtual std::string getUnits(const std::string &name) const override;

    /// @brief this methods retrieves the table(s) in the image and stores them in the casacore::Record
    /// @param[in] name - image name
    /// @param[in] tblName - name of the table to retrieve the data. if tblName = "All" then retrieve all
    ///                      the tables in the image
    /// @param[out] info - casacore::Record to contain the tables' data.
    virtual void getInfo(const std::string &name, const std::string& tableName, casacore::Record &info) override;

    /// @brief Get a particular keyword from the image metadata (A.K.A header)
    /// @details This reads a given keyword from the image metadata.
    /// @param[in] name Image name
    /// @param[in] keyword The name of the metadata keyword
    /// @return pair of strings - keyword value and comment
    virtual std::pair<std::string, std::string> getMetadataKeyword(const std::string &name, const std::string &keyword) const override;

    //////////////////
    // Writing methods
    //////////////////

    /// @brief create a new image
    /// @details A call to this method should preceed any write calls. The actual
    /// image may be created only upon the first write call. Details depend on the
    /// implementation.
    /// @param[in] name image name
    /// @param[in] shape full shape of the image
    /// @param[in] csys coordinate system of the full image
    virtual void create(const std::string &name, const casacore::IPosition &shape,
                        const casacore::CoordinateSystem &csys);

    /// @brief write full image
    /// @param[in] name image name
    /// @param[in] arr array with pixels
    virtual void write(const std::string &name, const casacore::Array<T> &arr) override;

    /// @brief write a slice of an image
    /// @param[in] name image name
    /// @param[in] arr array with pixels
    /// @param[in] where bottom left corner where to put the slice to (trc is deduced from the array shape)
    virtual void write(const std::string &name, const casacore::Array<T> &arr,
                       const casacore::IPosition &where) override;

    /// @brief write an image and mask
    /// @param[in] name image name
    /// @param[in] arr array with pixels
    /// @param[in] mask array with mask
    virtual void write(const std::string &name, const casacore::Array<T> &arr,
                       const casacore::Array<bool> &mask);


    /// @brief write a slice of an image and mask
    /// @param[in] name image name
    /// @param[in] arr array with pixels
    /// @param[in] mask array with mask
    /// @param[in] where bottom left corner where to put the slice to (trc is deduced from the array shape)
    virtual void write(const std::string &name, const casacore::Array<T> &arr,
                       const casacore::Array<bool> &mask, const casacore::IPosition &where) override;

    /// @brief write a slice of an image mask
    /// @param[in] name image name
    /// @param[in] arr array with pixels
    /// @param[in] where bottom left corner where to put the slice to (trc is deduced from the array shape)
    virtual void writeMask(const std::string &name, const casacore::Array<bool> &mask,
                           const casacore::IPosition &where) override;
    /// @brief write a slice of an image mask
    /// @param[in] name image name
    /// @param[in] arr array with pixels

    virtual void writeMask(const std::string &name, const casacore::Array<bool> &mask) override;


    /// @brief set brightness units of the image
    /// @details
    /// @param[in] name image name
    /// @param[in] units string describing brightness units of the image (e.g. "Jy/beam")
    virtual void setUnits(const std::string &name, const std::string &units) override;

    /// @brief set restoring beam info
    /// @details For the restored image we want to carry size and orientation of the restoring beam
    /// with the image. This method allows to assign this info.
    /// @param[in] name image name
    /// @param[in] maj major axis in radians
    /// @param[in] min minor axis in radians
    /// @param[in] pa position angle in radians
    virtual void setBeamInfo(const std::string &name, double maj, double min, double pa) override;

    /// @brief set restoring beam info for all channels
    /// @details For the restored image we want to carry size and orientation of the restoring beam
    /// with the image. This method allows to assign this info.
    /// @param[in] name image name
    /// @param[in] beamlist The list of beams
    virtual void setBeamInfo(const std::string &name, const BeamList & beamlist) override;

    /// @brief apply mask to image
    /// @details Deteails depend upon the implemenation - CASA images will have the pixel mask assigned
    /// but FITS images will have it applied to the pixels ... which is an irreversible process
    /// @param[in] name image name
    /// @param[in] the mask

    virtual void makeDefaultMask(const std::string &name) override;

    /// @brief Set a particular keyword for the metadata (A.K.A header)
    /// @details This adds a given keyword to the image metadata.
    /// @param[in] name Image name
    /// @param[in] keyword The name of the metadata keyword
    /// @param[in] value The value for the keyword, in string format
    /// @param[in] desc A description of the keyword
    virtual void setMetadataKeyword(const std::string &name, const std::string &keyword,
                                    const std::string value, const std::string &desc = "") override;

    /// @brief Set the keywords for the metadata (A.K.A header)
    /// @details This adds keywords to the image metadata.
    /// @param[in] name Image name
    /// @param[in] keywords A parset with keyword entries (KEYWORD = ["keyword value","keyword description","STRING"])
    virtual void setMetadataKeywords(const std::string &name, const LOFAR::ParameterSet &keywords) override;

    /// @brief Add HISTORY messages to the image metadata
    /// @details Adds a list of strings detailing the history of the image
    /// @param[in] name Image name
    /// @param[in] historyLines History comments to add
    virtual void addHistory(const std::string &name, const std::vector<std::string> &historyLines) override;

    /// @brief set info for image that can vary by e.g., channel
    /// @details Add arbitrary info to the image as either keywords or a binary table
    /// @param[in] name image name
    /// @param[in] info record with information
    virtual void setInfo(const std::string &name, const casacore::RecordInterface & info) override;

private:
    // PJE -- added some private functions to assist parallel IO 
    /// @brief check if we can do parallel I/O on the file
    bool canDoParallelIO(const std::string &name) const;

    // add the parallel comms and the name associated with it, and whether parallel can be done
    int itsParallel;
    askapparallel::AskapParallel& itsComms;
    std::string itsName;
    mutable std::unique_ptr<ADIOSImage<T>> imagePtr_p;
    std::string configname;
    size_t itsCommIndex;
      
};

} // namespace accessors
} // namespace askap

#include "CasaADIOSImageAccessParallel.tcc"
#endif
