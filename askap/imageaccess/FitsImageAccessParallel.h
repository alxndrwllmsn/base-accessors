/// @file FitsImageAccess.h
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
/// @author Mark Wieringa <mark.wieringa@csiro.au>
///

#ifndef ASKAP_ACCESSORS_FITS_IMAGE_ACCESS_PARALLEL_H
#define ASKAP_ACCESSORS_FITS_IMAGE_ACCESS_PARALLEL_H

#include <boost/shared_ptr.hpp>

#include <askap/askapparallel/AskapParallel.h>
#include <askap/imageaccess/FitsImageAccess.h>

namespace askap {
namespace accessors {

/// @brief Access casa image using collective MPI I/O
/// @details This class adds collective read/write functions for FITS images
/// At present it can deal with 2D images and 3D cubes. Both images and cubes can
/// have up to 4 dimensions, but only 3 can differ from 1. For cubes one of the last
/// last two dimensions must be 1.
/// Distribution over the 1st axis (iax=0) may be less efficient than over the 2nd or 3rd.
/// @ingroup imageaccess

// TODO: integrate this with the base classes or make it a standalone class

struct FitsImageAccessParallel : public FitsImageAccess {

    public:

        /// @brief read part of the image - collective MPI read
        /// @param[in] comms, MPI communicator
        /// @param[in] name image name
        /// @param[in] iax, axis to distribute over: use 1 for 2D images,  0, 1 or 2 for x, y, z, i.e., yz planes, xz planes, xy planes
        /// @return array with pixels for the section of the image read
        casacore::Array<float> read_all(askapparallel::AskapParallel& comms,
                                        const std::string &name, int iax) const;


        /// @brief write an image - collective MPI write.
        /// @details Note that the fits header must be written to disk before calling this.
        /// @param[in] comms, MPI communicator
        /// @param[in] name image name
        /// @param[in] arr array with pixels. Array dimension iax * #ranks must match the corresponding image dimension
        /// @param[in] iax, axis to distribute over: 0, 1 or 2 for x, y, z, i.e., yz planes, xz planes, xy planes
        void write_all(askapparallel::AskapParallel& comms,
                        const std::string &name, const casacore::Array<float> &arr, int iax) const;

        /// @brief copy the header of a fits image (i.e., copies the fits 'cards' preceeding the data)
        /// @param[in] infile, the input fits file
        /// @param[in] outfile, the output fits file (overwritten if it exists)
        void copy_header(const casa::String &infile, const casa::String& outfile) const;

    private:
        /// @brief determine image dimensions (up to 3 non degenerate axes) and headersize from file
        void decode_header(const casa::String& infile, casa::Int& nx, casa::Int& ny,
                            casa::Int& nz, casa::Long& headersize) const;

        /// @brief add padding to the fits file to make it complient
        void fits_padding(const casa::String& filename) const;

        /// @brief determine the file access pattern, offset to start reading or writing and the buffer shape needed
        void setFileAccess(askapparallel::AskapParallel& comms, const casa::String& fullname, casa::IPosition& bufshape,
            MPI_Offset& offset, MPI_Datatype&  filetype, casacore::Int iax) const;

};


} // namespace accessors
} // namespace askap

#endif