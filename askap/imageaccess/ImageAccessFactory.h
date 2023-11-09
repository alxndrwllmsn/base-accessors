/// @file
/// @brief Build an appropriate image access class
/// @details This file contains a factory method generating a shared pointer to the image
/// accessor from the parset file
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
/// @author Max Voronkov <maxim.voronkov@csiro.au>
///

#ifndef ASKAP_ACCESSORS_IMAGE_ACCESS_FACTORY_H
#define ASKAP_ACCESSORS_IMAGE_ACCESS_FACTORY_H

#include <Common/ParameterSet.h>
#include <boost/shared_ptr.hpp>

#include <askap/imageaccess/IImageAccess.h>
#include <askap/askapparallel/AskapParallel.h>

namespace askap {
namespace accessors {

/// @brief Build an appropriate image access class
/// @details This is a factory method generating a shared pointer to the image
/// accessor from the parset file
/// @param[in] parset parameters containing description of image accessor to be constructed
/// @return shared pointer to the image access object
/// @note CASA images are used by default
boost::shared_ptr<IImageAccess< casacore::Float > > imageAccessFactory(const LOFAR::ParameterSet &parset);

/// @brief Build an appropriate image access class
/// @details This is a factory method generating a shared pointer to the image
/// accessor from the parset file
/// @param[in] parset parameters containing description of image accessor to be constructed
/// @param[in] comms, MPI communicator, needed for parallel I/O
/// @return shared pointer to the image access object
/// @note CASA images are used by default
boost::shared_ptr<IImageAccess< casacore::Float > > imageAccessFactory(const LOFAR::ParameterSet &parset,
                                                                       askapparallel::AskapParallel &comms,
                                                                       size_t comm_index = 0);

} // namespace accessors
} // namespace askap

#endif
