/// @file
/// @brief A handler of  ANTENNA subtable
/// @details This class provides access to the ANTENNA subtable (which contains 
/// antenna mounts and positions for all antennas). It caches the whole table
/// in constructor and then returns cached values. 
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

#ifndef ASKAP_ACCESSORS_MEM_ANTENNA_SUBTABLE_HANLDER_H
#define ASKAP_ACCESSORS_MEM_ANTENNA_SUBTABLE_HANLDER_H

// casa includes
#include <casacore/tables/Tables/Table.h>
#include <casacore/casa/Arrays/Vector.h>

// own includes
#include <askap/dataaccess/IAntennaSubtableHandler.h>

namespace askap {

namespace accessors {

/// @brief A handler of  ANTENNA subtable
/// @details This class provides access to the ANTENNA subtable (which contains 
/// antenna mounts and positions for all antennas). It caches the whole table
/// in constructor and then returns cached values. 
/// @ingroup dataaccess_tab
struct MemAntennaSubtableHandler : virtual public IAntennaSubtableHandler {
  
  /// read all required information from the ANTENNA subtable
  /// @param[in] ms an input measurement set (a table which has an
  /// ANTENNA subtable)
  explicit MemAntennaSubtableHandler(const casacore::Table &ms);
  
  /// @brief obtain the position of the given antenna
  /// @details
  /// @param[in] antID antenna ID to return the position for
  /// @return a reference to the MPosition measure
  virtual const casacore::MPosition& getPosition(casacore::uInt antID) const;
  
  /// @brief obtain the mount type for the given antenna
  /// @details
  /// @param[in] antID antenna ID to return the position for
  /// @return a string describing the mount type
  virtual const casacore::String& getMount(casacore::uInt antID) const;
  
  /// @brief check whether all antennas are equatorialy mounted
  /// @details
  /// This method checks the mount type for all antennas to be 
  /// either EQUATORIAL or equatorial. This mount type doesn't require
  /// parallactic angle rotation and can be trated separately.
  /// @return true, if all antennas are equatorially mounted
  virtual bool allEquatorial() const throw();   
  
  /// @brief get the number of antennas
  /// @details
  /// This method returns the number of antennas (i.e. all antID indices
  /// are expected to be less than this number). Following the general
  /// assumptions about ANTENNA subtable, this number is assumed to be
  /// fixed.
  /// @return total number of antennas 
  virtual casacore::uInt getNumberOfAntennas() const;
    
private:
  /// a cache of antenna mounts
  casacore::Vector<casacore::String> itsMounts;
  
  /// a cache of antenna positions
  casacore::Vector<casacore::MPosition> itsPositions;
  
  /// an internal buffer for a flag showing whether all antennas are equatorially
  /// mounted
  bool itsAllEquatorial;   
};


} // namespace accessors

} // namespace askap

#endif // #ifndef MEM_ANTENNA_SUBTABLE_HANLDER_H
