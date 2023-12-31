/// @file 
/// @brief An implementation of the IMiscTableInfoHolder interface
/// @details The main idea of this class is to carry around 
/// additional information describing how the table is going to be
/// processed. The overall design is similar to ISubtableInfoHolder and
/// derived classes. Although this additional info can be assigned to 
/// either TableHolder or SubtableInfoHolder, making a separate tree of
/// classes seems to be a more structured approach.
/// Finally, having this miscelaneous information carried between classes 
/// the same way as the table itself and associated derived information,
/// allows to avoid creating multiple copies for data source and iterators
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


// own includes
#include <askap/dataaccess/MiscTableInfoHolder.h>
#include <askap/askap/AskapError.h>

using namespace askap;
using namespace askap::accessors;

/// @brief construct a holder of miscellaneous table processing information
/// @details The class just remembers the default column name passed in
/// this method.
/// @param[in] dataColumn the name of the data column used by default
MiscTableInfoHolder::MiscTableInfoHolder(const std::string &dataColumn) :
       itsDefaultDataColumnName(dataColumn) 
{
   ASKAPCHECK(dataColumn != "", "Attempted to set an empty string to the default data column");
}

/// @brief obtain the name of the data column to use by default
/// @details The code allows to read/write data not only from DATA
/// column of the measurement set, but from any other suitable column
/// as well. It is possible to change the name of the column via 
/// selector (ITableDataSelector or derived classes only, as this is
/// a table-specific operation. Hence a dynamic_cast may be required).
/// An alternative is to change the default column name via data source
/// constructor (it will be carried accross all required classes by 
/// TableManager).
/// @return a const reference to the default name of the data column 
const std::string & MiscTableInfoHolder::defaultDataColumnName() const
{
  return itsDefaultDataColumnName;
}
