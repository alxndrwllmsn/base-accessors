/// @file TableSelectionDiskCache.h
///
/// @brief A class to keep a disk cache of selected tables
/// @details Singleton helper class keep and delete temporary casacore tables.
/// We keep a disk cache of selected tables to avoid the cost of
/// repeated table selection and delete them at program end
/// @copyright (c) 2025 CSIRO
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
#ifndef ASKAP_ACCESSORS_TABLE_SELECTION_DISK_CACHE_H
#define ASKAP_ACCESSORS_TABLE_SELECTION_DISK_CACHE_H

#include <askap_accessors.h>

/// ASKAPsoft includes
#include "casacore/tables/TaQL/ExprNode.h"

namespace askap {
namespace accessors {

/// @brief A class to keep a disk cache of selected tables
/// @details Singleton helper class creates and deletes temporary casacore tables.
/// We keep a disk cache of selected tables to avoid the cost of
/// repeated table selection and delete them at program end
class TableSelectionDiskCache
{
  public:
    /// @brief create cached selection of table
    /// @details Apply the table selection and cache the result.
    /// Note that the exprNode keeps a reference to the table to use.
    /// @param[in] name, the name of the cache file
    /// @param[in] exprNode, the table selection to apply
    static void create(const std::string& name, const casacore::TableExprNode& exprNode);

    /// @brief select table or retrieve cached selection of table
    /// @details Apply the table selection or return the 
    /// cached result. Note that the exprNode keeps a reference to the table to use.
    /// @param[in] name, the name of the cache file
    /// @param[in] exprNode, the table selection to apply
    /// @returns the selected table
    static const casacore::Table table(const std::string name, const casacore::TableExprNode& exprNode);

    /// @brief clear the selection cache
    /// @details This will clear all currently cached selected tables. Normally this is done from 
    /// the destructor, but there may be reasons to do this explictly
    static void clear();

  private:
    // No instantiation
    TableSelectionDiskCache() {};

    // No copying. 
    TableSelectionDiskCache(const TableSelectionDiskCache&) = delete;
    TableSelectionDiskCache& operator= (const TableSelectionDiskCache&) = delete;

    // destructor deletes temporary tables at program end
    ~TableSelectionDiskCache() {
      TableSelectionDiskCache::clear();
    }

    // return the single instance of the class
    static TableSelectionDiskCache& getInstance();

    // the set of unique names for the cached selected tables
    std::set<std::string> itsNames;
};

} // end of namespace accessors

} // end of namespace askap
#endif // #ifndef ASKAP_ACCESSORS_TABLE_SELECTION_DISK_CACHE_H
