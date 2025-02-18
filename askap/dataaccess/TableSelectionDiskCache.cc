/// @file TableSelectionDiskCache.cc
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

// casa includes
#include <casacore/casa/OS/File.h>

/// ASKAPsoft includes
#include <askap/askap/AskapLogging.h>
#include "casacore/tables/TaQL/ExprNode.h"
ASKAP_LOGGER(logger, "dataaccess");
/// Local package
#include <askap/dataaccess/TableSelectionDiskCache.h>

using namespace casa;
using namespace askap;
using namespace askap::accessors;

/// @brief create cached selection of table
/// @details Apply the table selection and cache the result.
///  Note that the exprNode keeps a reference to the table to use.
/// @param[in] exprNode, the table selection to apply
void TableSelectionDiskCache::create(const std::string& name, const casacore::TableExprNode& exprNode) {
  Table table = exprNode.table();
  // do the selection, this may be a costly operation
  casacore::Table sel = table(exprNode);
  // save the selection to disk
  ASKAPLOG_DEBUG_STR(logger,"Creating selected table: "<<name);
  sel.rename(name,casacore::Table::New);
  sel.flush();
  // make sure table gets deleted eventually
  getInstance().itsNames.insert(name);
}

/// @brief select table or retrieve cached selection of table
/// @details Apply the table selection or return the cached selection.
///  Note that the exprNode keeps a reference to the table to use.
/// @param[in] exprNode, the table selection to apply
/// @returns the selected table
const casacore::Table TableSelectionDiskCache::table(const std::string name, const casacore::TableExprNode& exprNode) {
  if (casacore::File(name).exists()) {
    ASKAPLOG_DEBUG_STR(logger,"Found pre-existing selected table: "<<name);
    casacore::Table sel(name);
    return sel;
  } else {
    Table table = exprNode.table();
    // do the selection, this may be a costly operation
    casacore::Table sel = table(exprNode);
    return sel;
  }
}

/// @brief clear the selection cache
/// @details This will clear all currently cached selected tables. Normally this is done from 
/// the destructor, but there may be reasons to do this explictly
void TableSelectionDiskCache::clear(){
  for (auto name : getInstance().itsNames) {
    casacore::Table(name,casacore::Table::Delete);
  }
}

TableSelectionDiskCache& TableSelectionDiskCache::getInstance() {
  static TableSelectionDiskCache theInstance;
  return theInstance;
}

