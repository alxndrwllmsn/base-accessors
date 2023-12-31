/// @file
/// @brief table-based implementation of the calibration solution source
/// @details This implementation reads calibration solutions from and writes to a casa table
/// Main functionality is implemented in the corresponding TableCalSolutionFiller class.
/// This class creates an instance of the MemCalSolutionAccessor with the above mentioned filler
/// when a writeable accessor is requested. Read-only functionality is implemented in the
/// base class.
///
/// @copyright (c) 2011 CSIRO
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
/// @author Max Voronkov <Maxim.Voronkov@csiro.au>

// own includes
#include <askap/calibaccess/TableCalSolutionSource.h>
#include <askap/calibaccess/TableCalSolutionFiller.h>
#include <askap/calibaccess/MemCalSolutionAccessor.h>

// casa includes
#include <casacore/measures/TableMeasures/TableMeasDesc.h>
#include <casacore/measures/TableMeasures/TableMeasRefDesc.h>
#include <casacore/measures/TableMeasures/TableMeasValueDesc.h>
#include <casacore/tables/Tables/ScaColDesc.h>
#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/Tables/TableUtil.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>
#include <casacore/measures/Measures/MCEpoch.h>
#include <casacore/tables/Tables/TableError.h>
#include <casacore/tables/Tables/SetupNewTab.h>
#include <casacore/tables/Tables/TableDesc.h>
#include <casacore/casa/OS/RegularFile.h>
#include <casacore/casa/OS/Directory.h>
#include <casacore/casa/OS/File.h>



namespace askap {

namespace accessors {

/// @brief constructor using a table defined explicitly
/// @details
/// @param[in] tab table to work with
/// @param[in] nAnt maximum number of antennas
/// @param[in] nBeam maximum number of beams
/// @param[in] nChan maximum number of channels
TableCalSolutionSource::TableCalSolutionSource(const casacore::Table &tab, const casacore::uInt nAnt,
         const casacore::uInt nBeam, const casacore::uInt nChan) : TableHolder(tab),
   TableCalSolutionConstSource(tab), itsNAnt(nAnt), itsNBeam(nBeam), itsNChan(nChan) {}

/// @brief constructor using a file name
/// @details The table is opened for writing
/// @param[in] name table file name
/// @param[in] nAnt maximum number of antennas
/// @param[in] nBeam maximum number of beams
/// @param[in] nChan maximum number of channels
TableCalSolutionSource::TableCalSolutionSource(const std::string &name, const casacore::uInt nAnt,
         const casacore::uInt nBeam, const casacore::uInt nChan) :
   TableHolder(casacore::Table()), TableCalSolutionConstSource(table()), itsNAnt(nAnt), itsNBeam(nBeam), itsNChan(nChan)
{
  try {
     table() = casacore::Table(name,casacore::Table::Update);
  }
  catch (...) {
     // we couldn't just opened an existing table
     try {
        casacore::SetupNewTable maker(name, casacore::TableDesc(), casacore::Table::New);
        table() = casacore::Table(maker);
     }
     catch (const casacore::TableError &te) {
        ASKAPTHROW(DataAccessError,"Unable create a new table for calibration solutions with the name="<<name<<
                   ". AipsError: " << te.what());
     }
  }
}



/// @brief obtain a solution ID to store new solution
/// @details This method provides a solution ID for a new solution. It must
/// be called before any write operation (one needs a writable accessor to
/// write the actual solution and to get this accessor one needs an ID).
/// @param[in] time time stamp of the new solution in seconds since MJD of 0.
/// @return solution ID
long TableCalSolutionSource::newSolutionID(const double time) {
   if (!table().actualTableDesc().isColumn("TIME")) {
       // this is a new table, we need to create new TIME column
       casacore::ScalarColumnDesc<casacore::Double> timeColDesc("TIME",
           "Time stamp when the calibration solution was obtained");
       table().addColumn(timeColDesc);
       casacore::TableMeasRefDesc measRef(casacore::MEpoch::UTC);
       casacore::TableMeasValueDesc measVal(table().actualTableDesc(), "TIME");
       casacore::TableMeasDesc<casacore::MEpoch> mepochCol(measVal, measRef);
       mepochCol.write(table());
   }
   const casacore::rownr_t newRow = table().nrow();
   table().addRow(1);
   ASKAPDEBUGASSERT(newRow < table().nrow());
   casacore::ScalarMeasColumn<casacore::MEpoch> bufCol(table(),"TIME");
   const casacore::Quantity qTime(time,"s");
   const casacore::MEpoch epoch(qTime,casacore::MEpoch::UTC);
   bufCol.put(newRow, epoch);
   return static_cast<long>(newRow);
}

/// @brief obtain a writeable accessor for a given solution ID
/// @details This method returns a shared pointer to the solution accessor, which
/// can be used to both read the parameters and write them back. If a solution with
/// the given ID doesn't exist, an exception is thrown. Existing solutions with undefined
/// parameters are managed via validity flags of gains, leakages and bandpasses
/// @param[in] id solution ID to access
/// @return shared pointer to an accessor object
boost::shared_ptr<ICalSolutionAccessor> TableCalSolutionSource::rwSolution(const long id) const {
   ASKAPCHECK((id >= 0) && (static_cast<long>(table().nrow()) > id), "Requested solution id="<<id<<" is not in the table");
   boost::shared_ptr<TableCalSolutionFiller> filler(new TableCalSolutionFiller(table(),id,itsNAnt,
          itsNBeam, itsNChan));
   ASKAPDEBUGASSERT(filler);
   boost::shared_ptr<MemCalSolutionAccessor> acc(new MemCalSolutionAccessor(filler,false));
   ASKAPDEBUGASSERT(acc);
   return acc;
}

/// @brief helper method to remove an old table
/// @details It just deletes the given table, which allows to create a new one
/// from scratch (this functionality is used if one needs to overwrite the previous
/// solution).
/// @param[in] fname file name to delete
/// @param[in] removeIfNotTable if true, the file is removed even if it is not a table.
/// An exception is thrown in this case if this parameter is false.
void TableCalSolutionSource::removeOldTable(const std::string &fname, const bool removeIfNotTable)
{
  if (casacore::TableUtil::canDeleteTable(fname,false)) {
      casacore::TableUtil::deleteTable(fname, false);
  } else {
     // check that the table simply doesn't exist
     ASKAPCHECK(!tableExists(fname), "Unable to delete existing table "<<fname);
     casacore::File tmpFile(fname);
     if (tmpFile.exists()) {
         ASKAPCHECK(removeIfNotTable,
                    "TableCalSolutionSource::removeOldTable: File or directory "<<fname<<
                    " exists, but it is not a table - unable to remove");
         // we need to remove the file with the given name
         if (tmpFile.isDirectory()) {
             casacore::Directory dir(fname);
             dir.removeRecursive();
         } else {
             ASKAPASSERT(tmpFile.isRegular());
             casacore::RegularFile rf(fname);
             rf.remove();
         }
     }
  }
}


} // namespace accessors

} // namespace askap
