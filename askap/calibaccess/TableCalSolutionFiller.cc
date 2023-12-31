/// @file
/// @brief solution filler reading required cubes from casa table
/// @details This is an example of a class which knows how to fill buffers
/// of MemCalSolutionAccessor. The cubes with calibration information are read
/// from (and written to) a casa table. The table has the following columns:
/// TIME, GAIN, GAIN_VALID, LEAKAGE, LEAKAGE_VALID, BANDPASS and BANDPASS_VALID.
/// This class is initialised with the reference row, which corresponds to the time
/// requested by the user. If there are gains, leakages or bandpasses defined for
/// a given row, they are read. Otherwise, a backward search is performed to find
/// the first defined value. An exception is thrown if the top of the table is reached.
/// If a new entry needs to be created, the given numbers of antennas and beams are used.
///
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

#include <map>
#include <askap/calibaccess/TableCalSolutionFiller.h>

namespace askap {

namespace accessors {

/// @brief construct the object and link it to the given table
/// @details read-only operation is assumed
/// @param[in] tab  table to use
/// @param[in] row reference row
TableCalSolutionFiller::TableCalSolutionFiller(const casa::Table& tab, const long row) : TableHolder(tab),
       TableBufferManager(tab), itsNAnt(0), itsNBeam(0), itsNChan(0), itsRefRow(row), itsGainsRow(-1),
       itsLeakagesRow(-1), itsBandpassesRow(-1), itsBPLeakagesRow(-1), itsIonoParamsRow(-1)
{
  ASKAPCHECK((itsRefRow >= 0) && (itsRefRow <= static_cast<long>(table().nrow())), "Requested calibration solution ID = "<<itsRefRow<<" is outside calibration table");
  // this is the reading case, we can use either of itsNAnt, itsNBeam or itsNChan to assert this initial state (they should be all 0). The whole job is encapsulated in
  // the  isReadOnly method
  ASKAPDEBUGASSERT(isReadOnly());
}

/// @brief construct the object and link it to the given table
/// @details Maximum allowed numbers of antennas, beams and spectral channels are
/// set by this constructor which is essential for read-write operations (i.e. new
/// table entries may need to be created
/// @param[in] tab  table to use
/// @param[in] row reference row
/// @param[in] nAnt maximum number of antennas
/// @param[in] nBeam maximum number of beams
/// @param[in] nChan maximum number of channels
TableCalSolutionFiller::TableCalSolutionFiller(const casa::Table& tab, const long row, const casa::uInt nAnt,
          const casa::uInt nBeam, const casa::uInt nChan) : TableHolder(tab),
       TableBufferManager(tab), itsNAnt(nAnt), itsNBeam(nBeam), itsNChan(nChan), itsRefRow(row), itsGainsRow(-1),
       itsLeakagesRow(-1), itsBandpassesRow(-1), itsBPLeakagesRow(-1), itsIonoParamsRow(-1)
{
  ASKAPCHECK((itsRefRow >= 0) && (itsRefRow <= static_cast<long>(table().nrow())), "Requested calibration solution ID = "<<itsRefRow<<" is outside calibration table");
  // this is the writing case, so numbers of antennas, beams and channels should be positive
  ASKAPCHECK(itsNAnt > 0, "TableCalSolutionFiller needs to know the number of antennas to be able to setup new table rows");
  ASKAPCHECK(itsNBeam > 0, "TableCalSolutionFiller needs to know the number of beams to be able to setup new table rows");
  ASKAPCHECK(itsNChan > 0, "TableCalSolutionFiller needs to know the number of spectral channels to be able to setup new table rows");
  ASKAPDEBUGASSERT(!isReadOnly());
}

/// @brief helper method to check that the filler is initialised for read only access
/// @return true, if the filler is expected to do read only operations
/// @note Look back until the last defined record is only done for read-only access. Read-write access
/// overwrites whatever row is requested
bool TableCalSolutionFiller::isReadOnly() const
{
   ASKAPDEBUGASSERT((itsNAnt == 0) == (itsNBeam == 0));
   ASKAPDEBUGASSERT((itsNAnt == 0) == (itsNChan == 0));
   return (itsNAnt == 0);
}

/// @brief check for gain solution
/// @return true, if there is no gain solution, false otherwise
bool TableCalSolutionFiller::noGain() const
{
  return !columnExists("GAIN");
}

/// @brief check for leakage solution
/// @return true, if there is no leakage solution, false otherwise
bool TableCalSolutionFiller::noLeakage() const
{
  return !columnExists("LEAKAGE");
}

/// @brief check for bandpass solution
/// @return true, if there is no bandpass solution, false otherwise
bool TableCalSolutionFiller::noBandpass() const
{
  return !columnExists("BANDPASS");
}

/// @brief check for bpleakage solution
/// @return true, if there is no bpleakage solution, false otherwise
bool TableCalSolutionFiller::noBPLeakage() const
{
  return !columnExists("BPLEAKAGE");
}

/// @brief check for ionospheric solution
/// @return true, if there is no ionospheric solution, false otherwise
bool TableCalSolutionFiller::noIonosphere() const
{
  return !columnExists("IONOSPHERE");
}

/// @brief helper method to check that the given column exists
/// @param[in] name column name
/// @return true if the given column exists
bool TableCalSolutionFiller::columnExists(const std::string &name) const
{
  const std::map<std::string, bool>::const_iterator it = itsColumnExistsCache.find(name);
  if (it != itsColumnExistsCache.end()) {
      return it->second;
  } else {
      const bool exists = table().actualTableDesc().isColumn(name);
      itsColumnExistsCache[name] = exists;
      return exists;
  }
}

/// @brief gains filler
/// @details
/// @param[in] gains pair of cubes with gains and validity flags (to be resised to 2 x nAnt x nBeam)
void TableCalSolutionFiller::fillGains(std::pair<casacore::Cube<casacore::Complex>, casacore::Cube<casacore::Bool> > &gains) const
{
  // cellDefined should not be called if noGain returns true according to C++ evaluation rules.
  const bool needToCreateGains = noGain() || !cellDefined<casa::Complex>("GAIN", static_cast<casacore::rownr_t>(itsRefRow));
  if (!isReadOnly() && needToCreateGains) {
      ASKAPDEBUGASSERT(itsGainsRow < 0);
      gains.first.resize(2, itsNAnt, itsNBeam);
      gains.first.set(1.);
      gains.second.resize(2, itsNAnt, itsNBeam);
      gains.second.set(false);
      itsGainsRow = itsRefRow;
  } else {
     // this is the case wwhere either the table is read-only or there is a need to read data first

     if (itsGainsRow < 0) {
         itsGainsRow = findDefinedCube("GAIN");
     }
     ASKAPASSERT(itsGainsRow>=0);
     if (itsGainsRow != itsRefRow) {
         // backwards search should only be possible in the read-only mode
         ASKAPDEBUGASSERT(isReadOnly());
         ASKAPDEBUGASSERT(needToCreateGains);
     }
     ASKAPCHECK(cellDefined<casa::Bool>("GAIN_VALID", static_cast<casacore::rownr_t>(itsGainsRow)),
         "Wrong format of the calibration table: GAIN element should always be accompanied by GAIN_VALID");
     readCube(gains.first, "GAIN", static_cast<casacore::rownr_t>(itsGainsRow));
     readCube(gains.second, "GAIN_VALID", static_cast<casacore::rownr_t>(itsGainsRow));
  }
  ASKAPCHECK(gains.first.shape() == gains.second.shape(), "GAIN and GAIN_VALID cubes are expected to have the same shape");
}

/// @brief leakage filler
/// @details
/// @param[in] leakages pair of cubes with leakages and validity flags (to be resised to 2 x nAnt x nBeam)
void TableCalSolutionFiller::fillLeakages(std::pair<casacore::Cube<casacore::Complex>, casacore::Cube<casacore::Bool> > &leakages) const
{
  // cellDefined should not be called if noLeakage returns true according to C++ evaluation rules.
  const bool needToCreateLeakage = noLeakage() || !cellDefined<casa::Complex>("LEAKAGE", static_cast<casacore::rownr_t>(itsRefRow));
  if (!isReadOnly() && needToCreateLeakage) {
      ASKAPDEBUGASSERT(itsLeakagesRow < 0);
      leakages.first.resize(2, itsNAnt, itsNBeam);
      leakages.first.set(0.);
      leakages.second.resize(2, itsNAnt, itsNBeam);
      leakages.second.set(false);
      itsLeakagesRow = itsRefRow;
  } else {
     // this is the case wwhere either the table is read-only or there is a need to read data first
     if (itsLeakagesRow < 0) {
         itsLeakagesRow = findDefinedCube("LEAKAGE");
     }
     ASKAPASSERT(itsLeakagesRow>=0);
     if (itsLeakagesRow != itsRefRow) {
         // backwards search should only be possible in the read-only mode
         ASKAPDEBUGASSERT(isReadOnly());
         ASKAPDEBUGASSERT(needToCreateLeakage);
     }
     ASKAPCHECK(cellDefined<casa::Bool>("LEAKAGE_VALID", static_cast<casacore::rownr_t>(itsLeakagesRow)),
         "Wrong format of the calibration table: LEAKAGE element should always be accompanied by LEAKAGE_VALID");
     readCube(leakages.first, "LEAKAGE", static_cast<casacore::rownr_t>(itsLeakagesRow));
     readCube(leakages.second, "LEAKAGE_VALID", static_cast<casacore::rownr_t>(itsLeakagesRow));
  }
  ASKAPCHECK(leakages.first.shape() == leakages.second.shape(), "LEAKAGE and LEAKAGE_VALID cubes are expected to have the same shape");
}

/// @brief bandpass filler
/// @details
/// @param[in] bp pair of cubes with bandpasses and validity flags (to be resised to (2*nChan) x nAnt x nBeam)
void TableCalSolutionFiller::fillBandpasses(std::pair<casacore::Cube<casacore::Complex>, casacore::Cube<casacore::Bool> > &bp) const
{
  // cellDefined should not be called if noBandpass returns true according to C++ evaluation rules.
  const bool needToCreateBandpass = noBandpass() || !cellDefined<casa::Complex>("BANDPASS", static_cast<casacore::rownr_t>(itsRefRow));
  if (!isReadOnly() && needToCreateBandpass) {
      ASKAPDEBUGASSERT(itsBandpassesRow < 0);
      bp.first.resize(2 * itsNChan, itsNAnt, itsNBeam);
      bp.first.set(1.);
      bp.second.resize(2 * itsNChan, itsNAnt, itsNBeam);
      bp.second.set(false);
      itsBandpassesRow = itsRefRow;
  } else {
     // this is the case wwhere either the table is read-only or there is a need to read data first
     if (itsBandpassesRow < 0) {
         itsBandpassesRow = findDefinedCube("BANDPASS");
     }
     ASKAPASSERT(itsBandpassesRow>=0);
     if (itsBandpassesRow != itsRefRow) {
         // backwards search should only be possible in the read-only mode
         ASKAPDEBUGASSERT(isReadOnly());
         ASKAPDEBUGASSERT(needToCreateBandpass);
     }
     ASKAPCHECK(cellDefined<casa::Bool>("BANDPASS_VALID", static_cast<casacore::rownr_t>(itsBandpassesRow)),
         "Wrong format of the calibration table: BANDPASS element should always be accompanied by BANDPASS_VALID");
     readCube(bp.first, "BANDPASS", static_cast<casacore::rownr_t>(itsBandpassesRow));
     readCube(bp.second, "BANDPASS_VALID", static_cast<casacore::rownr_t>(itsBandpassesRow));
  }
  ASKAPCHECK(bp.first.shape() == bp.second.shape(), "BANDPASS and BANDPASS_VALID cubes are expected to have the same shape");
}

/// @brief bpleakage filler
/// @details
/// @param[in] bp pair of cubes with bandpasses and validity flags (to be resized to (2*nChan) x nAnt x nBeam)
void TableCalSolutionFiller::fillBPLeakages(std::pair<casacore::Cube<casacore::Complex>, casacore::Cube<casacore::Bool> > &bpleakages) const
{
  // cellDefined should not be called if noBPLeakage returns true according to C++ evaluation rules.
  const bool needToCreateBPLeakage = noBPLeakage() || !cellDefined<casa::Complex>("BPLEAKAGE", static_cast<casacore::rownr_t>(itsRefRow));
  if (!isReadOnly() && needToCreateBPLeakage) {
      ASKAPDEBUGASSERT(itsBPLeakagesRow < 0);
      bpleakages.first.resize(2 * itsNChan, itsNAnt, itsNBeam);
      bpleakages.first.set(0.);
      bpleakages.second.resize(2 * itsNChan, itsNAnt, itsNBeam);
      bpleakages.second.set(false);
      itsBPLeakagesRow = itsRefRow;
  } else {
     // this is the case wwhere either the table is read-only or there is a need to read data first
     if (itsBPLeakagesRow < 0) {
         itsBPLeakagesRow = findDefinedCube("BPLEAKAGE");
     }
     ASKAPASSERT(itsBPLeakagesRow>=0);
     if (itsBPLeakagesRow != itsRefRow) {
         // backwards search should only be possible in the read-only mode
         ASKAPDEBUGASSERT(isReadOnly());
         ASKAPDEBUGASSERT(needToCreateBPLeakage);
     }
     ASKAPCHECK(cellDefined<casa::Bool>("BPLEAKAGE_VALID", static_cast<casacore::rownr_t>(itsBPLeakagesRow)),
         "Wrong format of the calibration table: BPLEAKAGE element should always be accompanied by BPLEAKAGE_VALID");
     readCube(bpleakages.first, "BPLEAKAGE", static_cast<casacore::rownr_t>(itsBPLeakagesRow));
     readCube(bpleakages.second, "BPLEAKAGE_VALID", static_cast<casacore::rownr_t>(itsBPLeakagesRow));
  }
  ASKAPCHECK(bpleakages.first.shape() == bpleakages.second.shape(), "BPLEAKAGE and BPLEAKAGE_VALID cubes are expected to have the same shape");
}

/// @brief ionospheric parameters filler
/// @details
/// @param[in] pair of cubes with ionospheric parameters and validity flags (to be resised to 1 x nAnt x nBeam)
void TableCalSolutionFiller::fillIonoParams(std::pair<casacore::Cube<casacore::Complex>,
                                                      casacore::Cube<casacore::Bool> > &params) const
{
  // cellDefined should not be called if noIonosphere returns true according to C++ evaluation rules.
  const bool needToCreateIono = noIonosphere() || !cellDefined<casa::Complex>("IONOSPHERE", static_cast<casacore::rownr_t>(itsRefRow));
  if (!isReadOnly() && needToCreateIono) {
      ASKAPDEBUGASSERT(itsIonoParamsRow < 0);
      params.first.resize(1, itsNAnt, itsNBeam);
      params.first.set(0.);
      params.second.resize(1, itsNAnt, itsNBeam);
      params.second.set(false);
      itsIonoParamsRow = itsRefRow;
  } else {
     // this is the case wwhere either the table is read-only or there is a need to read data first

     if (itsIonoParamsRow < 0) {
         itsIonoParamsRow = findDefinedCube("IONOSPHERE");
     }
     ASKAPASSERT(itsIonoParamsRow>=0);
     if (itsIonoParamsRow != itsRefRow) {
         // backwards search should only be possible in the read-only mode
         ASKAPDEBUGASSERT(isReadOnly());
         ASKAPDEBUGASSERT(needToCreateIono);
     }
     ASKAPCHECK(cellDefined<casa::Bool>("IONOSPHERE_VALID", static_cast<casacore::rownr_t>(itsIonoParamsRow)),
         "Wrong format of the calibration table: IONOSPHERE element should always be accompanied by IONOSPHERE_VALID");
     readCube(params.first, "IONOSPHERE", static_cast<casacore::rownr_t>(itsIonoParamsRow));
     readCube(params.second, "IONOSPHERE_VALID", static_cast<casacore::rownr_t>(itsIonoParamsRow));
  }
  ASKAPCHECK(params.first.shape() == params.second.shape(), "IONOSPHERE and IONOSPHERE_VALID cubes are expected to have the same shape");
}

/// @brief gains writer
/// @details
/// @param[in] gains pair of cubes with gains and validity flags (should be 2 x nAnt x nBeam)
void TableCalSolutionFiller::writeGains(const std::pair<casacore::Cube<casacore::Complex>, casacore::Cube<casacore::Bool> > &gains) const
{
  ASKAPASSERT(itsGainsRow>=0);
  ASKAPCHECK(gains.first.shape() == gains.second.shape(), "The cubes with gains and validity flags are expected to have the same shape");
  writeCube(gains.first, "GAIN", static_cast<casacore::rownr_t>(itsGainsRow));
  writeCube(gains.second, "GAIN_VALID", static_cast<casacore::rownr_t>(itsGainsRow));
}

/// @brief leakage writer
/// @details
/// @param[in] leakages pair of cubes with leakages and validity flags (should be 2 x nAnt x nBeam)
void TableCalSolutionFiller::writeLeakages(const std::pair<casacore::Cube<casacore::Complex>, casacore::Cube<casacore::Bool> > &leakages) const
{
  ASKAPASSERT(itsLeakagesRow>=0);
  ASKAPCHECK(leakages.first.shape() == leakages.second.shape(), "The cubes with leakages and validity flags are expected to have the same shape");
  writeCube(leakages.first, "LEAKAGE", static_cast<casacore::rownr_t>(itsLeakagesRow));
  writeCube(leakages.second, "LEAKAGE_VALID", static_cast<casacore::rownr_t>(itsLeakagesRow));
}

/// @brief bandpass writer
/// @details
/// @param[in] bp pair of cubes with bandpasses and validity flags (should be (2*nChan) x nAnt x nBeam)
void TableCalSolutionFiller::writeBandpasses(const std::pair<casacore::Cube<casacore::Complex>, casacore::Cube<casacore::Bool> > &bp) const
{
  ASKAPASSERT(itsBandpassesRow>=0);
  ASKAPCHECK(bp.first.shape() == bp.second.shape(), "The cubes with bandpasses and validity flags are expected to have the same shape");
  writeCube(bp.first, "BANDPASS", static_cast<casacore::rownr_t>(itsBandpassesRow));
  writeCube(bp.second, "BANDPASS_VALID", static_cast<casacore::rownr_t>(itsBandpassesRow));
}

/// @brief bpleakage writer
/// @details
/// @param[in] bpleakages pair of cubes with bpleakages and validity flags (should be (2*nChan) x nAnt x nBeam)
void TableCalSolutionFiller::writeBPLeakages(const std::pair<casacore::Cube<casacore::Complex>, casacore::Cube<casacore::Bool> > &bpleakages) const
{
  ASKAPASSERT(itsBPLeakagesRow>=0);
  ASKAPCHECK(bpleakages.first.shape() == bpleakages.second.shape(), "The cubes with bandpass leakage and validity flags are expected to have the same shape");
  writeCube(bpleakages.first, "BPLEAKAGE", static_cast<casacore::rownr_t>(itsBPLeakagesRow));
  writeCube(bpleakages.second, "BPLEAKAGE_VALID", static_cast<casacore::rownr_t>(itsBPLeakagesRow));
}

/// @brief ionospheric parameters writer
/// @details
/// @param[in] pair of cubes with ionospheric parameters and validity flags (should be 1 x nAnt x nBeam)
void TableCalSolutionFiller::writeIonoParams(const std::pair<casacore::Cube<casacore::Complex>, casacore::Cube<casacore::Bool> > &params) const
{
  ASKAPASSERT(itsIonoParamsRow>=0);
  ASKAPCHECK(params.first.shape() == params.second.shape(), "The cubes with ionospheric parameters and validity flags are expected to have the same shape");
  writeCube(params.first, "IONOSPHERE", static_cast<casacore::rownr_t>(itsIonoParamsRow));
  writeCube(params.second, "IONOSPHERE_VALID", static_cast<casacore::rownr_t>(itsIonoParamsRow));
}

/// @brief find first defined cube searching backwards
/// @details This assumes that the table rows are given in the time order. If the cell at the reference row
/// doesn't have a cube defined, the search is continued up to the top of the table. An exception is thrown
/// if no defined cube has been found.
/// @param[in] name column name
/// @return row number for a defined cube
/// @note The code always returns non-negative number.
long TableCalSolutionFiller::findDefinedCube(const std::string &name) const
{
  for (long tempRow = itsRefRow; tempRow >= 0; --tempRow) {
       if (cellDefined<casacore::Complex>(name, static_cast<casacore::rownr_t>(tempRow))) {
           return tempRow;
       }
  }
  ASKAPTHROW(AskapError, "Unable to find valid element in column "<<name<<" at row "<<itsRefRow<<" or earlier");
}


} // namespace accessors

} // namespace askap
