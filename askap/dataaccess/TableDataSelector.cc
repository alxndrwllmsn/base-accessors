/// @file TableDataSelector.cc
/// @brief Implementation of IDataSelector is the table-based case
/// @details
/// TableBasedDataSelector: Class representing a selection of visibility
///                data according to some criterion. This is an
///                implementation of the part of the IDataSelector
///                interface, which can be done with the table selection
///                mechanism in the table based case
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

/// own include
#include <askap/dataaccess/TableDataSelector.h>
#include <askap/dataaccess/DataAccessError.h>
#include <askap/dataaccess/TableTimeStampSelectorImpl.h>

// casa includes
#include <casacore/casa/OS/Path.h>
#include <casacore/tables/DataMan/TiledStManAccessor.h>
#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".dataaccess.tabledataselector");

using namespace askap;
using namespace askap::accessors;
using namespace casa;


/// construct a table selector passing a table/derived info manager
/// via a smart pointer
/// @param[in] msManager a shared pointer to the manager of the measurement set
/// (a derivative of ISubtableInfoHolder)
TableDataSelector::TableDataSelector(const
       boost::shared_ptr<ITableManager const> &msManager) :
       TableInfoAccessor(msManager),
#ifndef ASKAP_DEBUG
       itsDataColumnName(msManager->defaultDataColumnName()),
#endif
       itsChannelSelection(-1,0),itsNFreq(-1)
{
  ASKAPDEBUGASSERT(msManager);
#ifdef ASKAP_DEBUG
  itsDataColumnName = msManager->defaultDataColumnName();
#endif
}

/// Choose a time range. Both start and stop times are given via
/// casacore::MVEpoch object. The reference frame is specified by
/// the DataSource object.
/// @param[in] start the beginning of the chosen time interval
/// @param[in] stop  the end of the chosen time interval
void TableDataSelector::chooseTimeRange(const casacore::MVEpoch &start,
          const casacore::MVEpoch &stop)
{
   itsEpochSelector.reset(new TableTimeStampSelectorImpl<casacore::MVEpoch>(table(),
                          start, stop));
}

/// Choose time range. This method accepts a time range with
/// respect to the origin defined by the DataSource object.
/// Both start and stop times are given as Doubles.
/// The reference frame is the same as for the version accepting
/// MVEpoch and is specified via the DataSource object.
/// @param[in] start the beginning of the chosen time interval
/// @param[in] stop the end of the chosen time interval
void TableDataSelector::chooseTimeRange(casacore::Double start,casacore::Double stop)
{
   itsEpochSelector.reset(new TableTimeStampSelectorImpl<casacore::Double>(table(),
                          start, stop));
}

/// Choose cycles. This is an equivalent of choosing the time range,
/// but the selection is done in integer cycle numbers
/// @param[in] start the number of the first cycle to choose
/// @param[in] stop the number of the last cycle to choose
void TableDataSelector::chooseCycles(casacore::uInt start, casacore::uInt stop)
{
   ASKAPTHROW(DataAccessLogicError, "chooseCycles has not yet been implemented, requested start="<<start
              <<"stop="<<stop);
}

/// Obtain a table expression node for selection. This method is
/// used in the implementation of the iterator to form a subtable
/// obeying the selection criteria specified by the user via
/// IDataSelector interface
///
/// @param[in] conv  a reference to the converter, which is used to sort
///              out epochs used in the selection
const casacore::TableExprNode& TableDataSelector::getTableSelector(const
                    boost::shared_ptr<IDataConverterImpl const> &conv) const
{
   if (itsEpochSelector) {
       /// epoch selection has been done, we need to narrow down the
       /// selection by updating the table expression
       itsEpochSelector->setConverter(conv);
       itsEpochSelector->updateTableExpression(rwTableSelector());
   }
   return rwTableSelector();
}

/// Choose a subset of spectral channels
/// @param[in] nChan a number of spectral channels wanted in the output
/// @param[in] start the number of the first spectral channel to choose
/// @param[in] nAvg a number of adjacent spectral channels to average
///             default is no averaging
void TableDataSelector::chooseChannels(casacore::uInt nChan, casacore::uInt start,
                             casacore::uInt nAvg)
{
   if (nAvg != 1) {
       throw DataAccessLogicError("not yet implemented");
   }
   ASKAPDEBUGASSERT((nChan>0) && (start>=0));
   itsChannelSelection.first = int(nChan);
   itsChannelSelection.second = int(start);
}

/// Choose a subset of frequencies. The reference frame is
/// defined by the DataSource object by default, but can be specified
/// @param[in] nChan a number of spectral channels wanted in the output
/// @param[in] start the frequency of the first spectral channel to
///        choose, given as casacore::MFrequency object, this
///        includes the frequency type and optionally the reference direction
/// @param[in] freqInc an increment in terms of the frequency in the
///        same reference frame as start. This parameter plays
///        the same role as nAvg for chooseChannels, i.e. twice
///        the frequency resolution would average two adjacent channels
void TableDataSelector::chooseFrequencies(casacore::uInt nChan,
         const casacore::MFrequency &start,
         const casacore::MVFrequency &freqInc)
{
   ASKAPDEBUGASSERT((nChan>0) && (start.getValue()>=0));
   ASKAPCHECK(freqInc.near(casacore::MVFrequency(0)), "Non zero frequency increment not yet implemented");
   itsNFreq = nChan;
   itsFreqStart = start;
   itsFreqInc = freqInc.getValue();
}

/// Choose a subset of radial velocities. The reference frame is
/// defined by the DataSource object
/// @param[in] nChan a number of spectral channels wanted in the output
/// @param[in] start the velocity of the first spectral channel to
///        choose (given as casacore::MVRadialVelocity object)
/// @param[in] velInc an increment in terms of the radial velocity in the
///        same reference frame as start. This parameter plays
///        the same role as nAvg for chooseChannels, i.e. twice
///        the velocity resolution would average two adjacent channels
void TableDataSelector::chooseVelocities(casacore::uInt nChan,
         const casacore::MVRadialVelocity &start,
	 const casacore::MVRadialVelocity &velInc)
{
   ASKAPTHROW(DataAccessLogicError, "chooseVelocities has not yet been implemented, requested nChan="<<nChan
              <<"start="<<start<<" velInc="<<velInc);
}

/// Choose polarization.
/// @param[in] pols a string describing the wanted polarization
/// in the output. Allowed values are: I, "IQUV","XXYY","RRLL"
void TableDataSelector::choosePolarizations(const casacore::String &pols)
{
   ASKAPTHROW(DataAccessLogicError, "choosePolarizations has not yet been implemented, requested pols="<<pols);
}

/// @brief work out the tiling of the column in the time/row direction
/// @param[in] string column : name of column to inspect (usually DATA or FLAG)
/// @param[out] tileShape : the tileShape of the column
/// @param[out] hypercubeShape: the shape of the tiled hypercube
/// @return uint: the number of tiles in the row direction, or 0 if invalid
uint TableDataSelector::getTiling(const std::string& column,
        casacore::IPosition& tileShape, casacore::IPosition& hypercubeShape) const
{
    // first check if data is tiled and only has a single non-empty hypercube
    casacore::TableDesc td = table().actualTableDesc();
    const casacore::ColumnDesc& cdesc = td[column];
    const casacore::String dataManType = cdesc.dataManagerType();
    const casacore::String dataManGroup = cdesc.dataManagerGroup();
    if (!dataManType.contains("Tiled")) {
        ASKAPLOG_WARN_STR(logger,"Column is not stored using a Tiled Storage Manager");
        return 0;
    }
    casacore::ROTiledStManAccessor tsm(table(), dataManGroup);
    // try to find first non-empty hypercube and count how many we have
    int hypercube = -1;
    uint nhyper = 0;
    for (uint i = 0; i < tsm.nhypercubes(); i++) {
        if (tsm.getTileShape(i).product() > 0 ) {
            if (hypercube == -1) {
                hypercube = i;
            }
            nhyper++;
        }
    }
    if (nhyper != 1) {
        ASKAPLOG_WARN_STR(logger,"Column does not have a single tiling scheme");
        return 0;
    }
    tileShape = tsm.getTileShape(hypercube);
    hypercubeShape = tsm.getHypercubeShape(hypercube);
    // we implicitly assume that column is tiled over pol, chan and row
    ASKAPASSERT(tileShape.nelements()==3 && hypercubeShape.nelements()==3);
    ASKAPCHECK(tileShape(2) > 0, "Expect a positive number in the tileShape, you have "<<tileShape);
    const uint nRowTiles = hypercubeShape(2) / tileShape(2) +
        (hypercubeShape(2) % tileShape(2) ? 1 : 0);
    return nRowTiles;
}


/// Choose Data Tiles (in time/row direction)
/// @param[in] nTiles number of tiles to select
/// @param[in] start Starting tile number (>=0)
void TableDataSelector::chooseDataTiles(uint nTiles, uint start)
{
    casacore::IPosition tileShape, hyperCubeShape;
    const uint nRowTiles = getTiling(getDataColumnName(),tileShape,hyperCubeShape);
    ASKAPLOG_INFO_STR(logger,"Tile selection nTiles = " << nTiles << ", start = " << start << " -> " <<
                 "row selection nRows = " << nTiles*tileShape(2)<<", start = "<<start*tileShape(2));
    chooseRows(nTiles * tileShape(2), start * tileShape(2));
}

/// @brief choose data column
/// @details This method allows to choose any table column as the visibility
/// data column (e.g. DATA, CORRECTED_DATA, etc). Because this is a
/// table-specific operation, this method is defined in a table-specific
/// selector interface and is not present in IDataSelector (therefore,
/// a dynamic_pointer_cast is likely required).
/// @param[in] dataColumn column name, which contains visibility data
void TableDataSelector::chooseDataColumn(const std::string &dataColumn)
{
   itsDataColumnName = dataColumn;
}

/// @brief clone a selector
/// @details The same selector can be used to create a number of iterators.
/// Selector stores a name of the data column to use and, therefore, it can
/// be changed after some iterators are created. To avoid bugs due to this
/// reference semantics, the iterator will clone selector in its constructor.
/// @note This functionality is not exposed to the end user, which
/// normally interacts with the IDataSelector class only. This is because
/// cloning is done at the low level (e.g. inside the iterator)
boost::shared_ptr<ITableDataSelectorImpl const> TableDataSelector::clone() const
{
  return boost::shared_ptr<ITableDataSelectorImpl const>(new TableDataSelector(*this));
}

/// @brief obtain the name of data column
/// @details This method returns the current name of the data column, set
/// either in the constructor or by the chooseDataColumn method
/// @return the name of the data column
const std::string& TableDataSelector::getDataColumnName() const throw()
{
  return itsDataColumnName;
}

/// @brief check whether channel selection has been done
/// @details By default all channels are selected. However, if chooseChannels
/// has been called, less channels are returned. This method returns true if
/// this is the case and false otherwise.
/// @return true, if a subset of channels has been selected
bool TableDataSelector::channelsSelected() const throw()
{
  return itsChannelSelection.first >= 0;
}

/// @brief obtain channel selection
/// @details By default all channels are selected. However, if chooseChannels
/// has been called, less channels are returned by the accessor. This method
/// returns the number of channels and the first channel (in the full sample)
/// selected. If the first element of the pair is negative, no channel-based
/// selection has been done. This is also checked by channelsSelected method,
/// which is probably a prefered way to do this check to retain the code clarity.
/// @return a pair, the first element gives the number of channels selected and
/// the second element gives the start channel (0-based)
std::pair<int,int> TableDataSelector::getChannelSelection() const throw()
{
  return itsChannelSelection;
}

/// @brief check whether frequency selection has been done
/// @details By default all channels are selected. However, if chooseFrequencies
/// has been called, less channels are returned. This method returns true if
/// this is the case and false otherwise.
/// @return true, if a subset of frequencies has been selected
bool TableDataSelector::frequenciesSelected() const throw()
{
  return itsNFreq > 0;
}

/// @brief obtain frequency selection
/// @details By default all channels are selected. However, if chooseFrequencies
/// has been called, less channels are returned by the accessor.
/// @return the number of channels, the start frequency and increment (Hz) and the reference type
std::tuple<int,casacore::MFrequency,double> TableDataSelector::getFrequencySelection() const throw()
{
    return std::tuple<int,casacore::MFrequency,double>(itsNFreq,itsFreqStart,itsFreqInc);
}

/// @brief specify the name for caching the selection
/// @details The (row) selection can be cached for repeated use, specify the name
/// of the cache. If unspecified or empty, no caching will be done
/// @param[in] string name : name of the cache file
void  TableDataSelector::setSelectionCacheName(const std::string& name) {
    itsSelectionCacheName = name;
}

/// @brief return the name for caching the selection
/// @details The (row) selection can be cached for repeated use, return the name
/// of the cache. If empty, no caching is active
/// @return string name : name of the cache file
const std::string& TableDataSelector::getSelectionCacheName() const {
    return itsSelectionCacheName;
}

  /// @brief return the table name
  /// @details The table name may be useful in setting the selection cache name.
  /// This function returns the base name of the table.
  /// @return string: the basename of the table (without directories)
  std::string TableDataSelector::getTableName() const {
    return casacore::Path(table().tableName()).baseName();
  }
