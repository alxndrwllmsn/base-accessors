/// @file TableConstDataIterator.cc
///
/// @brief Implementation of IConstDataIterator in the table-based case
/// @details
/// TableConstDataIterator: Allow read-only iteration across preselected data. Each
/// iteration step is represented by the IConstDataAccessor interface.
/// This is an implementation in the table-based case.
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

#include <askap_accessors.h>

/// ASKAPsoft includes
#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapError.h>
#include <askap/profile/AskapProfiler.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>
#include <casacore/scimath/Mathematics/SquareMatrix.h>
#include <casacore/measures/Measures/MeasFrame.h>
#include <casacore/measures/Measures/MCFrequency.h>
#include <casacore/casa/Arrays/Slicer.h>
#include <casacore/casa/Arrays/IPosition.h>
#include <casacore/casa/Arrays/MatrixMath.h>
#include <casacore/casa/Arrays/ArrayMath.h>

/// Local package
#include <askap/dataaccess/TableConstDataIterator.h>
#include <askap/dataaccess/DataAccessError.h>
#include <askap/dataaccess/DirectionConverter.h>
#include <askap/dataaccess/TableSelectionDiskCache.h>

ASKAP_LOGGER(logger, "dataaccess");

using namespace casa;
using namespace askap;
using namespace askap::accessors;

/// @param[in] msManager a manager of the measurement set to use
/// @param[in] sel shared pointer to selector
/// @param[in] conv shared pointer to converter
/// @param[in] cacheSize a number of uvw machines in the cache (default is 1)
/// @param[in] tolerance pointing direction tolerance in radians, exceeding which leads
/// to initialisation of a new UVW Machine
/// @param[in] maxChunkSize maximum number of rows per accessor
TableConstDataIterator::TableConstDataIterator(
            const boost::shared_ptr<ITableManager const> &msManager,
            const boost::shared_ptr<ITableDataSelectorImpl const> &sel,
            const boost::shared_ptr<IDataConverterImpl const> &conv,
            size_t cacheSize, double tolerance,
            casacore::uInt maxChunkSize) :
        TableInfoAccessor(msManager),
        // it is essential that accessor is initialised after cache parameters!
	    itsUVWCacheSize(cacheSize), itsUVWCacheTolerance(tolerance),
	    itsAccessor(*this),
#ifndef ASKAP_DEBUG
        itsSelector(sel->clone()),
	    itsConverter(conv->clone()),
#endif
	    itsMaxChunkSize(maxChunkSize)
{
  ASKAPDEBUGASSERT(conv);
  ASKAPDEBUGASSERT(sel);
  #ifdef ASKAP_DEBUG
    itsConverter = conv->clone();
    itsSelector  = sel->clone();
  #endif
  init();

}

/// Restart the iteration from the beginning
void TableConstDataIterator::init()
{
  ASKAPTRACE("TableConstDataIterator::init");

  itsCurrentTopRow=0;
  itsCurrentDataDescID=-100; // this value can't be in the table,
                              // therefore it is a flag of a new data descriptor
  itsCurrentFieldID = -100; // this value can't be in the table,
                            // therefore it is a flag of a new field ID
  // by default use FIELD_ID column if it exists, otherwise use time to select
  // pointings
  itsUseFieldID = table().tableDesc().isColumn("FIELD_ID");

  if (itsTabIterator.isNull()) {
    const casacore::TableExprNode &exprNode =
                itsSelector->getTableSelector(itsConverter);
    if (exprNode.isNull()) {
      itsTabIterator=casacore::TableIterator(table(),"TIME",
      casacore::TableIterator::Ascending,casacore::TableIterator::NoSort);
    } else {
      const std::string & name = itsSelector->getSelectionCacheName();
      const Table sel = name.empty() ? table()(exprNode) : TableSelectionDiskCache::table(name, exprNode);
      itsTabIterator = casacore::TableIterator(sel, "TIME", casacore::TableIterator::Ascending,
        casacore::TableIterator::NoSort);
    }
  } else {
    itsTabIterator.reset();
  }
  itsChannelsSelected = false;
  itsFlagData = false;
  setUpIteration();
}

/// operator* delivers a reference to data accessor (current chunk)
/// @return a reference to the current chunk
const IConstDataAccessor& TableConstDataIterator::operator*() const
{
  return itsAccessor;
}

/// Checks whether there are more data available.
/// @return True if there are more data available
casacore::Bool TableConstDataIterator::hasMore() const throw()
{
  if (!itsTabIterator.pastEnd()) {
      return true;
  }
  if (itsCurrentTopRow+itsNumberOfRows<itsCurrentIteration.nrow()) {
      return true;
  }
  return false;
}

/// advance the iterator one step further
/// @return True if there are more data (so constructions like
///         while(it.next()) {} are possible)
casacore::Bool TableConstDataIterator::next()
{
  ASKAPTRACE("TableConstDataIterator::next");
  itsCurrentTopRow+=itsNumberOfRows;
  if (itsCurrentTopRow>=itsCurrentIteration.nrow()) {
      ASKAPDEBUGASSERT(!itsTabIterator.pastEnd());
      itsCurrentTopRow=0;
      // need to advance table iterator further
      itsTabIterator.next();
      if (!itsTabIterator.pastEnd()) {
          setUpIteration();
      }
  } else {
      casacore::rownr_t remainder=itsCurrentIteration.nrow()-itsCurrentTopRow;
      itsNumberOfRows=remainder<=itsMaxChunkSize ?
                      remainder : itsMaxChunkSize;
      itsAccessor.invalidateIterationCaches();
      // itsDirectionCache don't need invalidation because the time is the same
      // as for the previous iteration

      // determine whether DATA_DESC_ID is uniform in the whole chunk
      // and reduce itsNumberOfRows if necessary
      makeUniformDataDescID();

      // determine whether FIELD_ID is uniform in the whole chunk
      // and reduce itsNumberOfRows if necessary
      // invalidate direction cache if necessary.
      // do nothing if itsUseFieldID is false
      makeUniformFieldID();
  }
  return hasMore();
}

/// setup accessor for a new iteration of the table iterator
void TableConstDataIterator::setUpIteration()
{
  itsCurrentIteration=itsTabIterator.table();
  itsAccessor.invalidateIterationCaches();

  itsNumberOfRows=itsCurrentIteration.nrow()<=itsMaxChunkSize ?
                  itsCurrentIteration.nrow() : itsMaxChunkSize;

  if ((itsDirectionCache.isValid() || itsParallacticAngleCache.isValid())
       && itsCurrentDataDescID>=0) {
      // extra checks make sense if the cache is valid (and this means it
      // has been used before)
      const casacore::MEpoch epoch = currentEpoch();
      const casacore::uInt spWindow = currentSpWindowID();
      const bool newField = itsUseFieldID ? false : subtableInfo().getField().newField(epoch);
      // a case where fieldID changes is dealt with separately.
      const IFeedSubtableHandler &feedSubtable = subtableInfo().getFeed();
      if ( newField || !subtableInfo().getAntenna().allEquatorial()) {
          itsParallacticAngleCache.invalidate();
      }
      if ( newField || ((!subtableInfo().getAntenna().allEquatorial() ||
           feedSubtable.newBeamDetails(epoch,spWindow)) &&
           !feedSubtable.allBeamOffsetsZero(epoch,spWindow))) {
              itsDirectionCache.invalidate();
              // rotated uvw depend on the direction (phase centres)
              itsAccessor.invalidateRotatedUVW();
              // itsDishPointingCache doesn't depend on feeds
              if (newField) {
                  itsDishPointingCache.invalidate();
              }
      }
  }
  // retrieve the number of channels and polarizations from the table
  if (itsNumberOfRows) {
      // determine whether DATA_DESC_ID is uniform in the whole chunk
      // and reduce itsNumberOfRows if necessary
      // set up visibility cube shape if necessary
      makeUniformDataDescID();

      // determine whether FIELD_ID is uniform in the whole chink
      // and reduce itsNumberOfRows if necessary
      // invalidate direction cache if necessary.
      // do nothing if itsUseFieldID is false
      makeUniformFieldID();
  } else {
      itsNumberOfChannels = 0;
      itsNumberOfPols = 0;
      itsCurrentDataDescID = -100;
      itsCurrentFieldID = -100;
      itsDirectionCache.invalidate();
      // rotated uvw depends on the direction (phase centres)
      itsAccessor.invalidateRotatedUVW();
      itsParallacticAngleCache.invalidate();
      itsDishPointingCache.invalidate();
  }
}

/// @brief method ensures that the chunk has uniform DATA_DESC_ID
/// @details This method reduces itsNumberOfRows to achieve a
/// uniform DATA_DESC_ID reading for all rows in the current chunk.
/// The resulting itsNumberOfRows will be 1 or more.
/// theirAccessor's spectral axis cache is reset if new DATA_DESC_ID is
/// different from itsCurrentDataDescID
/// This method also sets up itsNumberOfPols and itsNumberOfChannels
/// when DATA_DESC_ID changes (and therefore at the first run as well)
void TableConstDataIterator::makeUniformDataDescID()
{
  ASKAPDEBUGASSERT(itsNumberOfRows);
  ASKAPDEBUGASSERT(itsCurrentTopRow+itsNumberOfRows<=
                    itsCurrentIteration.nrow());

  ROScalarColumn<Int> dataDescCol(itsCurrentIteration,"DATA_DESC_ID");
  const Int newDataDescID=dataDescCol(itsCurrentTopRow);
  ASKAPDEBUGASSERT(newDataDescID>=0);
  if (itsCurrentDataDescID!=newDataDescID) {
      itsAccessor.invalidateSpectralCaches();
      itsCurrentDataDescID=newDataDescID;
      if (itsDirectionCache.isValid()) {
          // if-statement, because it is pointless to do further checks in the
          // case when the cache is already invalid due to
          // the time change. In addition, checks require an access to the table,
          // which we want to avoid if, e.g., we don't need pointing direction
          // at all
          const casacore::uInt spWindow = currentSpWindowID();
          const casacore::MEpoch epoch = currentEpoch();
          const IFeedSubtableHandler &feedSubtable = subtableInfo().getFeed();
          if (!feedSubtable.allBeamOffsetsZero(epoch,spWindow)) {
              if (feedSubtable.newBeamDetails(epoch,spWindow)) {
                  itsDirectionCache.invalidate();
                  // rotated uvw depend on the direction (phase centres)
                  itsAccessor.invalidateRotatedUVW();
              }
          }
      }

      // determine the shape of the visibility cube
      ROArrayColumn<Complex> visCol(itsCurrentIteration,getDataColumnName());
      const casacore::IPosition shape=visCol.shape(itsCurrentTopRow);
      ASKAPASSERT(shape.size() && (shape.size()<3));
      itsNumberOfPols=shape[0];
      itsNumberOfChannels=shape.size()>1?shape[1]:1;
      ASKAPDEBUGASSERT(itsSelector);
      if (itsSelector->channelsSelected()) {
          // validity checks that selection doesn't extend beyond the channels available
          const std::pair<int,int> chanSelection = itsSelector->getChannelSelection();
          ASKAPCHECK(itsNumberOfChannels >= casacore::uInt(chanSelection.first + chanSelection.second),
               "Channel selection from "<<chanSelection.second+1<<" to "<<chanSelection.first+
               chanSelection.second<<" (1-based) extends beyond "<<itsNumberOfChannels<<
               " channel(s) available in  the dataset");
      }
  }

  // if we're selecting by frequency we need to redo the channel selection when time or DATA_DESC_ID changes
  if (itsSelector->frequenciesSelected()) {
      // reset the channel selection
      itsChannelsSelected = false;
  }

  for (uInt row=1;row<itsNumberOfRows;++row) {
       if (dataDescCol(row+itsCurrentTopRow)!=itsCurrentDataDescID) {
           itsNumberOfRows=row;
	   break;
       }
  }
}

/// @brief method ensures that the chunk has a uniform FIELD_ID
/// @details This method reduces itsNumberOfRows until FIELD_ID is
/// the same for all rows in the current chunk. The resulting
/// itsNumberOfRows will be 1 or more. If itsUseFieldID is false,
/// the method returns without doing anything. itsAccessor's direction
/// cache is reset if new FIELD_ID is different from itsCurrentFieldID
/// (and it sets it up at the first run as well)
void TableConstDataIterator::makeUniformFieldID()
{
  if (itsUseFieldID) {
      ASKAPDEBUGASSERT(itsNumberOfRows);
      ASKAPDEBUGASSERT(itsCurrentTopRow+itsNumberOfRows<=
                       itsCurrentIteration.nrow());

      ROScalarColumn<Int> fieldIDCol(itsCurrentIteration,"FIELD_ID");
      const Int newFieldID=fieldIDCol(itsCurrentTopRow);
      ASKAPDEBUGASSERT(newFieldID>=0);
      if (newFieldID != itsCurrentFieldID) {
          itsCurrentFieldID = newFieldID;
          itsDirectionCache.invalidate();
          // rotated uvw depend on the direction (phase centres)
          itsAccessor.invalidateRotatedUVW();
          itsParallacticAngleCache.invalidate();
          itsDishPointingCache.invalidate();
      }
      // break the iteration if necessary
      for (uInt row=1;row<itsNumberOfRows;++row) {
           if (fieldIDCol(row+itsCurrentTopRow)!=itsCurrentFieldID) {
               itsNumberOfRows=row;
	           break;
           }
      }
  }
}

/// @brief read an array column of the table into a cube
/// @details populate the buffer provided with the information
/// read in the current iteration. This method is templated and can be
/// used for both visibility and flag data fillers.
/// @param[in] cube a reference to the nPol x nChannel x nRow buffer
///            cube to fill with the information from table
/// @param[in] columnName a name of the column to read
template<typename T>
void TableConstDataIterator::fillCube(casacore::Cube<T> &cube,
               const std::string &columnName) const
{
  // Setup a slicer to extract the specified channel range only
  const Slicer chanSlicer(Slice(), Slice(startChannel(), nChannel()));
  ROArrayColumn<T> tableCol(itsCurrentIteration, columnName);
  if (itsNumberOfRows == itsCurrentIteration.nrow()) {
    tableCol.getColumn(chanSlicer, cube, True);
  } else {
    const Slicer rowSlicer(Slice(itsCurrentTopRow, itsNumberOfRows));
    tableCol.getColumnRange(rowSlicer, chanSlicer, cube, True);
  }
}

/// populate the buffer of visibilities with the values of current
/// iteration
/// @param[out] vis a reference to the nPol x nChannel x nRow buffer
///            cube to fill with the complex visibility data
void TableConstDataIterator::fillVisibility(casacore::Cube<casacore::Complex> &vis) const
{
  fillCube(vis,getDataColumnName());
}

/// @brief read flagging information
/// @details populate the buffer of flags with the information
/// read in the current iteration
/// @param[in] flag a reference to the nPol x nChannel x nRow buffer
///            cube to fill with the flag information (each element has
///            bool type)
void TableConstDataIterator::fillFlag(casacore::Cube<casacore::Bool> &flag) const
{
  fillCube(flag,"FLAG");

  if (itsFlagData) {
    flag = true;
  } else {
    if (itsCurrentIteration.tableDesc().isColumn("FLAG_ROW")) {
      ROScalarColumn<Bool> flagRowCol(itsCurrentIteration,"FLAG_ROW");
      for (uInt row=0; row<itsNumberOfRows; ++row) {
        // flag the row if row flag is set
        if (flagRowCol(row + itsCurrentTopRow)) {
          flag.xyPlane(row) = true;
        }
      }
    }
  }
}

/// populate the buffer of noise figures with the values of current
/// iteration
/// @param[in] noise a reference to the nPol x nChannel x nRow buffer
///            cube to be filled with the noise figures
void TableConstDataIterator::fillNoise(casacore::Cube<casacore::Complex> &noise) const
{
  ASKAPDEBUGASSERT(itsSelector);
  const casacore::uInt nChan = nChannel();
  const casacore::uInt startChan = startChannel();

  // default action first - just resize the cube and assign 1.
  noise.resize(itsNumberOfPols, nChan, itsNumberOfRows);
  noise.set(casacore::Complex(1.,1.));
  // if the sigma spectrum exists, use those sigmas to fill the noise cube
  static bool once=true;
  if (table().tableDesc().isColumn("SIGMA_SPECTRUM")) {
      // noise is given per channel and polarisation
      // Setup a slicer to extract the specified channel range only
      if (once) {
        ASKAPLOG_INFO_STR(logger,"Using SIGMA_SPECTRUM column for noise");
        once=false;
      }
      const Slicer chanSlicer(Slice(),Slice(startChan,nChan));
      casa::Matrix<Float> buf(itsNumberOfPols,nChan);
      ROArrayColumn<Float> sigmaCol(itsCurrentIteration,"SIGMA_SPECTRUM");
      for (uInt row = 0; row<itsNumberOfRows; ++row) {
#ifdef ASKAP_DEBUG
           const casa::IPosition shape = sigmaCol.shape(row);
           ASKAPDEBUGASSERT(shape.size()==2);
           ASKAPDEBUGASSERT((shape[0] == casa::Int(itsNumberOfPols)) &&
                            (shape[1] == casa::Int(itsNumberOfChannels)));
#endif
           sigmaCol.getSlice(row+itsCurrentTopRow,chanSlicer,buf,False);

           // SIGMA_SPECTRUM is ordered (pol,chan)
           for (casa::uInt chan=0; chan<nChan; chan++) {
                for (casa::uInt pol=0; pol<itsNumberOfPols; pol++) {
                     // same noise for both real and imaginary parts
                     const casa::Float val = buf(pol,chan);
                     noise(pol,chan,row) = casa::Complex(val,val);
                }
           }
      } // loop over rows
  } // if-statement checking that SIGMA_SPECTRUM column is present
  else if (table().tableDesc().isColumn("SIGMA")) {
      if (once) {
        ASKAPLOG_INFO_STR(logger,"Using SIGMA column for noise");
        once=false;
      }
      ROArrayColumn<Float> sigmaCol(itsCurrentIteration,"SIGMA");
      ASKAPDEBUGASSERT((sigmaCol.shape(0).size()==1) && sigmaCol.shape(0)[0] == casacore::Int(itsNumberOfPols));
      const Slicer rowSlicer(Slice(itsCurrentTopRow, itsNumberOfRows));
      const casacore::Matrix<Float> buf = sigmaCol.getColumnRange(rowSlicer);
      for (uInt row = 0; row<itsNumberOfRows; ++row) {
          // noise is given per polarisation, same for all spectral channels
          for (casacore::uInt pol=0; pol<itsNumberOfPols; ++pol) {
              // same noise for both real and imaginary parts
              const casacore::Float val = buf(pol,row);
              // nChan is 1 for imaging, so indexing is faster than slicing
              for (uInt chan = 0; chan< nChan; ++chan) {
                  noise(pol,chan,row) = casacore::Complex(val,val);
              }
          }
      } // loop over rows
  } // if-statement checking that SIGMA column is present
}

/// populate the buffer with uvw
/// @param[in] uvw a reference to vector of rigid vectors (3 elemets,
///            u,v and w for each row) to fill
void TableConstDataIterator::fillUVW(casacore::Vector<casacore::RigidVector<casacore::Double, 3> >&uvw) const
{
  uvw.resize(itsNumberOfRows);

  ROArrayColumn<Double> uvwCol(itsCurrentIteration,"UVW");
  // temporary buffer
  Vector<Double> buf(3);
  for (uInt row=0;row<itsNumberOfRows;++row) {
#ifdef ASKAP_DEBUG
       const casacore::IPosition shape=uvwCol.shape(row);
       ASKAPDEBUGASSERT(shape.size()==1);
       ASKAPDEBUGASSERT(shape[0]==3);
#endif // ASKAP_DEBUG
       // extract data record for this row, no resizing
       uvwCol.get(row+itsCurrentTopRow,buf,False);
       uvw(row) = buf;
  }
}

/// @brief obtain a current spectral window ID
/// @details This method obtains a spectral window ID corresponding to the
/// current data description ID and tests its validity
/// @return current spectral window ID
casacore::uInt TableConstDataIterator::currentSpWindowID() const
{
  ASKAPDEBUGASSERT(itsCurrentDataDescID>=0);
  const int spWindowIndex = subtableInfo().getDataDescription().
                            getSpectralWindowID(itsCurrentDataDescID);
  if (spWindowIndex<0) {
      ASKAPTHROW(DataAccessError,"A negative spectral window index ("<<
              spWindowIndex<<") is encountered for Data Description ID="<<
	      itsCurrentDataDescID);
  }
  return static_cast<const uInt>(spWindowIndex);
}

/// @brief obtain a current polarisation ID
/// @details This method obtains a polarisation ID corresponding to the current
/// data description ID and tests its validity
/// @return current polarisation ID
casacore::uInt TableConstDataIterator::currentPolID() const
{
  ASKAPDEBUGASSERT(itsCurrentDataDescID>=0);
  const int polIndex = subtableInfo().getDataDescription().
                            getPolarizationID(itsCurrentDataDescID);
  if (polIndex<0) {
      ASKAPTHROW(DataAccessError,"A negative polarisation index ("<<
              polIndex<<") is encountered for Data Description ID="<<
	      itsCurrentDataDescID);
  }
  return static_cast<const uInt>(polIndex);
}

/// @brief obtain a reference direction for the current iteration
/// @details Currently we assume that the dish pointing centre stays
/// fixed for the whole chunk. We break the iteration, if necessary
/// to achieve this. This helper method extracts the reference direction
/// from the FIELD subtable using either FIELD_ID, or current time if
/// the former is not supported by the main table.
/// @return a reference to direction measure
const casacore::MDirection& TableConstDataIterator::getCurrentReferenceDir() const
{
  const IFieldSubtableHandler &fieldSubtable = subtableInfo().getField();
  if (itsUseFieldID) {
      ASKAPCHECK(itsCurrentFieldID>=0, "Elements of FIELD_ID column should be 0 or positive. You have "<<
                 itsCurrentFieldID);
      return fieldSubtable.getReferenceDir(itsCurrentFieldID);
  }
  const casacore::MEpoch epoch = currentEpoch();
  return fieldSubtable.getReferenceDir(epoch);
}

/// @brief obtain selected range of channels
/// @details A subset spectral channels can be selected for this iterator to work with.
/// This method returns the number of channels and the first selected channel.
/// @return a pair, first element is the number of channels, second is the first channel
/// in the full cube
std::pair<casacore::uInt, casacore::uInt> TableConstDataIterator::getChannelRange() const
{
  ASKAPDEBUGASSERT(itsSelector);
  if (!itsChannelsSelected) {
      itsFlagData = true;

      if (itsSelector->frequenciesSelected()) {
          const std::tuple<int,casacore::MFrequency,double> freqSel = itsSelector->getFrequencySelection();
          // cannot do multiple channels yet
          ASKAPCHECK(std::get<0>(freqSel)<=1, "Can only do a single channel in frequency selection mode");
          // convert frequency in requested frame to MS frame
          // Using antenna 0 and antenna pointing (= field direction) as reference (or direction ref in MFrequency)
          // Note this differs from imager which uses current phase centre direction in freq conversion
          const casacore::MFrequency freqMeas = std::get<1>(freqSel);
          casacore::MeasRef<casacore::MFrequency> freqRef = freqMeas.getRef();
          const casacore::Measure *pMeas = freqRef.getFrame().direction();
          // If the MFrequency in freqSel has a reference direction use that, otherwise use pointing
          casacore::MDirection velDir = (pMeas ? MDirection(pMeas) : getCurrentReferenceDir());
          casacore::MeasFrame frame(MEpoch(currentEpoch()),subtableInfo().getAntenna().getPosition(0),velDir);
          const ITableSpWindowHolder& spWindowSubtable=subtableInfo().getSpWindow();
          const casacore::MFrequency::Types dataType =
            casacore::MFrequency::castType(spWindowSubtable.getReferenceFrame(currentSpWindowID()).getType());
          casacore::MFrequency::Types selType = MFrequency::castType(freqRef.getType());
          if (selType == casacore::MFrequency::Undefined) {
              selType = dataType;
          }

          const casacore::MFrequency::Ref refin(dataType,frame); // the frame of the input channels
          const casacore::MFrequency::Ref refout(selType,frame); // the frame desired
          casacore::MFrequency::Convert backw(refout,refin); // from desired to input

          const MVFrequency requiredFreq = backw(std::get<1>(freqSel)).getValue();
          // Now find corresponding channel
          const casacore::Vector<casacore::Double> dataFreqs(spWindowSubtable.getFrequencies(currentSpWindowID()));
          // assuming linear freq scale
          itsNumberOfChannelsSelected = 1;
          itsStartChannelSelected = 0;
          const uint nFreq = dataFreqs.nelements();
          if (nFreq > 1) {
              const double freqInc = dataFreqs(1) - dataFreqs(0);
              ASKAPDEBUGASSERT(freqInc != 0);
              ASKAPCHECK(abs((dataFreqs(nFreq-1)-dataFreqs(0))/((nFreq-1)*freqInc)-1)<0.001,
                "Frequency axis non-linear, cannot do frequency selection with current code");

              const double channel = (requiredFreq - dataFreqs(0)) / freqInc;
              // for now just use nearest channel, but could do linear interpolation between nearest two
              const int nearestChannel = std::lrint(channel);
              if (nearestChannel >= 0 && nearestChannel < nFreq) {
                  itsStartChannelSelected = (uint)nearestChannel;
                  itsFlagData = false;
              } else {
                  if (nearestChannel >= nFreq) {
                      itsStartChannelSelected = nFreq - 1;
                  }
              }
          }
      } else {
          const std::pair<int, int> chanSelection = itsSelector->getChannelSelection();
          itsNumberOfChannelsSelected = itsSelector->channelsSelected() ?
                                   casacore::uInt(chanSelection.first) : itsNumberOfChannels;
          itsStartChannelSelected = itsSelector->channelsSelected() ?
                                   casacore::uInt(chanSelection.second) : 0;
          ASKAPDEBUGASSERT(itsNumberOfChannelsSelected + itsStartChannelSelected <= itsNumberOfChannels);
          itsFlagData = false;
      }
      itsChannelsSelected = true;
  }

  return std::pair<casacore::uInt, casacore::uInt>(itsNumberOfChannelsSelected,itsStartChannelSelected);
}

/// @brief test if given row is selected with supplied TableExprNode
/// @details This method tests row selection on the underlying table, adjusting
/// the specified row number for the current state of iteration
/// @return true if row is selected
casacore::Bool TableConstDataIterator::isSelected(const casacore::TableExprNode& ten, casacore::uInt row)
{
    ASKAPASSERT(!ten.isNull());
    ASKAPASSERT(ten.nrow() == table().nrow());
    ASKAPASSERT(row + itsCurrentTopRow < itsCurrentIteration.nrow());
    const rownr_t baseRow = itsCurrentIteration.rowNumbers()(row + itsCurrentTopRow);
    return ten.getBool(baseRow);
}

/// @brief fill the buffer with the polarisation types
/// @param[in] stokes a reference to a vector to be filled
void TableConstDataIterator::fillStokes(casacore::Vector<casacore::Stokes::StokesTypes> &stokes) const
{
  const ITablePolarisationHolder& polSubtable = subtableInfo().getPolarisation();

  ASKAPDEBUGASSERT(itsCurrentDataDescID>=0);
  const casacore::uInt polID = currentPolID();
  ASKAPASSERT(polSubtable.nPol(polID) == nPol());
  stokes = polSubtable.getTypes(polID).copy();
}

/// populate the buffer with frequencies
/// @param[in] freq a reference to a vector to fill
void TableConstDataIterator::fillFrequency(casacore::Vector<casacore::Double> &freq) const
{
  ASKAPDEBUGASSERT(itsConverter);
  const ITableSpWindowHolder& spWindowSubtable=subtableInfo().getSpWindow();
  ASKAPDEBUGASSERT(itsCurrentDataDescID>=0);
  const casacore::uInt spWindowID = currentSpWindowID();

  const casacore::uInt nChan = nChannel();
  const casacore::uInt startChan = startChannel();

  // for the time being we don't do the short-cut if a subset of channels
  // is selected without any conversion. In principle it is possible, but
  // we need to take care of constness as taking a slice is not a const
  // operation.
  if (itsConverter->isVoid(spWindowSubtable.getReferenceFrame(spWindowID),
	                   spWindowSubtable.getFrequencyUnit()) && !itsSelector->channelsSelected()
                       && !itsSelector->frequenciesSelected()) {
      // the conversion is void, i.e. table units/frame are exactly what
      // we need for output. This simplifies things a lot.
      freq.reference(spWindowSubtable.getFrequencies(spWindowID));
      if (itsNumberOfChannels!=freq.nelements()) {
          ASKAPTHROW(DataAccessError,"The measurement set has bad or corrupted "<<
	      "SPECTRAL_WINDOW subtable. The number of spectral channels for data "<<
	      itsNumberOfChannels<<" doesn't match the number of channels in the "<<
	      "frequency axis ("<<freq.nelements()<<")");
      }
  } else {
      // have to process element by element as a conversion is required
      const casacore::MEpoch epoch=currentEpoch();
      // always use the dish pointing centre, rather than a pointing centre
      // of each individual feed for frequency conversion. The error is not
      // huge. If this code will ever work for SKA, this may need to be changed.
      // Currently use the FIELD table, not the actual pointing. It is probably
      // correct to use the phase centre for conversion as opposed to the
      // pointing centre.
      const casacore::MDirection& antReferenceDir = getCurrentReferenceDir();
      // currently use the position of the first antenna for convertion.
      // we may need some average position + a check that they are close
      // enough to throw an exception if someone gives a VLBI measurement set.
      itsConverter->setMeasFrame(casacore::MeasFrame(epoch, subtableInfo().
                     getAntenna().getPosition(0), antReferenceDir));
      freq.resize(nChan);
      for (uInt ch=0;ch<nChan;++ch) {
           freq[ch]=itsConverter->frequency(spWindowSubtable.getFrequency(
	                                    spWindowID,ch+startChan));
      }
  }
}

/// @return the time stamp
casacore::Double TableConstDataIterator::getTime() const
{
  // add additional checks in debug mode
  #ifdef ASKAP_DEBUG
   ROScalarColumn<Double> timeCol(itsCurrentIteration,"TIME");
   Double time=timeCol(itsCurrentTopRow);
    Vector<Double> allTimes=timeCol.getColumnRange(Slicer(IPosition(1,
                       itsCurrentTopRow),IPosition(1,itsNumberOfRows)));
    for (uInt row=0;row<allTimes.nelements();++row)
         if (allTimes[row]!=time) {
             throw DataAccessLogicError("Time column is not homogeneous for each "
	                            "DataAccessor. This shouldn't happend");
	 }
  #endif
  // end of additional checks

  ROScalarMeasColumn<MEpoch> timeMeasCol(itsCurrentIteration,"TIME");
  return itsConverter->epoch(timeMeasCol(itsCurrentTopRow));
}


/// populate the buffer with IDs of the first antenna
/// @param[in] ids a reference to a vector to fill
void TableConstDataIterator::fillAntenna1(casacore::Vector<casacore::uInt>& ids) const
{
  fillVectorOfIDs(ids,"ANTENNA1");
}

/// populate the buffer with IDs of the second antenna
/// @param[in] ids a reference to a vector to fill
void TableConstDataIterator::fillAntenna2(casacore::Vector<casacore::uInt> &ids) const
{
  fillVectorOfIDs(ids,"ANTENNA2");
}

/// populate the buffer with IDs of the first feed
/// @param[in] ids a reference to a vector to fill
void TableConstDataIterator::fillFeed1(casacore::Vector<casacore::uInt> &ids) const
{
  fillVectorOfIDs(ids,"FEED1");
}

/// populate the buffer with IDs of the second feed
/// @param[in] ids a reference to a vector to fill
void TableConstDataIterator::fillFeed2(casacore::Vector<casacore::uInt> &ids) const
{
  fillVectorOfIDs(ids,"FEED2");
}


/// @brief a helper method to read a column with IDs of some sort
/// @details It reads the column of casacore::Int and fills a Vector of
/// casacore::uInt. A check to ensure all numbers are non-negative is done
/// in the debug mode.
/// @param[in] ids a reference to a vector to fill
/// @param[in] name a name of the column to read
void TableConstDataIterator::fillVectorOfIDs(casacore::Vector<casacore::uInt> &ids,
                     const casacore::String &name) const
{
  ROScalarColumn<Int> col(itsCurrentIteration,name);
  ids.resize(itsNumberOfRows);
  Vector<Int> buf=col.getColumnRange(Slicer(IPosition(1,
                      itsCurrentTopRow),IPosition(1,itsNumberOfRows)));
  ASKAPDEBUGASSERT(buf.nelements()==ids.nelements());
  // need a copy because the type is different. There are no
  // appropriate cast operators for casacore::Vectors
  Vector<Int>::const_iterator ci=buf.begin();
  Vector<uInt>::iterator it=ids.begin();
  for (; ci!=buf.end() && it!=ids.end() ; ++ci,++it) {
       ASKAPDEBUGASSERT(*ci>=0);
	   *it=static_cast<uInt>(*ci);
  }
}

/// @brief an alternative way to get the time stamp
/// @details This method uses the accessor to get cached time stamp. It
/// is returned as an epoch measure.
casacore::MEpoch TableConstDataIterator::currentEpoch() const
{
  ASKAPDEBUGASSERT(itsConverter);
  return itsConverter->epochMeasure(itsAccessor.time());
}

/// @brief Fill internal buffer with parallactic angles
/// @details This buffer holds parallactic angles for all antennas. The buffer
/// is invalidated when the time changes for an alt-az array, for an equatorial
/// array it happens only if the pointing changes.
/// @param[in] angles a reference to a vector to be filled
void TableConstDataIterator::fillParallacticAngleCache(casacore::Vector<casacore::Double> &angles) const
{
  angles.resize(subtableInfo().getAntenna().getNumberOfAntennas());
  ASKAPDEBUGASSERT(angles.size());
  if (subtableInfo().getAntenna().allEquatorial()) {
      angles.set(0.);
  } else {

  const casacore::MEpoch epoch=currentEpoch();

  // we need a separate converter for parallactic angle calculations
  DirectionConverter dirConv((casacore::MDirection::Ref(casacore::MDirection::AZEL)));
  dirConv.setMeasFrame(epoch);

  // we currently use FIELD table to get the pointing direction. This table
  // does not depend on the antenna.
  const casacore::MDirection& antReferenceDir = getCurrentReferenceDir();


  for (casacore::uInt ant = 0; ant<angles.size(); ++ant) {
       const casacore::String &antMount=subtableInfo().getAntenna().
                                          getMount(ant);

       if (antMount == "ALT-AZ" || antMount == "alt-az")  {
           casacore::MDirection celestialPole;
           celestialPole.set(MDirection::Ref(MDirection::HADEC));
           dirConv.setMeasFrame(casacore::MeasFrame(subtableInfo().getAntenna().
                              getPosition(ant),epoch));
           const casacore::Double pAngle=dirConv(antReferenceDir).
                               positionAngle(dirConv(celestialPole).getValue());
           angles[ant] = pAngle;
       } else if (antMount == "FIXED"  ||  antMount == "fixed") {
           // LOFAR has a fixed antenna mount.
           angles[ant] = 0.;
       } else if (antMount == "X-Y"  ||  antMount == "x-y") {
           angles[ant] = 0.;
       } else if (antMount != "EQUATORIAL" && antMount != "equatorial") {
           ASKAPTHROW(DataAccessError,"Unknown mount type "<<antMount<<
              " for antenna "<<ant);
       }
    }
  } // if all equatorial
}


/// @brief Fill internal buffer with the pointing directions
/// @details  The layout of this buffer is the same as the layout of
/// the FEED subtable for current time and spectral window.
/// getAntennaIDs and getFeedIDs methods of the
/// subtable handler can be used to unwrap this 1D array.
/// The buffer can be invalidated if the time changes (i.e. for an alt-az array),
/// for an equatorial array this happends only if the FEED or FIELD subtable
/// are time-dependent or if FIELD_ID changes
/// @param[in] dirs a reference to a vector to fill
void TableConstDataIterator::fillDirectionCache(casacore::Vector<casacore::MVDirection> &dirs) const
{
  // the code fills both pointing directions and position angles. For ASKAP, it would
  // probably be a bit faster if we split these two operations between two methods, as
  // position angle will be fixed and will not need as much updating as the pointing.

  // threshold (in radians) beyond which a proper parallactic angle rotation is done
  // it is hard coded at the moment
  const double parallacticAngleThreshold = 1e-9;

  const casacore::Vector<casacore::Double> &parallacticAngles = itsParallacticAngleCache.value(*this,
                 &TableConstDataIterator::fillParallacticAngleCache);

  const IFeedSubtableHandler &feedSubtable = subtableInfo().getFeed();

  const casacore::MEpoch epoch=currentEpoch();
  ASKAPDEBUGASSERT(itsCurrentDataDescID>=0);
  const casacore::uInt spWindowID = currentSpWindowID();
  // antenna and feed IDs here are those in the FEED subtable, rather than
  // in the current accessor
  const casacore::Vector<casacore::Int> &antIDs = feedSubtable.getAntennaIDs(epoch,spWindowID);

  dirs.resize(antIDs.nelements());

  // we currently use FIELD table to get the pointing direction. This table
  // does not depend on the antenna.
  const casacore::MDirection& antReferenceDir = getCurrentReferenceDir();

  const casacore::Vector<casacore::RigidVector<casacore::Double, 2> > &offsets =
               feedSubtable.getAllBeamOffsets(epoch,spWindowID);

  for (casacore::uInt element=0;element<antIDs.nelements();++element) {
       const casacore::uInt ant=antIDs[element];

       // if we decide to be paranoid about performance, we can add a method
       // to the converter to test whether antenna position and/or epoch are
       // really required to the requested convertion. Because the antenna
       // position are cached, the overhead of the present straightforward
       // approach should be relatively minor.
       itsConverter->setMeasFrame(casacore::MeasFrame(epoch,subtableInfo().
                     getAntenna().getPosition(ant)));

       casacore::RigidVector<casacore::Double, 2> offset = offsets[element];
       ASKAPDEBUGASSERT(ant<parallacticAngles.nelements());
       const casacore::Double posAngle = parallacticAngles[ant];

       if (std::abs(posAngle)>parallacticAngleThreshold)  {
           // need to do a proper parallactic angle rotation
           casacore::SquareMatrix<casacore::Double, 2>
                       rotMatrix(casacore::SquareMatrix<casacore::Double, 2>::General);
           const casacore::Double cpa=cos(posAngle);
           const casacore::Double spa=sin(posAngle);
           rotMatrix(0,0)=cpa;
           rotMatrix(0,1)=-spa;
           rotMatrix(1,0)=spa;
           rotMatrix(1,1)=cpa;
           offset*=rotMatrix;
       }
       casacore::MDirection feedPointingCentre(antReferenceDir);
       // x direction is flipped to convert az-el type frame to ra-dec
       feedPointingCentre.shift(casacore::MVDirection(-offset(0),
                             offset(1)),casacore::True);
       itsConverter->direction(feedPointingCentre,dirs[element]);
  }
}


/// fill the buffer with the pointing directions of the first antenna/feed
/// @param[in] dirs a reference to a vector to fill
void TableConstDataIterator::fillPointingDir1(
                          casacore::Vector<casacore::MVDirection> &dirs) const
{
  const casacore::Vector<casacore::uInt> &feedIDs=itsAccessor.feed1();
  const casacore::Vector<casacore::uInt> &antIDs=itsAccessor.antenna1();
  fillVectorOfPointings(dirs,antIDs,feedIDs);
}

/// fill the buffer with the pointing directions of the second antenna/feed
/// @param[in] dirs a reference to a vector to fill
void TableConstDataIterator::fillPointingDir2(
                          casacore::Vector<casacore::MVDirection> &dirs) const
{
  const casacore::Vector<casacore::uInt> &feedIDs=itsAccessor.feed2();
  const casacore::Vector<casacore::uInt> &antIDs=itsAccessor.antenna2();
  fillVectorOfPointings(dirs,antIDs,feedIDs);
}

/// fill the buffer with the position angles of the first antenna/feed
/// @param[in] angles a reference to vector to be filled
void TableConstDataIterator::fillFeed1PA(casacore::Vector<casacore::Float> &angles) const
{
  const casacore::Vector<casacore::uInt> &feedIDs=itsAccessor.feed1();
  const casacore::Vector<casacore::uInt> &antIDs=itsAccessor.antenna1();
  fillVectorOfPositionAngles(angles,antIDs,feedIDs);
}

/// fill the buffer with the position angles of the second antenna/feed
/// @param[in] angles a reference to vector to be filled
void TableConstDataIterator::fillFeed2PA(casacore::Vector<casacore::Float> &angles) const
{
  const casacore::Vector<casacore::uInt> &feedIDs=itsAccessor.feed2();
  const casacore::Vector<casacore::uInt> &antIDs=itsAccessor.antenna2();
  fillVectorOfPositionAngles(angles,antIDs,feedIDs);
}


/// @brief a helper method to get dish pointings
/// @details fillDishPointing1 and fillDishPointing2 methods do very
/// similar operations, which differ only by the antennaIDs used.
/// This method encapsulated these common operations.
/// @note fillVectorOfPointings computes pointing directions for
/// individual feeds, not for the centre of the dish as this method
/// does
/// @param[in] dirs a reference to a vector to fill
/// @param[in] antIDs a vector with antenna IDs
void TableConstDataIterator::fillVectorOfDishPointings(casacore::Vector<casacore::MVDirection> &dirs,
               const casacore::Vector<casacore::uInt> &antIDs) const
{
  ASKAPDEBUGASSERT(itsNumberOfRows == antIDs.nelements());
  const casacore::Vector<casacore::MVDirection> &dishPointingCache = itsDishPointingCache.
                      value(*this,&TableConstDataIterator::fillDishPointingCache);
  dirs.resize(itsNumberOfRows);
  for (casacore::uInt row=0; row<itsNumberOfRows; ++row) {
       ASKAPDEBUGASSERT(antIDs[row] < dishPointingCache.nelements());
       dirs[row] = dishPointingCache[antIDs[row]];
  }
}

/// @brief fill the buffer with the pointing directions for the first antenna centre
/// @details The difference from fillPointingDir1 is that no feed offset is applied.
/// @param[in] dirs a reference to a vector to fill
void TableConstDataIterator::fillDishPointing1(casacore::Vector<casacore::MVDirection> &dirs) const
{
  const casacore::Vector<casacore::uInt> &antIDs=itsAccessor.antenna1();
  fillVectorOfDishPointings(dirs,antIDs);
}

/// @brief fill the buffer with the pointing directions for the second antenna centre
/// @details The difference from fillPointingDir2 is that no feed offset is applied.
/// @param[in] dirs a reference to a vector to fill
void TableConstDataIterator::fillDishPointing2(casacore::Vector<casacore::MVDirection> &dirs) const
{
  const casacore::Vector<casacore::uInt> &antIDs=itsAccessor.antenna2();
  fillVectorOfDishPointings(dirs,antIDs);
}

/// @brief fill the buffer with the dish pointing directions
/// @details The difference from fillDirectionCache is that
/// this method computes the pointing directions for the dish centre, not for
/// individual feeds (or synthetic beams, strictly speaking). The number of elements
/// in the buffer equals to the number of antennas. This is also different from
/// fillDirectionCache, which projects feeds to the same 1D array as well.
/// @note At this stage we use FIELD subtable to get the pointing directions.
/// Therefore, these directions do not depend on antenna/feed. This method writes
/// the same value for all elements of the array. It will be used for both antennas
/// in the pair.
/// @param[in] dirs a reference to a vector to fill
void TableConstDataIterator::fillDishPointingCache(casacore::Vector<casacore::MVDirection> &dirs) const
{
  ASKAPDEBUGASSERT(itsConverter);
  const casacore::MEpoch epoch = currentEpoch();

  dirs.resize(subtableInfo().getAntenna().getNumberOfAntennas());

  // we currently use FIELD table to get the pointing direction. This table
  // does not depend on the antenna. However, the reference frame can introduce such
  // a dependence (i.e. a large array and AZEL frame requested)
  const casacore::MDirection& antReferenceDir = getCurrentReferenceDir();

  for (casacore::uInt ant = 0; ant<dirs.nelements(); ++ant) {
       // if we decide to be paranoid about performance, we can add a method
       // to the converter to test whether antenna position and/or epoch are
       // really required to the requested convertion. Because the antenna
       // position are cached, the overhead of the present straightforward
       // approach should be relatively minor.
       itsConverter->setMeasFrame(casacore::MeasFrame(epoch,subtableInfo().
                  getAntenna().getPosition(ant)));
       itsConverter->direction(antReferenceDir,dirs[ant]);
  }
}

/// @brief A helper method to fill a given vector with pointing directions.
/// @details fillPointingDir1 and fillPointingDir2 methods do very similar
/// operations, which differ only by the feedIDs and antennaIDs used.
/// This method encapsulates these common operations
/// @param[in] dirs a reference to a vector to fill
/// @param[in] antIDs a vector with antenna IDs
/// @param[in] feedIDs a vector with feed IDs
void TableConstDataIterator::fillVectorOfPointings(
               casacore::Vector<casacore::MVDirection> &dirs,
               const casacore::Vector<casacore::uInt> &antIDs,
               const casacore::Vector<casacore::uInt> &feedIDs) const
{
  ASKAPDEBUGASSERT(antIDs.nelements() == feedIDs.nelements());
  const casacore::Vector<casacore::MVDirection> &directionCache =
      itsDirectionCache.value(*this,&TableConstDataIterator::fillDirectionCache);
  const casacore::Matrix<casacore::Int> &directionCacheIndices =
                 subtableInfo().getFeed().getIndices();
  dirs.resize(itsNumberOfRows);

  for (casacore::uInt row=0; row<itsNumberOfRows; ++row) {
       if ((feedIDs[row]>=directionCacheIndices.ncolumn()) ||
           (antIDs[row]>=directionCacheIndices.nrow())) {
              ASKAPTHROW(DataAccessError, "antID="<<antIDs[row]<<
                   " and/or feedID="<<feedIDs[row]<<
                   " are beyond the range of the FEED table");
           }
       if (directionCacheIndices(antIDs[row],feedIDs[row])<0) {
           ASKAPTHROW(DataAccessError, "The pair antID="<<antIDs[row]<<
                   " feedID="<<feedIDs[row]<<" doesn't have beam parameters defined");
       }
       const casacore::uInt index = static_cast<casacore::uInt>(
                    directionCacheIndices(antIDs[row],feedIDs[row]));
       ASKAPDEBUGASSERT(index < directionCache.nelements());
       dirs[row]=directionCache[index];
  }
}

/// @brief A helper method to fill a given vector with position angles.
/// @details fillFeedPA1 and fillFeedPA2 method do very similar operations,
/// which differ only by the feedIDs and antennaIDs used.
/// This method encapsulated these common operations
/// @param[in] angles a reference to vector to fill
/// @param[in] antIDs a vector with antenna IDs
/// @param[in] feedIDs a vector with feed IDs
/// @note There are some similarities between the code of this method and that
/// of fillVectorOfPointings. They are different with just a command called within
/// the loop. Theoretically, we can combine these two methods together, it would just
/// invovle some coding to make it look nice and probably some minor performance penalty.
void TableConstDataIterator::fillVectorOfPositionAngles(casacore::Vector<casacore::Float> &angles,
              const casacore::Vector<casacore::uInt> &antIDs,
              const casacore::Vector<casacore::uInt> &feedIDs) const
{
  ASKAPDEBUGASSERT(antIDs.nelements() == feedIDs.nelements());
  const casacore::Vector<casacore::Double> &parallacticAngles = itsParallacticAngleCache.value(*this,
                 &TableConstDataIterator::fillParallacticAngleCache);

  const casacore::MEpoch epoch=currentEpoch();
  ASKAPDEBUGASSERT(itsCurrentDataDescID>=0);
  const casacore::uInt spWindowID = currentSpWindowID();

  const IFeedSubtableHandler& feedSubtable = subtableInfo().getFeed();

  angles.resize(itsNumberOfRows);

  for (casacore::uInt row=0; row<itsNumberOfRows; ++row) {
       ASKAPDEBUGASSERT(antIDs[row] < parallacticAngles.nelements());
       angles[row] = casacore::Float(feedSubtable.getBeamPA(epoch, spWindowID, antIDs[row],
                               feedIDs[row]) + parallacticAngles[antIDs[row]]);
  }
}

/// @brief obtain the name of the data column
/// @details The visibility data can be taken not only from the DATA column,
/// but from any other appropriate column, e.g. CORRECTED_DATA. This method
/// returns the name of the column used to store such data. We need it in
/// derived classes to perform writing
/// @return name of the table column with visibility data
const std::string& TableConstDataIterator::getDataColumnName() const throw()
{
  ASKAPDEBUGASSERT(itsSelector);
  return itsSelector->getDataColumnName();
}

/// @brief obtain a current field ID
/// @details This method obtains a field ID corresponding to the
/// current iteration, if field ID column is present (and used). Otherwise
/// zero is always returned
/// @return current field ID
casacore::uInt TableConstDataIterator::currentFieldID() const
{
  return itsUseFieldID ? itsCurrentFieldID : 0u;
}

/// @brief obtain a current scan ID
/// @details This method obtains a scan number corresponding to the
/// current iteration. At this stage, this funcionality is not exposed via the
/// generic interface and is for use in test code only. In addition, there are
/// no measures taken to ensure that all rows of the iteration correspond to the
/// same scan ID (although realistically it should be the case because all
/// chunk corresponds to the same time stamp), although the MS standard allows it.
/// This method does the checks and throws an exception if scan number varies across
/// the chunk.
/// @return current scan ID
casacore::uInt TableConstDataIterator::currentScanID() const
{
  casacore::Vector<casacore::uInt> ids;
  fillVectorOfIDs(ids,"SCAN_NUMBER");
  ASKAPCHECK(ids.nelements()>0, "An attempt to extract scan ID for empty iteration");
  const casacore::uInt scanID = ids[0];
  // do cross-check
  for (casacore::uInt row=1; row<ids.nelements(); ++row) {
       ASKAPCHECK(scanID == ids[row], "Scan ID seem to differ for row="<<row<<" of the current iteration; was "<<scanID<<
                  " now "<<ids[row]);
  }
  return scanID;
}
