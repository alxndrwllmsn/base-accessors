/// @file TableConstDataIterator.h
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

#ifndef ASKAP_ACCESSORS_TABLE_CONST_DATA_ITERATOR_H
#define ASKAP_ACCESSORS_TABLE_CONST_DATA_ITERATOR_H

// std includes
#include <string>
#include <utility>

// boost includes
#include <boost/shared_ptr.hpp>

// casa includes
#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/Tables/TableIter.h>
#include <casacore/measures/Measures/Stokes.h>


// own includes
#include <askap/dataaccess/IConstDataIterator.h>
#include <askap/dataaccess/IDataConverterImpl.h>
#include <askap/dataaccess/ITableDataSelectorImpl.h>
#include <askap/dataaccess/TableConstDataAccessor.h>
#include <askap/dataaccess/TableInfoAccessor.h>
#include <askap/dataaccess/ITableManager.h>
#include <askap/dataaccess/CachedAccessorField.tcc>

namespace askap {

namespace accessors {

/// @brief Implementation of IConstDataIterator in the table-based case
/// @details
/// TableConstDataIterator: Allow read-only iteration across preselected data. Each
/// iteration step is represented by the IConstDataAccessor interface.
/// This is an implementation in the table-based case.
/// @ingroup dataaccess_tab
class TableConstDataIterator : virtual public IConstDataIterator,
                               virtual public TableInfoAccessor
{
public:
  /// @brief constructor of the const iterator
  /// @param[in] msManager a manager of the measurement set to use
  /// @param[in] sel shared pointer to selector
  /// @param[in] conv shared pointer to converter
  /// @param[in] cacheSize a number of uvw machines in the cache (default is 1)
  /// @param[in] tolerance pointing direction tolerance in radians, exceeding which leads
  /// to initialisation of a new UVW Machine
  /// @param[in] maxChunkSize maximum number of rows per accessor
  TableConstDataIterator(const boost::shared_ptr<ITableManager const>
              &msManager,
              const boost::shared_ptr<ITableDataSelectorImpl const> &sel,
	      const boost::shared_ptr<IDataConverterImpl const> &conv,
	      size_t cacheSize = 1, double tolerance = 1e-6,
	      casacore::uInt maxChunkSize = INT_MAX);

  /// Restart the iteration from the beginning
  virtual void init();

  /// operator* delivers a reference to data accessor (current chunk)
  /// @return a reference to the current chunk
  virtual const IConstDataAccessor& operator*() const;

  /// Checks whether there are more data available.
  /// @return True if there are more data available
  virtual casacore::Bool hasMore() const throw();

  /// advance the iterator one step further
  /// @return True if there are more data (so constructions like
  ///         while(it.next()) {} are possible)
  virtual casacore::Bool next();

  /// methods used in the accessor.

  /// @return number of rows in the current accessor
  casacore::uInt inline nRow() const throw() { return itsNumberOfRows;}

  /// temporary - access to number of channels and polarizations
  /// it will be determined by the selector, get them from the table
  /// for now

  /// @return number of channels in the current accessor
  casacore::uInt inline nChannel() const throw() { return getChannelRange().first;}

  /// @return number of channels in the current accessor
  casacore::uInt inline nPol() const throw() { return itsNumberOfPols;}

  /// populate the buffer of visibilities with the values of current
  /// iteration
  /// @param[in] vis a reference to the nRow x nChannel x nPol buffer
  ///            cube to fill with the complex visibility data
  void fillVisibility(casacore::Cube<casacore::Complex> &vis) const;

  /// populate the buffer of noise figures with the values of current
  /// iteration
  /// @param[in] noise a reference to the nRow x nChannel x nPol buffer
  ///            cube to be filled with the noise figures
  void fillNoise(casacore::Cube<casacore::Complex> &noise) const;

  /// @brief read flagging information
  /// @details populate the buffer of flags with the information
  /// read in the current iteration
  /// @param[in] flag a reference to the nRow x nChannel x nPol buffer
  ///            cube to fill with the flag information (each element has
  ///            bool type)
  void fillFlag(casacore::Cube<casacore::Bool> &flag) const;

  /// populate the buffer with uvw
  /// @param[in] uvw a reference to vector of rigid vectors (3 elemets,
  ///            u,v and w for each row) to fill
  void fillUVW(casacore::Vector<casacore::RigidVector<casacore::Double, 3> >&uvw) const;

  /// populate the buffer with frequencies
  /// @param[in] freq a reference to a vector to fill
  void fillFrequency(casacore::Vector<casacore::Double> &freq) const;

  /// @return the time stamp in the table's native frame/units
  /// @note this method doesn't do any caching. It reads the table each
  /// time it is called. It is intended for use from the accessor only, where
  /// caching is done.
  casacore::Double getTime() const;

  /// @brief an alternative way to get the time stamp
  /// @details This method uses the accessor to get cached time stamp. It
  /// is returned as an epoch measure.
  casacore::MEpoch currentEpoch() const;

  /// populate the buffer with IDs of the first antenna
  /// @param[in] ids a reference to a vector to fill
  void fillAntenna1(casacore::Vector<casacore::uInt> &ids) const;

  /// populate the buffer with IDs of the second antenna
  /// @param[in] ids a reference to a vector to fill
  void fillAntenna2(casacore::Vector<casacore::uInt> &ids) const;

  /// populate the buffer with IDs of the first feed
  /// @param[in] ids a reference to a vector to fill
  void fillFeed1(casacore::Vector<casacore::uInt> &ids) const;

  /// populate the buffer with IDs of the second feed
  /// @param[in] ids a reference to a vector to fill
  void fillFeed2(casacore::Vector<casacore::uInt> &ids) const;

  /// fill the buffer with the pointing directions of the first antenna/feed
  /// @param[in] dirs a reference to a vector to fill
  void fillPointingDir1(casacore::Vector<casacore::MVDirection> &dirs) const;

  /// fill the buffer with the pointing directions of the second antenna/feed
  /// @param[in] dirs a reference to a vector to fill
  void fillPointingDir2(casacore::Vector<casacore::MVDirection> &dirs) const;

  /// fill the buffer with the position angles of the first antenna/feed
  /// @param[in] angles a reference to vector to be filled
  void fillFeed1PA(casacore::Vector<casacore::Float> &angles) const;

  /// fill the buffer with the position angles of the second antenna/feed
  /// @param[in] angles a reference to vector to be filled
  void fillFeed2PA(casacore::Vector<casacore::Float> &angles) const;

  /// @brief fill the buffer with the pointing directions for the first antenna centre
  /// @details The difference from fillPointingDir1 is that no feed offset is applied.
  /// @param[in] dirs a reference to a vector to fill
  void fillDishPointing1(casacore::Vector<casacore::MVDirection> &dirs) const;

  /// @brief fill the buffer with the pointing directions for the second antenna centre
  /// @details The difference from fillPointingDir2 is that no feed offset is applied.
  /// @param[in] dirs a reference to a vector to fill
  void fillDishPointing2(casacore::Vector<casacore::MVDirection> &dirs) const;

  /// @brief fill the buffer with the polarisation types
  /// @param[in] stokes a reference to a vector to be filled
  void fillStokes(casacore::Vector<casacore::Stokes::StokesTypes> &stokes) const;

  /// @brief UVW machine cache size
  /// @return size of the uvw machine cache
  inline size_t uvwMachineCacheSize() const {return itsUVWCacheSize;}

  /// @brief direction tolerance used for UVW machine cache
  /// @return direction tolerance used for UVW machine cache (in radians)
  inline double uvwMachineCacheTolerance() const {return itsUVWCacheTolerance;}

  /// @brief obtain a current field ID
  /// @details This method obtains a field ID corresponding to the
  /// current iteration, if field ID column is present (and used). Otherwise
  /// zero is always returned
  /// @return current field ID
  casacore::uInt currentFieldID() const;

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
  casacore::uInt currentScanID() const;

protected:
  /// @brief obtain selected range of channels
  /// @details A subset spectral channels can be selected for this iterator to work with.
  /// This method returns the number of channels and the first selected channel.
  /// @return a pair, first element is the number of channels, second is the first channel
  /// in the full cube
  std::pair<casacore::uInt, casacore::uInt> getChannelRange() const;

  /// @brief a short cut to get the first channel in the full cube
  /// @return the number of the first channel in the full cube
  inline casacore::uInt startChannel() const { return getChannelRange().second;}

  /// @brief read an array column of the table into a cube
  /// @details populate the buffer provided with the information
  /// read in the current iteration. This method is templated and can be
  /// used for both visibility and flag data fillers.
  /// @param[in] cube a reference to the nRow x nChannel x nPol buffer
  ///            cube to fill with the information from table
  /// @param[in] columnName a name of the column to read
  template<typename T>
  void fillCube(casacore::Cube<T> &cube, const std::string &columnName) const;

  /// @brief A helper method to fill a given vector with pointing directions.
  /// @details fillPointingDir1 and fillPointingDir2 methods do very similar
  /// operations, which differ only by the feedIDs and antennaIDs used.
  /// This method encapsulates these common operations
  /// @param[in] dirs a reference to vector to fill
  /// @param[in] antIDs a vector with antenna IDs
  /// @param[in] feedIDs a vector with feed IDs
  void fillVectorOfPointings(casacore::Vector<casacore::MVDirection> &dirs,
               const casacore::Vector<casacore::uInt> &antIDs,
               const casacore::Vector<casacore::uInt> &feedIDs) const;

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
  void fillVectorOfPositionAngles(casacore::Vector<casacore::Float> &angles,
                const casacore::Vector<casacore::uInt> &antIDs,
                const casacore::Vector<casacore::uInt> &feedIDs) const;

  /// @brief a helper method to get dish pointings
  /// @details fillDishPointing1 and fillDishPointing2 methods do very
  /// similar operations, which differ only by the antennaIDs used.
  /// This method encapsulated these common operations.
  /// @note fillVectorOfPointings computes pointing directions for
  /// individual feeds, not for the centre of the dish as this method
  /// does
  /// @param[in] dirs a reference to a vector to fill
  /// @param[in] antIDs a vector with antenna IDs
  void fillVectorOfDishPointings(casacore::Vector<casacore::MVDirection> &dirs,
               const casacore::Vector<casacore::uInt> &antIDs) const;

  /// @brief a helper method to read a column with IDs of some sort
  /// @details It reads the column of casacore::Int and fills a Vector of
  /// casacore::uInt. A check to ensure all numbers are non-negative is done
  /// in the debug mode.
  /// @param[in] ids a reference to a vector to fill
  /// @param[in] name a name of the column to read
  void fillVectorOfIDs(casacore::Vector<casacore::uInt> &ids,
                       const casacore::String &name) const;

  /// setup accessor for a new iteration
  void setUpIteration();

  /// @brief method ensures that the chunk has a uniform DATA_DESC_ID
  /// @details This method reduces itsNumberOfRows to achieve
  /// uniform DATA_DESC_ID reading for all rows in the current chunk.
  /// The resulting itsNumberOfRows will be 1 or more.
  /// itsAccessor's spectral axis cache is reset if new DATA_DESC_ID is
  /// different from itsCurrentDataDescID
  /// This method also sets up itsNumberOfPols and itsNumberOfChannels
  /// when DATA_DESC_ID changes (and therefore at the first run as well)
  void makeUniformDataDescID();

  /// @brief method ensures that the chunk has a uniform FIELD_ID
  /// @details This method reduces itsNumberOfRows until FIELD_ID is
  /// the same for all rows in the current chunk. The resulting
  /// itsNumberOfRows will be 1 or more. If itsUseFieldID is false,
  /// the method returns without doing anything. itsAccessor's direction
  /// cache is reset if new FIELD_ID is different from itsCurrentFieldID
  /// (and it sets it up at the first run as well)
  void makeUniformFieldID();

  /// obtain a reference to the accessor (for derived classes)
  inline const TableConstDataAccessor& getAccessor() const throw()
  { return itsAccessor;}

  /// @brief Fill an internal buffer with the pointing directions
  /// @details  The layout of this buffer is the same as the layout of
  /// the FEED subtable for current time and spectral window.
  /// getAntennaIDs and getFeedIDs methods of the
  /// subtable handler can be used to unwrap this 1D array.
  /// The buffer can be invalidated if the time changes (i.e. for an alt-az array),
  /// for an equatorial array this happends only if the FEED or FIELD subtable
  /// are time-dependent (i.e. when the pointing changes)
  /// @param[in] dirs a reference to a vector to be filled
  void fillDirectionCache(casacore::Vector<casacore::MVDirection> &dirs) const;

  /// @brief Fill internal buffer with parallactic angles
  /// @details This buffer holds parallactic angles for all antennae. The buffer
  /// is invalidated when the time changes for an alt-az array, for an equatorial
  /// array it happens only if the pointing changes.
  /// @param[in] angles a reference to a vector to be filled
  void fillParallacticAngleCache(casacore::Vector<casacore::Double> &angles) const;

  /// @brief fill the buffer with the dish pointing directions
  /// @details The difference from fillDirectionCache is that
  /// this method computes the pointing directions for the dish centre, not for
  /// individual feeds (or synthetic beams, strictly speaking). The number of elements
  /// in the buffer equals to the number of antennae. This is also different from
  /// fillDirectionCache, which projects feeds to the same 1D array as well.
  /// @note At this stage we use FIELD subtable to get the pointing directions.
  /// Therefore, these directions do not depend on antenna/feed. This method writes
  /// the same value for all elements of the array. It will be used for both antennae
  /// in the pair.
  /// @param[in] dirs a reference to a vector to fill
  void fillDishPointingCache(casacore::Vector<casacore::MVDirection> &dirs) const;

  /// @brief obtain a current spectral window ID
  /// @details This method obtains a spectral window ID corresponding to the
  /// current data description ID and tests its validity
  /// @return current spectral window ID
  casacore::uInt currentSpWindowID() const;

  /// @brief obtain a current polarisation ID
  /// @details This method obtains a polarisation ID corresponding to the current
  /// data description ID and tests its validity
  /// @return current polarisation ID
  casacore::uInt currentPolID() const;

  /// @brief obtain the current iteration of the table iterator
  /// @details This class uses TableIterator behind the scene. This method
  /// returns the current iteration, which can be used in derived classes
  /// (e.g. for read-write access)
  /// @return a const reference to table object representing the current iteration
  inline const casacore::Table& getCurrentIteration() const throw()
       {return itsCurrentIteration;}

  /// @brief obtain the current top row
  /// @details This class uses TableIterator behind the scence. One iteration
  /// of the table iterator may cover more than one iteration of the iterator
  /// represented by this class. The result of this method is a row number,
  /// where current data accessor starts.
  /// @return row number in itsCurrentItrertion corresponding to row 0 of the
  /// accessor at this iteration
  inline casacore::rownr_t getCurrentTopRow() const throw() {return itsCurrentTopRow;}

  /// @brief obtain the name of the data column
  /// @details The visibility data can be taken not only from the DATA column,
  /// but from any other appropriate column, e.g. CORRECTED_DATA. This method
  /// returns the name of the column used to store such data. We need it in
  /// derived classes to perform writing
  /// @return name of the table column with visibility data
  const std::string& getDataColumnName() const throw();

  /// @brief obtain a reference direction for the current iteration
  /// @details Currently we assume that the dish pointing centre stays
  /// fixed for the whole chunk. We break the iteration, if necessary
  /// to achieve this. This helper method extracts the reference direction
  /// from the FIELD subtable using either FIELD_ID, or current time if
  /// the former is not supported by the main table.
  /// @return a reference to direction measure
  const casacore::MDirection& getCurrentReferenceDir() const;

private:
  // note, it is essential that itsUVWCacheSize and itsUVWCacheTolerance are initialised
  // prior to itsAccessor (accessor uses them in its setup)

  /// @brief a number of uvw machines in the cache (default is 1)
  /// @details To speed up mosaicing it is possible to cache any number of uvw machines
  /// as it takes time to setup the transformation which depends on the phase centre.
  /// A change to this parameter applies to all iterators created afterwards.
  size_t itsUVWCacheSize;

  /// @brief pointing direction tolerance in radians (for uvw machine cache)
  /// @details Exceeding this tolerance leads to initialisation of a new UVW Machine in the cache
  double itsUVWCacheTolerance;

  /// accessor (a chunk of data)
  /// although the accessor type can be different
  TableConstDataAccessor itsAccessor;

  boost::shared_ptr<ITableDataSelectorImpl const>  itsSelector;
  boost::shared_ptr<IDataConverterImpl>  itsConverter;
  /// the maximum allowed number of rows in the accessor.
  casacore::uInt itsMaxChunkSize;
  casacore::TableIterator itsTabIterator;
  /// current group of data returned by itsTabIterator
  casacore::Table itsCurrentIteration;
  /// current row in the itsCurrentIteration projected to the row 0
  /// of the data accessor
  casacore::rownr_t itsCurrentTopRow;
  /// number of rows in the current chunk
  casacore::uInt itsNumberOfRows;
  /// next two data members show the number of channels and
  /// polarisations in the actual table. Selector controls what
  /// is sent out
  casacore::uInt itsNumberOfChannels;
  /// see above
  casacore::uInt itsNumberOfPols;

  /// current DATA_DESC_ID, the iteration is broken if this
  /// ID changes
  casacore::Int itsCurrentDataDescID;

  /// @brief current FIELD_ID.
  /// @details This ID is tracked if FIELD_ID column is present in the
  /// table. The iteration is broken if this ID changes.
  casacore::Int itsCurrentFieldID;

  /// @brief a flag showing that FIELD_ID column should be used.
  /// @details There are two ways to discriminate between different pointings:
  /// use FIELD_ID column, if it is present in the table, and check times. If
  /// this flag is set, the iterator will check FIELD_ID column (an exception is
  /// raised if the column doesn't exist). The constructor checks the presence
  /// of the FIELD_ID column and set this flag if it exists. The flag is introduced
  /// to allow in the future to force the code to use time instead of FIELD_ID, even
  /// if the latter is present.
  bool itsUseFieldID;

  /// @brief cache of pointing directions  for each feed
  /// @details This is an internal buffer for pointing
  /// directions for the whole current cache of the Feed subtable handler
  CachedAccessorField<casacore::Vector<casacore::MVDirection> > itsDirectionCache;

  /// @brief cache of parallactic angles for each antenna
  /// @details This is an internal buffer for parallactic angles (in radians) for the whole
  /// current cache of the antenna subtable handler
  CachedAccessorField<casacore::Vector<casacore::Double> > itsParallacticAngleCache;

  /// internal buffer for dish pointings for all antennae
  CachedAccessorField<casacore::Vector<casacore::MVDirection> > itsDishPointingCache;

  /// currently selected number of channels
  mutable uint itsNumberOfChannelsSelected;
  /// currently selected start channel
  mutable uint itsStartChannelSelected;
  /// selection initialised?
  mutable bool itsChannelsSelected;
  /// selection invalid?
  mutable bool itsFlagData;
  /// are we at the start?
  mutable bool itsAtStart;
};


} // end of namespace accessors

} // end of namespace askap
#endif // #ifndef TABLE_CONST_DATA_ITERATOR_H
