/// @file TableConstDataSource.h
/// @brief Implementation of IConstDataSource in the table-based case
/// @details
/// TableConstDataSource: Allow read-only access to the data stored in the
/// measurement set. This class implements IConstDataSource interface.
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

#ifndef ASKAP_ACCESSORS_TABLE_CONST_DATA_SOURCE_H
#define ASKAP_ACCESSORS_TABLE_CONST_DATA_SOURCE_H

// boost includes
#include <boost/shared_ptr.hpp>

// casa includes
#include <casacore/tables/Tables/Table.h>

// own includes
#include <askap/dataaccess/IConstDataSource.h>
#include <askap/dataaccess/TableInfoAccessor.h>

// std includes
#include <string>

namespace askap {

namespace accessors {

/// @brief Implementation of IConstDataSource in the table-based case
/// @details
/// TableConstDataSource: Allow read-only access to the data stored in the
/// measurement set. This class implements IConstDataSource interface.
/// @ingroup dataaccess_tab
class TableConstDataSource : virtual public IConstDataSource,
                             virtual protected TableInfoAccessor
{
public:
  /// @brief construct a read-only data source object
  /// @details All iterators obtained from this object will be read-only
  /// iterators.
  /// @param[in] fname file name of the measurement set to use
  /// @param[in] dataColumn a name of the data column used by default
  ///                       (default is DATA)
  explicit TableConstDataSource(const std::string &fname, 
                                const std::string &dataColumn = "DATA");
  
  /// create a converter object corresponding to this type of the
  /// DataSource. The user can change converting policies (units,
  /// reference frames) by appropriate calls to this converter object
  /// and pass it back to createConstIterator(...). The data returned by
  /// the iteratsr will automatically be in the requested frame/units
  ///
  /// @return a shared pointer to a new DataConverter object
  ///
  /// The method acts as a factory by creating a new DataConverter.
  /// The lifetime of this converter is the same as the lifetime of the
  /// DataSource object. Therefore, it can be reused multiple times,
  /// if necessary. However, the behavior of iterators created
  /// with a particular DataConverter is undefined, if you change
  /// the DataConverter after the creation of an iterator, unless you
  /// call init() of the iterator (and start a new iteration loop).
  virtual IDataConverterPtr createConverter() const;
  
  /// get iterator over a selected part of the dataset represented
  /// by this DataSource object with an explicitly specified conversion
  /// policy. This is the most general createConstIterator(...) call, 
  /// which is used as a default implementation for all less general 
  /// cases (although they can be overriden in the derived classes, if it 
  /// will be necessary because of the performance issues)
  ///
  /// @param[in] sel a shared pointer to the selector object defining 
  ///            which subset of the data is used
  /// @param[in] conv a shared pointer to the converter object defining
  ///            reference frames and units to be used
  /// @return a shared pointer to DataIterator object
  ///
  /// The method acts as a factory by creating a new DataIterator.
  /// The lifetime of this iterator is the same as the lifetime of
  /// the DataSource object. Therefore, it can be reused multiple times,
  /// if necessary. Call init() to rewind the iterator.
  virtual boost::shared_ptr<IConstDataIterator> createConstIterator(const
             IDataSelectorConstPtr &sel,
             const IDataConverterConstPtr &conv) const;
  
  // we need this to get access to the overloaded syntax in the base class 
  using IConstDataSource::createConstIterator;
 
  /// create a selector object corresponding to this type of the
  /// DataSource
  ///
  /// @return a shared pointer to the DataSelector corresponding to
  /// this type of DataSource. DataSource acts as a factory and
  /// creates a selector object of the appropriate type
  ///
  /// This method acts as a factory by creating a new DataSelector
  /// appropriate to the given DataSource. The lifetime of the
  /// DataSelector is the same as the lifetime of the DataSource 
  /// object. Therefore, it can be reused multiple times,
  /// if necessary. However, the behavior of iterators already obtained
  /// with this DataSelector is undefined, if one changes the selection
  /// unless the init method is called for the iterator (and the new
  /// iteration loop is started).
  virtual IDataSelectorPtr createSelector() const;
  
  /// @brief configure caching of the uvw-machines
  /// @details A number of uvw machines can be cached at the same time. This can
  /// result in a significant performance improvement in the mosaicing case. By default
  /// only single machine is cached and this method should be called to change it. 
  /// All subsequent iterators will be created with the parameters set in this method until
  /// it is called again. Call this method without parameters to revert to default settings.
  /// @note This method is a feature of this implementation and is not available via the 
  /// general interface (intentionally)
  /// @param[in] cacheSize a number of uvw machines in the cache (default is 1)
  /// @param[in] tolerance pointing direction tolerance in radians, exceeding which leads 
  /// to initialisation of a new UVW Machine
  void configureUVWMachineCache(size_t cacheSize = 1, double tolerance = 1e-6);

  /// @brief configure restriction on the chunk size
  /// @param[in] maxNumRows maximum number of rows wanted
  /// @note The new restriction will apply to any iterator created in the future, but will not
  /// affect iterators already created
  void configureMaxChunkSize(casacore::uInt maxNumRows);

  /// @brief obtain the position of the given antenna
  /// @details
  /// @param[in] antID antenna index to use, matches indices in the data table
  /// @return const reference to MPosition measure for the chosen antenna
  /// @note This method is deliberately not exposed via the IConstDataSource interface because
  /// it is table-specific and cannot be implemented in general in the streaming model where such metadata
  /// should be provided some other way (i.e. not in the stream). However, in the table-based case it can
  /// be used directly as the type is created explicitly at some point (or one could dynamic cast and
  /// test whether this operation is supported). Same information can be extracted manually via the 
  /// getTableManager() method of table-based iterators, essentially the same code used in this
  /// shortcut, but this method is public for iterators.
  const casacore::MPosition& getAntennaPosition(casacore::uInt antID) const;

  /// @brief obtain the number of antennas
  /// @details This is another method specific to the table-based implementation (in the streaming approach
  /// this has to be provided some other way, through configuration). Therefore, similarly to getAntennaPosition,
  /// it is not exposed via the IConstDataSource interface making its use more explicit in the code. In principle,
  /// the number of antennas should rarely be needed in the user code as only valid indices are returned by the 
  /// the accessor. 
  /// @note Strictly speaking, this is not the number of antennas, in general, but the number of entries in the 
  /// ANTENNA table of the measurement set which may not match (and the actual data may only use a subset of indices -
  /// this is yet another indication that ideally the user-level code should avoid this implementation-specific information).
  /// @return number of antennas in the measurement set (all antenna indices are less than this number)
  casacore::uInt getNumberOfAntennas() const;
  
protected:
  /// construct a part of the read only object for use in the
  /// derived classes
  TableConstDataSource();
  
  /// @brief UVW machine cache size
  /// @return size of the uvw machine cache
  inline size_t uvwMachineCacheSize() const {return itsUVWCacheSize;}
  
  /// @brief direction tolerance used for UVW machine cache
  /// @return direction tolerance used for UVW machine cache (in radians)
  inline double uvwMachineCacheTolerance() const {return itsUVWCacheTolerance;}   

  /// @brief current restriction on the chunk size
  /// @return maximum number of rows in the accessor (the current setting, affects future iterators)
  inline casacore::uInt maxChunkSize() const {return itsMaxChunkSize;}
  
private:
  /// @brief a number of uvw machines in the cache (default is 1)
  /// @details To speed up mosaicing it is possible to cache any number of uvw machines
  /// as it takes time to setup the transformation which depends on the phase centre. 
  /// A change to this parameter applies to all iterators created afterwards. 
  size_t itsUVWCacheSize;
  
  /// @brief pointing direction tolerance in radians (for uvw machine cache)
  /// @details Exceeding this tolerance leads to initialisation of a new UVW Machine in the cache
  double itsUVWCacheTolerance;


  /// @brief maximum number of rows per accessor
  /// @details By default, it is initialised with INT_MAX, which essentially means no restrictions.
  /// However, the maximum number of rows can be constrained to some value. This would provide more
  /// iterations but with smaller chunks. It can make sense if there are multiple adapters in the
  /// processing chain which do data copy (usually in the temporary code/hacks which technically shouldn't
  /// stay long term in the ideal case).
  casacore::uInt itsMaxChunkSize;
};
 
} // namespace accessors

} // namespace askap

#endif // #ifndef TABLE_CONST_DATA_SOURCE_H
