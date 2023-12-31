/// @file TableConstDataSource.cc
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

/// boost includes
#include <boost/shared_ptr.hpp>

/// own includes
#include <askap/dataaccess/TableConstDataSource.h>
#include <askap/dataaccess/TableConstDataIterator.h>
#include <askap/dataaccess/TableDataSelector.h>
#include <askap/dataaccess/BasicDataConverter.h>
#include <askap/dataaccess/DataAccessError.h>

using namespace askap;
using namespace askap::accessors;
using namespace casa;

/// @brief construct a read-only data source object
/// @details All iterators obtained from this object will be read-only
/// iterators.
/// @param[in] fname file name of the measurement set to use
/// @param[in] dataColumn a name of the data column used by default
///                       (default is DATA)
TableConstDataSource::TableConstDataSource(const std::string &fname,
               const std::string &dataColumn) :
         TableInfoAccessor(casacore::Table(fname), false, dataColumn),
         itsUVWCacheSize(1), itsUVWCacheTolerance(1e-6),
         itsMaxChunkSize(INT_MAX) {}

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
const casacore::MPosition& TableConstDataSource::getAntennaPosition(casacore::uInt antID) const
{
  // the validity of indicies and initialisation of caches is checked inside these methods, 
  // but only in the debug mode
  return subtableInfo().getAntenna().getPosition(antID);
}

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
casacore::uInt TableConstDataSource::getNumberOfAntennas() const
{
  // the validity of indicies and initialisation of caches is checked inside these methods, 
  // but only in the debug mode
  return subtableInfo().getAntenna().getNumberOfAntennas();
}

/// @brief configure restriction on the chunk size
/// @param[in] maxNumRows maximum number of rows wanted
/// @note The new restriction will apply to any iterator created in the future, but will not
/// affect iterators already created
void TableConstDataSource::configureMaxChunkSize(casacore::uInt maxNumRows)
{
   ASKAPCHECK(maxNumRows > 0, "Restriction on the number of rows should be a positive number");
   itsMaxChunkSize = maxNumRows;
}

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
void TableConstDataSource::configureUVWMachineCache(size_t cacheSize, double tolerance)
{
  itsUVWCacheSize = cacheSize;
  itsUVWCacheTolerance = tolerance;
}

/// construct a part of the read only object for use in the
/// derived classes
/// @note Due to virtual inheritance, TableInfoAccessor will be initialized
/// in the concrete derived class. This empty constructor is added to make
/// the compiler happy
TableConstDataSource::TableConstDataSource() :
         TableInfoAccessor(boost::shared_ptr<ITableManager const>()),
         itsUVWCacheSize(1), itsUVWCacheTolerance(1e-6),
         itsMaxChunkSize(INT_MAX) {} 

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
IDataConverterPtr TableConstDataSource::createConverter() const
{
  return IDataConverterPtr(new BasicDataConverter);
}

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
boost::shared_ptr<IConstDataIterator>
TableConstDataSource::createConstIterator(const IDataSelectorConstPtr &sel,
              const IDataConverterConstPtr &conv) const
{
   // cast input selector to "implementation" interface
   boost::shared_ptr<ITableDataSelectorImpl const> implSel=
           boost::dynamic_pointer_cast<ITableDataSelectorImpl const>(sel);
   boost::shared_ptr<IDataConverterImpl const> implConv=
           boost::dynamic_pointer_cast<IDataConverterImpl const>(conv);
   	   
   if (!implSel || !implConv) {
       ASKAPTHROW(DataAccessLogicError, "Incompatible selector and/or "<<
                 "converter are received by the createConstIterator method");
   }
   return boost::shared_ptr<IConstDataIterator>(new TableConstDataIterator(
                getTableManager(),implSel,implConv,uvwMachineCacheSize(), uvwMachineCacheTolerance(),
                maxChunkSize()));
}

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
IDataSelectorPtr TableConstDataSource::createSelector() const
{
  return IDataSelectorPtr(new TableDataSelector(getTableManager()));
}

