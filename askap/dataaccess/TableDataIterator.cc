/// @file
///
/// @brief Implementation of IDataIterator in the table-based case
/// @details
/// TableConstDataIterator: Allow read-only iteration across preselected data. Each
/// iteration step is represented by the IConstDataAccessor interface.
/// TableDataIterator extends the interface further to read-write operations.
/// Each iteration step is represented by the IDataAccessor interface in this
/// case.
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

// stl includes
#include <algorithm>
#include <functional>
#include <utility>

// own includes
#include <askap/dataaccess/TableDataIterator.h>
#include <askap/dataaccess/TableBufferDataAccessor.h>
#include <askap/dataaccess/TableDataAccessor.h>
#include <askap/dataaccess/TableInfoAccessor.h>
#include <askap/dataaccess/IBufferManager.h>
#include <askap/dataaccess/DataAccessError.h>

// casa includes
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>


namespace askap {

namespace accessors {

/// @brief an adapter to use stl functions with maps.
/// @details
/// STL provides adapters to use pointers to member functions as an object
/// function. However, they won't work with maps, because map has an iterator
/// returning a pair rather than the value itself. This adapter implements the
/// same thing as one of the mem_fun STL adapters, but for maps.
/// It can be moved in its own file, if found useful in other parts of the code
/// @ingroup dataaccess_hlp
template<typename X>
struct MapMemFun : public std::unary_function<X*, void> {
  /// construct adapter for member function in
  /// @param[in] in member function
  explicit MapMemFun(void (X::*in)()) : func(in) {}

  /// @param[in] in a reference to the pair-like object with a pointer-like
  /// object stored in the second parameter.
  template<typename P>
  void operator()(const P &in) const
    {
      ASKAPDEBUGASSERT(in.second);
      return ((*in.second).*func)();
    }
private:
  /// @brief buffered pointer to a member function
  void (X::*func)();
};

/// a helper method to autodetect the input type
/// @param[in] in a member function to call for all elements in the map
template<typename X>
MapMemFun<X> mapMemFun(void (X::*in)()) {
  return MapMemFun<X>(in);
}

}
}

using namespace askap;
using namespace askap::accessors;

/// @param[in] msManager a manager of the measurement set to use
/// @param[in] sel shared pointer to selector
/// @param[in] conv shared pointer to converter
/// @param[in] maxChunkSize maximum number of rows per accessor
TableDataIterator::TableDataIterator(
            const boost::shared_ptr<ITableManager const> &msManager,
            const boost::shared_ptr<ITableDataSelectorImpl const> &sel,
            const boost::shared_ptr<IDataConverterImpl const> &conv,
            size_t cacheSize, double tolerance,
            casacore::uInt maxChunkSize) :
         TableInfoAccessor(msManager),
           TableConstDataIterator(msManager,sel,conv,cacheSize, tolerance, maxChunkSize),
	      itsOriginalVisAccessor(new TableDataAccessor(*this)),
	      itsIterationCounter(0)
{
  itsActiveBufferPtr=itsOriginalVisAccessor;
}

/// @brief operator* delivers a reference to data accessor (current chunk)
/// @details
/// @return a reference to the current chunk
/// @note
/// constness of the return type is changed to allow read/write
/// operations.
///
IDataAccessor& TableDataIterator::operator*() const
{
  ASKAPDEBUGASSERT(itsActiveBufferPtr);
  return *itsActiveBufferPtr;
}

/// @brief Switch the output of operator* and operator-> to one of
/// the buffers.
/// @details This is meant to be done to provide the same
/// interface for a buffer access as exists for the original
/// visibilities (e.g. it->visibility() to get the cube).
/// It can be used for an easy substitution of the original
/// visibilities to ones stored in a buffer, when the iterator is
/// passed as a parameter to mathematical algorithms.
/// The operator* and operator-> will refer to the chosen buffer
/// until a new buffer is selected or the chooseOriginal() method
/// is executed to revert operators to their default meaning
/// (to refer to the primary visibility data).
/// @param[in] bufferID  the name of the buffer to choose
void TableDataIterator::chooseBuffer(const std::string &bufferID)
{
  std::map<std::string,
      boost::shared_ptr<TableBufferDataAccessor> >::const_iterator bufferIt =
                      itsBuffers.find(bufferID);

  if (bufferIt==itsBuffers.end()) {
      // deal with new buffer
      itsActiveBufferPtr = itsBuffers[bufferID] =
            boost::shared_ptr<TableBufferDataAccessor>(
	          new TableBufferDataAccessor(bufferID,*this));
  } else {
      itsActiveBufferPtr=bufferIt->second;
  }
}

/// Switch the output of operator* and operator-> to the original
/// state (present after the iterator is just constructed)
/// where they point to the primary visibility data. This method
/// is indended to cancel the results of chooseBuffer(std::string)
///
void TableDataIterator::chooseOriginal()
{
  itsActiveBufferPtr=itsOriginalVisAccessor;
}

/// @brief obtain any associated buffer for read/write access.
/// @details The buffer is identified by its bufferID. The method
/// ignores a chooseBuffer/chooseOriginal setting.
/// @param[in] bufferID the name of the buffer requested
/// @return a reference to writable data accessor to the
///         buffer requested
IDataAccessor& TableDataIterator::buffer(const std::string &bufferID) const
{
  std::map<std::string,
      boost::shared_ptr<TableBufferDataAccessor> >::const_iterator bufferIt =
                      itsBuffers.find(bufferID);
  if (bufferIt!=itsBuffers.end()) {
      // this buffer already exists
      return *(bufferIt->second);
  }
  // this is a request for a new buffer
  return *(itsBuffers[bufferID] =
      boost::shared_ptr<TableBufferDataAccessor>(new TableBufferDataAccessor(
                                       bufferID,*this)));
}

/// Restart the iteration from the beginning
void TableDataIterator::init()
{
  // call sync() member function for all accessors in itsBuffers
  std::for_each(itsBuffers.begin(),itsBuffers.end(),
           mapMemFun(&TableBufferDataAccessor::sync));
  ASKAPDEBUGASSERT(itsOriginalVisAccessor);
  itsOriginalVisAccessor->sync();

  TableConstDataIterator::init();
  itsIterationCounter=0;

  // call notifyNewIteration() member function for all accessors
  // in itsBuffers
  std::for_each(itsBuffers.begin(),itsBuffers.end(),
           mapMemFun(&TableBufferDataAccessor::notifyNewIteration));
  // original visibilities will be read on-demand by the code in
  // TableConstDataAccessor in a usual way
}


/// advance the iterator one step further
/// @return True if there are more data (so constructions like
///         while(it.next()) {} are possible)
casacore::Bool TableDataIterator::next()
{
  // call sync() member function for all accessors in itsBuffers
  std::for_each(itsBuffers.begin(),itsBuffers.end(),
           mapMemFun(&TableBufferDataAccessor::sync));
  ASKAPDEBUGASSERT(itsOriginalVisAccessor);
  itsOriginalVisAccessor->sync();

  ++itsIterationCounter;

  // call notifyNewIteration() member function for all accessors
  // in itsBuffers
  std::for_each(itsBuffers.begin(),itsBuffers.end(),
           mapMemFun(&TableBufferDataAccessor::notifyNewIteration));
  // original visibilities will be read on-demand by the code in
  // TableConstDataAccessor in a usual way

  return TableConstDataIterator::next();
}

/// populate the cube with the data stored in the given buffer
/// @param[in] vis a reference to the nRow x nChannel x nPol buffer
///            cube to fill with the complex visibility data
/// @param[in] name a name of the buffer to work with
void TableDataIterator::readBuffer(casacore::Cube<casacore::Complex> &vis,
                        const std::string &name) const
{
  const IBufferManager &bufManager=subtableInfo().getBufferManager();
  const TableConstDataAccessor &accessor=getAccessor();
  const casacore::IPosition requiredShape(3, accessor.nRow(),
          accessor.nChannel(), accessor.nPol());
  if (bufManager.bufferExists(name,itsIterationCounter)) {
      bufManager.readBuffer(vis,name,itsIterationCounter);
      if (vis.shape()!=requiredShape) {
	  // this is an old buffer with a different shape. Can't be used
	  vis.resize(requiredShape);
      }
  } else {
      vis.resize(requiredShape);
  }
}

/// write the cube back to the given buffer
/// @param[in] vis a reference to the nRow x nChannel x nPol buffer
///            cube to fill with the complex visibility data
/// @param[in] name a name of the buffer to work with
void TableDataIterator::writeBuffer(const casacore::Cube<casacore::Complex> &vis,
                         const std::string &name) const
{
  subtableInfo().getBufferManager().writeBuffer(vis,name,itsIterationCounter);
}

/// destructor required to sync buffers on the last iteration
TableDataIterator::~TableDataIterator()
{
  // call sync() member function for all accessors in itsBuffers
  std::for_each(itsBuffers.begin(),itsBuffers.end(),
           mapMemFun(&TableBufferDataAccessor::sync));
  if (itsOriginalVisAccessor) {
      // there is no much point to throw an exception here if itsOriginalVisAccesor
      // doesn't point to a valid instance for some reason (it shouldn't happend)
      itsOriginalVisAccessor->sync();
  }
}

/// @brief helper templated method to write back a cube to main table column
/// @details For now, it is only used in writeOriginalVis/Flag methods
/// and therefore can be kept in cc rather than tcc file (it is private, so
/// can be used by this class only). This can easily be changed in the future,
/// if need arises. This method encapsulates handling of channel selection
/// @param[in] cube Cube to work with, type should match the column type. Should be
///                 of the appropriate shape
/// @param[in] colName Name of the column
template<typename T>
void TableDataIterator::writeCube(const casacore::Cube<T> &cube,
                                  const std::string &colName) const
{
  const casacore::uInt nChan = nChannel();
  const casacore::uInt startChan = startChannel();
  // Setup a slicer to extract the specified channel range only
  const casacore::Slicer chanSlicer(casacore::Slice(),casacore::Slice(startChan,nChan));

  // no change of shape is permitted
  ASKAPASSERT(cube.nrow() == nRow() &&
              cube.ncolumn() == nChan &&
              cube.nplane() == nPol());
  casacore::ArrayColumn<T> visCol(getCurrentIteration(), colName);
  ASKAPDEBUGASSERT(getCurrentIteration().nrow() >= getCurrentTopRow()+
                    nRow());
  casacore::rownr_t tableRow = getCurrentTopRow();
  casacore::Matrix<T> buf(nPol(),nChan);
  for (casacore::uInt row=0;row<cube.nrow();++row,++tableRow) {
       const casacore::IPosition shape = visCol.shape(row);
       ASKAPDEBUGASSERT(shape.size() && (shape.size()<3));
       const casacore::uInt thisRowNumberOfPols = shape[0];
       const casacore::uInt thisRowNumberOfChannels = shape.size()>1 ? shape[1] : 1;
       if (thisRowNumberOfPols != cube.nplane()) {
           ASKAPTHROW(DataAccessError, "Current implementation of the writing to original "
                "visibilities does not support partial selection of the data");
       }
       if (thisRowNumberOfChannels < nChan + startChan) {
           ASKAPTHROW(DataAccessError, "Channel selection doesn't fit into exisiting visibility array");
       }
       // for now just copy
       bool useSlicer = (startChan!=0) || (startChan+nChan!=thisRowNumberOfChannels);

       for (casacore::uInt chan=0; chan<nChan; ++chan) {
            for (casacore::uInt pol=0; pol<thisRowNumberOfPols; ++pol) {
                 buf(pol,chan) = cube(row,chan,pol);
            }
       }
       if (useSlicer) {
           visCol.putSlice(tableRow,chanSlicer,buf);
       } else {
           visCol.put(tableRow, buf);
       }
  }
}

/// @brief write back the original visibilities
/// @details The write operation is possible if the shape of the
/// visibility cube stays the same as the shape of the data in the
/// table. The method uses DataAccessor to obtain a reference to the
/// visibility cube (hence no parameters).
void TableDataIterator::writeOriginalVis() const
{
   writeCube(getAccessor().visibility(), getDataColumnName());
}

/// @brief write back flags
/// @details The write operation is possible if the shape of the
/// flag cube stays the same as the shape of the data in the
/// table. The method uses DataAccessor to obtain a reference to the
/// visibility cube (hence no parameters).
/// @note This operation is specific to table (i.e MS) based implementaton
/// of the interface
void TableDataIterator::writeOriginalFlag() const
{
   const bool rowBasedFlagUsed = getCurrentIteration().tableDesc().isColumn("FLAG_ROW");
   const casacore::Cube<casacore::Bool>& flags = getAccessor().flag();
   if (rowBasedFlagUsed) {
       // check that updated flag doesn't contradict row-based flag
       casacore::ROScalarColumn<casacore::Bool> rowFlagCol(getCurrentIteration(), "FLAG_ROW");
       const casacore::Vector<casacore::Bool> rowBasedFlag = rowFlagCol.getColumn();
       const casacore::rownr_t topRow = getCurrentTopRow();
       ASKAPDEBUGASSERT(static_cast<casacore::rownr_t>(rowBasedFlag.nelements()) >= topRow+flags.nrow());
       for (casacore::uInt row = 0; row < flags.nrow(); ++row) {
            if (rowBasedFlag[row + topRow]) {
                bool oneUnflagged = false;
                casacore::Matrix<casacore::Bool> thisRow = flags.yzPlane(row);
                for (casacore::Matrix<casacore::Bool>::const_iterator ci = thisRow.begin();
                     ci != thisRow.end(); ++ci) {
                     if (!(*ci)) {
                         oneUnflagged = true;
                         break;
                     }
                }
                //std::cout<<row<<" "<<rowBasedFlag[row + topRow]<<" "<<oneUnflagged<<std::endl;
                ASKAPCHECK(!oneUnflagged, "Flag modification attempted to unflag data for the row ("<<
                      row<<") which is flagged via row-based flagging mechanism. This is not supported");
            }
       }

   }
   writeCube(flags, "FLAG");
}


/// @brief check whether one can write to the main table
/// @details Buffers held in subtables are not covered by this method.
/// @return true if write operation is allowed
bool TableDataIterator::mainTableWritable() const throw()
{
  try {
    return getCurrentIteration().isWritable();
  }
  catch (...) {}
  return false;
}
