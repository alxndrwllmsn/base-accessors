/// @file 
/// $brief Tests of the multi-chunk iterator adapter
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

#ifndef STACKABLE_ADAPTER_TEST_H
#define STACKABLE_ADAPTER_TEST_H

// boost includes
#include <boost/shared_ptr.hpp>

// cppunit includes
#include <cppunit/extensions/HelperMacros.h>
// own includes
#include <askap/dataaccess/TableDataSource.h>
#include <askap/dataaccess/TableConstDataAccessor.h>
#include <askap/dataaccess/IConstDataSource.h>
#include <askap/dataaccess/TimeChunkIteratorAdapter.h>
#include <askap/dataaccess/MemBufferDataAccessorStackable.h>
#include <askap/askap/AskapError.h>
#include "TableTestRunner.h"
#include <askap/askap/AskapUtil.h>

#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/Tables/TableIter.h>
#include <casacore/casa/Containers/Block.h>
#include <casacore/casa/BasicSL/String.h>


namespace askap {

namespace accessors {

class StackableAdapterTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(StackableAdapterTest);
  CPPUNIT_TEST(testInput);
  CPPUNIT_TEST(testInstantiate);
  CPPUNIT_TEST(testConstInstantiate);
  CPPUNIT_TEST(testStack);
  CPPUNIT_TEST(testCompare);
  CPPUNIT_TEST_SUITE_END();
protected:
  static size_t countSteps(const IConstDataSharedIter &it) {
     size_t counter;
     for (counter = 0; it!=it.end(); ++it,++counter) {}
     return counter;     
  }
  static bool testAcc(MemBufferDataAccessorStackable &test) {
    test.setAccessorIndex(2);
    // testing channel 2 (0 based) and baseline (2) (0 based)
    // The UVW for this baseline and this input data set is:
    // [-218.044021106325, 975.585041111335 , 826.584555325564]
    // the data complex number should be
    // (0.351497501134872,0.0155263254418969)
    // UVW test
    CPPUNIT_ASSERT_DOUBLES_EQUAL(double(-218.044021106325),double(test.uvw()[2](0)),1e-9);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(double(975.585041111335),double(test.uvw()[2](1)),1e-9);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(double(826.584555325564),double(test.uvw()[2](2)),1e-9);
    // Data value test
    CPPUNIT_ASSERT_DOUBLES_EQUAL(double(0.351497501134872),double(test.rwVisibility().at(2,2,0).real()),1e-9);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(double(0.0155263254418969),double(test.rwVisibility().at(2,2,0).imag()),1e-9);
  }
  
  
public:
  
  void testInstantiate() {
    // Test the instatiation and the auto stacking in the constructor
    // The input data source
    TableDataSource ds(TableTestRunner::msName());
    // An iterator into the data source
    IDataSharedIter it = ds.createIterator();

    MemBufferDataAccessorStackable adapter(it);
    
    // This should have buffered all the input visibilities
   
    testAcc(adapter);
    
    
  }
  void testConstInstantiate() {
    // Test the instatiation and the auto stacking in the constructor
    // The input data source
    // Also tests some pointer semantics
    TableConstDataSource ds(TableTestRunner::msName());
    // An iterator into the data source
    IConstDataSharedIter it = ds.createConstIterator();

    boost::shared_ptr<MemBufferDataAccessorStackable> adapter(new MemBufferDataAccessorStackable(it));
    
    testAcc(*adapter);
    
    
  }
  void testStack() {
    
    // The input data source
    TableConstDataSource ds(TableTestRunner::msName());
    // An iterator into the data source
    IConstDataSharedIter it = ds.createConstIterator();
    size_t count=0;
    vector<MemBufferDataAccessor> theStack;
    boost::shared_ptr<MemBufferDataAccessorStackable> adapter(new MemBufferDataAccessorStackable(*it));
    for (;it != it.end(); it.next()) {
      MemBufferDataAccessor acc(*it);
      acc.rwVisibility() = it->visibility();
      count++;
      adapter->append(acc);
    }
      
    testAcc(*adapter);
  }
  void testCompare() {
    TableConstDataSource ds(TableTestRunner::msName());
    IConstDataSharedIter it = ds.createConstIterator();
    MemBufferDataAccessorStackable adapter(it);
    // compare the contents.
    int index = 0;
    for (it.init() ;it != it.end(); it.next(),index++) {
      adapter.setAccessorIndex(index);
      
      for (casacore::uInt row=0;row<it->nRow();++row) {
        
        const casacore::RigidVector<casacore::Double, 3> &uvw = it->uvw()(row);
        const casacore::Double uvDist = sqrt(casacore::square(uvw(0))+
                                       casacore::square(uvw(1)));
        
        const casacore::RigidVector<casacore::Double, 3> &adapteruvw = adapter.uvw()(row);
        const casacore::Double adapteruvDist = sqrt(casacore::square(adapteruvw(0))+
        casacore::square(adapteruvw(1)));
        CPPUNIT_ASSERT_DOUBLES_EQUAL(double(adapteruvDist),double(uvDist),1e-7);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(double(it->time()),double(adapter.time()),1e-1);
        
        // lets just test if the rotatedUVW are in the correct spot.
        const casacore::MDirection fakeTangent(it->dishPointing1()[0], casacore::MDirection::J2000);
        const casacore::Vector<casacore::RigidVector<casacore::Double, 3> >& Ruvw = it->rotatedUVW(fakeTangent);
        const casacore::Vector<casacore::RigidVector<casacore::Double, 3> >& adapterRuvw = adapter.rotatedUVW(fakeTangent);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(double(Ruvw(row)(0)), double(adapterRuvw(row)(0)), 1e-9);
        
        
        
      }
    }
  }
  void testInput() {
     TableConstDataSource ds(TableTestRunner::msName());
     IDataConverterPtr conv=ds.createConverter();
     conv->setEpochFrame(); // ensures seconds since 0 MJD
     CPPUNIT_ASSERT_EQUAL(size_t(420), countSteps(ds.createConstIterator(conv)));
     boost::shared_ptr<TimeChunkIteratorAdapter> it(new TimeChunkIteratorAdapter(ds.createConstIterator(conv)));
     CPPUNIT_ASSERT_EQUAL(size_t(420), countSteps(it));
     it.reset(new TimeChunkIteratorAdapter(ds.createConstIterator(conv),599));
     size_t counter = 0;
     for (;it->moreDataAvailable();++counter) {
          CPPUNIT_ASSERT_EQUAL(size_t(1),countSteps(it));
          if (it->moreDataAvailable()) {
              it->resume();
          }
     }
     CPPUNIT_ASSERT_EQUAL(size_t(420), counter);     
     // now trying bigger chunks
     it.reset(new TimeChunkIteratorAdapter(ds.createConstIterator(conv),5990));
     for (counter = 0; it->moreDataAvailable(); ++counter) {
          CPPUNIT_ASSERT_EQUAL(size_t(10),countSteps(it));
          if (it->moreDataAvailable()) {
              it->resume();
          }
     }
     CPPUNIT_ASSERT_EQUAL(size_t(42), counter);     
     
  }
/*
  void testReadOnlyBuffer() {
     TableConstDataSource ds(TableTestRunner::msName());
     boost::shared_ptr<TimeChunkIteratorAdapter> it(new TimeChunkIteratorAdapter(ds.createConstIterator()));
     // this should generate an exception
     it->buffer("TEST");
  }

  void testReadOnlyAccessor() {
     TableConstDataSource ds(TableTestRunner::msName());
     boost::shared_ptr<TimeChunkIteratorAdapter> it(new TimeChunkIteratorAdapter(ds.createConstIterator()));
     boost::shared_ptr<IDataAccessor> acc;
     try {
        boost::shared_ptr<IDataAccessor> staticAcc(it->operator->(),utility::NullDeleter());
        ASKAPASSERT(staticAcc);
        acc = staticAcc;
     }
     catch (const AskapError &) {
        // just to ensure no exception is thrown from the try-block
        CPPUNIT_ASSERT(false);
     }
     CPPUNIT_ASSERT(acc);
     // this should generate an exception
     acc->rwVisibility();
  }
  
  void testNoResume() {
     TableConstDataSource ds(TableTestRunner::msName());
     IDataConverterPtr conv=ds.createConverter();
     conv->setEpochFrame(); // ensures seconds since 0 MJD
     boost::shared_ptr<TimeChunkIteratorAdapter> it(new TimeChunkIteratorAdapter(ds.createConstIterator(conv),5990));     
     try {
        // this code shouldn't throw AskapError 
        CPPUNIT_ASSERT(it->hasMore());        
        boost::shared_ptr<IConstDataIterator> cit = boost::dynamic_pointer_cast<IConstDataIterator>(it);
        CPPUNIT_ASSERT(cit);
        cit->next();
        CPPUNIT_ASSERT(cit->hasMore());        
        // just to access some data field
        *(*cit);
        (*cit)->antenna1();
        // access those fields directly
        *(*it);
        (*it)->antenna1();
        CPPUNIT_ASSERT_EQUAL(size_t(9),countSteps(it));
        CPPUNIT_ASSERT(it->moreDataAvailable());
        CPPUNIT_ASSERT(!it->hasMore());
     }
     catch (const AskapError &) {
        CPPUNIT_ASSERT(false);
     }
     // the following should throw AskapError
     CPPUNIT_ASSERT(it->next());
  }
*/
};

} // namespace accessors

} // namespace askap

#endif // #ifndef STACKABLE_ADAPTER_TEST_H


