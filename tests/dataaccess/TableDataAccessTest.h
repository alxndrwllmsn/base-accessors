/// @file TableDataAccessTest.h
/// $brief Tests of the table-based Data Accessor classes
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
#ifndef TABLE_DATA_ACCESS_TEST_H
#define TABLE_DATA_ACCESS_TEST_H

// boost includes
#include <boost/shared_ptr.hpp>

// casa includes
#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/Tables/TableError.h>
#include <casacore/casa/OS/EnvVar.h>

// std includes
#include <string>
#include <vector>

// cppunit includes
#include <cppunit/extensions/HelperMacros.h>

// own includes
#include <askap/dataaccess/DataAccessError.h>
#include <askap/dataaccess/TableInfoAccessor.h>
#include <askap/dataaccess/TableDataSource.h>
#include <askap/dataaccess/IConstDataSource.h>
#include <askap/dataaccess/TableConstDataIterator.h>
#include "TableTestRunner.h"

namespace askap {

namespace accessors {

class TableDataAccessTest : public CppUnit::TestFixture {

  CPPUNIT_TEST_SUITE(TableDataAccessTest);
  CPPUNIT_TEST(corrTypeSelectionTest);
  CPPUNIT_TEST(userDefinedIndexSelectionTest);
  CPPUNIT_TEST(uvDistanceSelectionTest);
  CPPUNIT_TEST(nonZeroMinUVSelectionTest);
  CPPUNIT_TEST(antennaSelectionTest);
  CPPUNIT_TEST_EXCEPTION(bufferManagerExceptionTest,casacore::TableError);
  CPPUNIT_TEST(bufferManagerTest);
  CPPUNIT_TEST(dataDescTest);
  CPPUNIT_TEST(spWindowTest);
  CPPUNIT_TEST(polarisationTest);
  CPPUNIT_TEST(feedTest);
  CPPUNIT_TEST(fieldTest);
  CPPUNIT_TEST(antennaTest);
  CPPUNIT_TEST(antennaPositionShortcutTest);
  CPPUNIT_TEST(originalVisRewriteTest);
  CPPUNIT_TEST(originalFlagRewriteTest);
  CPPUNIT_TEST(readOnlyTest);
  CPPUNIT_TEST(channelSelectionTest);
  CPPUNIT_TEST(freqSelectionTest);
  CPPUNIT_TEST(chunkSizeTest);
  CPPUNIT_TEST_SUITE_END();
public:

  /// set up the test suite
  void setUp();
  /// destruct the test suite
  void tearDown();
  /// test of correlation type selection
  void corrTypeSelectionTest();
  /// test of user-defined index selection
  void userDefinedIndexSelectionTest();
  /// test of selection based on uv-distance
  void uvDistanceSelectionTest();
  /// test of selection based on non-zero min uv-distance
  void nonZeroMinUVSelectionTest();
  /// test of selection based on antenna index
  void antennaSelectionTest();
  /// test of read only operations of the whole table-based implementation
  void readOnlyTest();
  /// test exception if disk-based buffers are requested for a read-only table
  void bufferManagerExceptionTest();
  /// extensive test of buffer operations
  void bufferManagerTest();
  /// test access to data description subtable
  void dataDescTest();
  /// test access to spectral window subtable
  void spWindowTest();
  /// test access to polarisation subtable
  void polarisationTest();
  /// test access to the feed subtable
  void feedTest();
  /// test access to the field subtable
  void fieldTest();
  /// test access to the antenna subtable
  void antennaTest();
  /// test access to antenna positions via a shortcut method
  void antennaPositionShortcutTest();
  /// test to rewrite original visibilities
  void originalVisRewriteTest();
  /// test to rewrite original flags
  void originalFlagRewriteTest();
  /// test read/write with channel selection
  void channelSelectionTest();
  /// test read/write with frequency selection
  void freqSelectionTest();
  /// test restriction of the chunk size
  void chunkSizeTest();
protected:
  void doBufferTest() const;
private:
  boost::shared_ptr<ITableInfoAccessor> itsTableInfoAccessor;
}; // class TableDataAccessTest

void TableDataAccessTest::setUp()
{
}

void TableDataAccessTest::tearDown()
{
  itsTableInfoAccessor.reset();
}

/// test of read only operations of the whole table-based implementation
void TableDataAccessTest::readOnlyTest()
{
  TableConstDataSource ds(TableTestRunner::msName());

  IDataConverterPtr conv=ds.createConverter();
  conv->setFrequencyFrame(casacore::MFrequency::Ref(casacore::MFrequency::BARY),"MHz");
  conv->setEpochFrame(casacore::MEpoch(casacore::Quantity(50257.29,"d"),
                      casacore::MEpoch::Ref(casacore::MEpoch::UTC)),"s");
  conv->setDirectionFrame(casacore::MDirection::Ref(casacore::MDirection::AZEL));

  int maxiter=4; // we don't need to read the whole dataset as it
                 // may take a long time. A few iterations should be sufficient.
                 // It is however useful to check that iteration finishes
                 // properly at the end of the measurement set. Hence, we
                 // continue iterating through the dataset without reading.
  casacore::MDirection testDir(casacore::MVDirection(0.12345,-0.12345), casacore::MDirection::J2000);
  casacore::MDirection testDir2(casacore::MVDirection(-0.12345,0.12345), casacore::MDirection::J2000);

  for (IConstDataSharedIter it=ds.createConstIterator(conv);it!=it.end();++it) {
       if (maxiter<0) {
           continue;
       }
       --maxiter;
       // just call several accessor methods to ensure that no exception is
       // thrown
       CPPUNIT_ASSERT(it->visibility().nrow() == it->nRow());
       CPPUNIT_ASSERT(it->visibility().ncolumn() == it->nChannel());
       CPPUNIT_ASSERT(it->visibility().nplane() == it->nPol());
       CPPUNIT_ASSERT(it->frequency().nelements() == it->nChannel());
       CPPUNIT_ASSERT(it->flag().shape() == it->visibility().shape());
       CPPUNIT_ASSERT(it->pointingDir2().nelements() == it->nRow());
       CPPUNIT_ASSERT(it->antenna1().nelements() == it->nRow());
       it->time();
       CPPUNIT_ASSERT(it->feed1PA().nelements() == it->nRow());
       CPPUNIT_ASSERT(it->noise().shape() == it->visibility().shape());
       CPPUNIT_ASSERT(it->rotatedUVW(testDir).nelements() == it->nRow());
       CPPUNIT_ASSERT(it->uvwRotationDelay(testDir,testDir2).nelements() == it->nRow());
       CPPUNIT_ASSERT(it->stokes().nelements() == it->nPol());
       CPPUNIT_ASSERT(it->nPol() == 2);
       CPPUNIT_ASSERT(it->stokes()[0] == casacore::Stokes::XX);
       CPPUNIT_ASSERT(it->stokes()[1] == casacore::Stokes::YY);

       // checks specific to table-based implementation
       boost::shared_ptr<TableConstDataIterator> actualIt = it.dynamicCast<TableConstDataIterator>();
       CPPUNIT_ASSERT(actualIt);
       CPPUNIT_ASSERT_EQUAL(0u, actualIt->currentFieldID());
       CPPUNIT_ASSERT_EQUAL(0u, actualIt->currentScanID());
  }
}

void TableDataAccessTest::userDefinedIndexSelectionTest()
{
  TableConstDataSource ds(TableTestRunner::msName());
  IDataSelectorPtr sel = ds.createSelector();
  sel->chooseUserDefinedIndex("ANTENNA1",1);
  for (IConstDataSharedIter it=ds.createConstIterator(sel);it!=it.end();++it) {
       for (casacore::uInt row=0;row<it->nRow();++row) {
            CPPUNIT_ASSERT(it->antenna1()[row] == 1);
       }
  }
  sel = ds.createSelector();
  sel->chooseCrossCorrelations();
  sel->chooseUserDefinedIndex("ANTENNA1",1);
  for (IConstDataSharedIter it=ds.createConstIterator(sel);it!=it.end();++it) {
       for (casacore::uInt row=0;row<it->nRow();++row) {
            CPPUNIT_ASSERT(it->antenna1()[row] == 1);
            CPPUNIT_ASSERT(it->antenna2()[row] != 1);
       }
  }
}

/// test restriction of the chunk size
void TableDataAccessTest::chunkSizeTest()
{
   TableConstDataSource ds(TableTestRunner::msName());
   IDataSelectorPtr sel = ds.createSelector();
   sel->chooseCrossCorrelations();
   const casacore::uInt nAnt = 6; // we have 6 antennas in the test dataset
   const casacore::uInt nBeams = 1; // we have 1 beam in the test dataset
   const casacore::uInt nRowsExpected = nBeams * nAnt * (nAnt - 1) / 2;
   casacore::uInt nIterOrig = 0;
   for (IConstDataSharedIter it=ds.createConstIterator(sel);it!=it.end();++it,++nIterOrig) {
        CPPUNIT_ASSERT_EQUAL(nRowsExpected, it->nRow());
        CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(nRowsExpected), it->visibility().nrow());
   }

   // restrict the chunk size for the following iterators
   ds.configureMaxChunkSize(nRowsExpected / 2);

   casacore::uInt count = 0;
   for (IConstDataSharedIter it=ds.createConstIterator(sel);it!=it.end();++it,++count) {
        const casacore::uInt nRowsExpectedThisIteration = ((nRowsExpected % 2 == 0) || (count % 3 != 2)) ? nRowsExpected / 2 : 1;
        if (count / 3 < nIterOrig) {
            // exclude the last iteration from the check as binning may be different
            CPPUNIT_ASSERT_EQUAL(nRowsExpectedThisIteration, it->nRow());
            CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(nRowsExpectedThisIteration), it->visibility().nrow());
        }
   }
}


/// test of correlation type selection
void TableDataAccessTest::corrTypeSelectionTest()
{
  TableConstDataSource ds(TableTestRunner::msName());
  IDataSelectorPtr sel = ds.createSelector();
  sel->chooseAutoCorrelations();
  for (IConstDataSharedIter it=ds.createConstIterator(sel);it!=it.end();++it) {
       for (casacore::uInt row=0;row<it->nRow();++row) {
            CPPUNIT_ASSERT(it->antenna1()[row] == it->antenna2()[row]);
            CPPUNIT_ASSERT(it->feed1()[row] == it->feed2()[row]);
       }
  }
  sel = ds.createSelector();
  sel->chooseCrossCorrelations();
  for (IConstDataSharedIter it=ds.createConstIterator(sel);it!=it.end();++it) {
       for (casacore::uInt row=0;row<it->nRow();++row) {
            CPPUNIT_ASSERT((it->antenna1()[row] != it->antenna2()[row]) ||
                           (it->feed1()[row] != it->feed2()[row]));
       }
  }
}
/// test of selection based on non-zero min uv-distance
void TableDataAccessTest::nonZeroMinUVSelectionTest()
{
  TableConstDataSource ds(TableTestRunner::msName());
  IDataSelectorPtr sel = ds.createSelector();
  sel->chooseMinNonZeroUVDistance(1000.);
  for (IConstDataSharedIter it=ds.createConstIterator(sel);it!=it.end();++it) {
       for (casacore::uInt row=0;row<it->nRow();++row) {
            const casacore::RigidVector<casacore::Double, 3> &uvw = it->uvw()(row);
            const casacore::Double uvDist = sqrt(casacore::square(uvw(0))+
                                             casacore::square(uvw(1)));
            // it's ok to compare doubles with zero here because the specific case with zero uvw
            // arises from direct assignment of 0. to each uvw coordinate and selection rule
            // also checks for match
            CPPUNIT_ASSERT(uvDist>=1000. || (uvw(0) == 0. && uvw(1) == 0. && uvw(2) == 0.));
       }
  }
  // explicit selection of auto-correlations to ensure zero uv gets through
  sel->chooseAutoCorrelations();
  casacore::uInt counter = 0;
  for (IConstDataSharedIter it=ds.createConstIterator(sel);it!=it.end();++it) {
       counter += it->nRow();
       for (casacore::uInt row=0;row<it->nRow();++row) {
            const casacore::RigidVector<casacore::Double, 3> &uvw = it->uvw()(row);
            const casacore::Double uvDist = sqrt(casacore::square(uvw(0))+
                                             casacore::square(uvw(1)));
            CPPUNIT_ASSERT(uvDist < 1e-6);
       }
  }
  CPPUNIT_ASSERT(counter > 0);
}

/// test of selection based on the minimum/maximum uv distance
void TableDataAccessTest::uvDistanceSelectionTest()
{
  TableConstDataSource ds(TableTestRunner::msName());
  IDataSelectorPtr sel = ds.createSelector();
  sel->chooseMinUVDistance(1000.);
  for (IConstDataSharedIter it=ds.createConstIterator(sel);it!=it.end();++it) {
       for (casacore::uInt row=0;row<it->nRow();++row) {
            const casacore::RigidVector<casacore::Double, 3> &uvw = it->uvw()(row);
            const casacore::Double uvDist = sqrt(casacore::square(uvw(0))+
                                             casacore::square(uvw(1)));
            CPPUNIT_ASSERT(uvDist>=1000.);
       }
  }
  sel = ds.createSelector();
  sel->chooseCrossCorrelations();
  sel->chooseMaxUVDistance(3000.);
  for (IConstDataSharedIter it=ds.createConstIterator(sel);it!=it.end();++it) {
       for (casacore::uInt row=0;row<it->nRow();++row) {
            const casacore::RigidVector<casacore::Double, 3> &uvw = it->uvw()(row);
            const casacore::Double uvDist = sqrt(casacore::square(uvw(0))+
                                             casacore::square(uvw(1)));
            CPPUNIT_ASSERT(uvDist<=3000.);
       }
  }
}

/// test of selection based on antenna index
void TableDataAccessTest::antennaSelectionTest()
{
  TableConstDataSource ds(TableTestRunner::msName());
  IDataSelectorPtr sel = ds.createSelector();
  sel->chooseAntenna(2u);
  for (IConstDataSharedIter it=ds.createConstIterator(sel);it!=it.end();++it) {
       for (casacore::uInt row=0;row<it->nRow();++row) {
            CPPUNIT_ASSERT((it->antenna1()[row] == 2u) || (it->antenna2()[row] == 2u));
       }
  }
  // and now checking that chaining different selectors work as expected
  sel = ds.createSelector();
  sel->chooseCrossCorrelations();
  sel->chooseAntenna(2u);
  for (IConstDataSharedIter it=ds.createConstIterator(sel);it!=it.end();++it) {
       casacore::uInt cntFirst = 0u;
       casacore::uInt cntSecond = 0u;
       for (casacore::uInt row=0;row<it->nRow();++row) {
            CPPUNIT_ASSERT((it->antenna1()[row] != it->antenna2()[row]) ||
                           (it->feed1()[row] != it->feed2()[row]));
            CPPUNIT_ASSERT((it->antenna1()[row] == 2u) || (it->antenna2()[row] == 2u));
            if (it->antenna1()[row] == 2u) {
                ++cntFirst;
            } else {
                ++cntSecond;
            }
       }
       // test dataset has 6 antennas, so 5 cross-correlations with antenna 2
       // first index is slow varying, so 3 cross-correlations will have antenna 2 as
       // the first index and 2 as the second
       CPPUNIT_ASSERT_EQUAL(5u, cntFirst + cntSecond);
       CPPUNIT_ASSERT_EQUAL(3u, cntFirst);
       CPPUNIT_ASSERT_EQUAL(2u, cntSecond);
  }
}


void TableDataAccessTest::bufferManagerExceptionTest()
{
  // test with the disk buffers, and leave the table read only. This
  // should throw a TableError
  itsTableInfoAccessor.reset(new TableInfoAccessor(
          casacore::Table(TableTestRunner::msName()),false));
  doBufferTest();
}

void TableDataAccessTest::bufferManagerTest()
{
  // first test with memory buffers
  itsTableInfoAccessor.reset(new TableInfoAccessor(
              casacore::Table(TableTestRunner::msName()),true));
  doBufferTest();
  // now test with the disk buffers
  itsTableInfoAccessor.reset(new TableInfoAccessor(
            casacore::Table(TableTestRunner::msName(),casacore::Table::Update), false));
  doBufferTest();
}

/// test access to data description subtable
void TableDataAccessTest::dataDescTest()
{
  // because we're not accessing the buffers here, it shouldn't really
  // matter whether we open it with memory buffers or with disk buffers
  // and read-only table should be enough.
  itsTableInfoAccessor.reset(new TableInfoAccessor(
              casacore::Table(TableTestRunner::msName()),false));
  ASKAPASSERT(itsTableInfoAccessor);
  const ITableDataDescHolder &dataDescription=itsTableInfoAccessor->
                    subtableInfo().getDataDescription();
  CPPUNIT_ASSERT(dataDescription.getSpectralWindowID(0)==0);
  CPPUNIT_ASSERT(dataDescription.getPolarizationID(0)==0);
  CPPUNIT_ASSERT(dataDescription.getDescIDsForSpWinID(0).size()==1);
  CPPUNIT_ASSERT(dataDescription.getDescIDsForSpWinID(1).size()==0);
}

/// test access to spectral window subtable
void TableDataAccessTest::spWindowTest()
{
  // because we're not accessing the buffers here, it shouldn't really
  // matter whether we open it with memory buffers or with disk buffers
  // and read-only table should be enough.
  itsTableInfoAccessor.reset(new TableInfoAccessor(
              casacore::Table(TableTestRunner::msName()),false));
  ASKAPASSERT(itsTableInfoAccessor);
  const ITableSpWindowHolder &spWindow=itsTableInfoAccessor->
                    subtableInfo().getSpWindow();
  CPPUNIT_ASSERT(spWindow.getReferenceFrame(0).getType() ==
                 casacore::MFrequency::TOPO);
  CPPUNIT_ASSERT(spWindow.getFrequencyUnit().getName() == "Hz");
  CPPUNIT_ASSERT(spWindow.getFrequencies(0).size() == 13);
  for (casacore::uInt chan=0;chan<13;++chan) {
     CPPUNIT_ASSERT(spWindow.getFrequencies(0)[chan] ==
              spWindow.getFrequency(0,chan).getValue().getValue());
  }
  CPPUNIT_ASSERT(fabs(spWindow.getFrequencies(0)[0]-1.4e9)<1e-5);
}

/// test access to polarisation subtable
void TableDataAccessTest::polarisationTest()
{
  // because we're not accessing the buffers here, it shouldn't really
  // matter whether we open it with memory buffers or with disk buffers
  // and read-only table should be enough.
  itsTableInfoAccessor.reset(new TableInfoAccessor(
              casacore::Table(TableTestRunner::msName()),false));
  ASKAPASSERT(itsTableInfoAccessor);
  const ITablePolarisationHolder &polHandler = itsTableInfoAccessor->
                      subtableInfo().getPolarisation();
  CPPUNIT_ASSERT(polHandler.nPol(0) == 2);
  casacore::Vector<casacore::Stokes::StokesTypes> polTypes = polHandler.getTypes(0);
  CPPUNIT_ASSERT(polHandler.nPol(0) == polTypes.nelements());
  for (casacore::uInt pol=0; pol<polHandler.nPol(0); ++pol) {
       CPPUNIT_ASSERT(polHandler.getType(0,pol) == polTypes[pol]);
  }
  CPPUNIT_ASSERT(polTypes[0] == casacore::Stokes::XX);
  CPPUNIT_ASSERT(polTypes[1] == casacore::Stokes::YY);
}

/// test access to the feed subtable
void TableDataAccessTest::feedTest()
{
  // because we're not accessing the buffers here, it shouldn't really
  // matter whether we open it with memory buffers or with disk buffers
  // and read-only table should be enough.
  itsTableInfoAccessor.reset(new TableInfoAccessor(
              casacore::Table(TableTestRunner::msName()),false));
  CPPUNIT_ASSERT(itsTableInfoAccessor);
  const IFeedSubtableHandler &feedSubtable=itsTableInfoAccessor->
                      subtableInfo().getFeed();
  casacore::MEpoch time(casacore::MVEpoch(casacore::Quantity(50257.29,"d")),
                    casacore::MEpoch::Ref(casacore::MEpoch::UTC));
  for (casacore::uInt feed=0; feed<5; ++feed) {
       for (casacore::uInt ant=1; ant<6; ++ant) {
            CPPUNIT_ASSERT(fabs(feedSubtable.getBeamOffset(time,0,ant,feed)(0)-
                       feedSubtable.getBeamOffset(time,0,0,feed)(0))<1e-7);
            CPPUNIT_ASSERT(fabs(feedSubtable.getBeamOffset(time,0,ant,feed)(1)-
                       feedSubtable.getBeamOffset(time,0,0,feed)(1))<1e-7);
            CPPUNIT_ASSERT(fabs(feedSubtable.getBeamPA(time,0,ant,feed)-
                       feedSubtable.getBeamPA(time,0,0,feed))<1e-7);
       }
       if (feed!=4) {
           CPPUNIT_ASSERT(fabs(feedSubtable.getBeamOffset(time,0,0,
                                             feed)(0))*206265-900<1e-5);
           CPPUNIT_ASSERT(fabs(feedSubtable.getBeamOffset(time,0,0,
                                             feed)(1))*206265-900<1e-5);
       } else {
           CPPUNIT_ASSERT(fabs(feedSubtable.getBeamOffset(time,0,0,
                                             feed)(0))<1e-5);
           CPPUNIT_ASSERT(fabs(feedSubtable.getBeamOffset(time,0,0,
                                             feed)(1))<1e-5);
       }
       CPPUNIT_ASSERT(fabs(feedSubtable.getBeamPA(time,0,0,feed))<1e-5);
  }
}

/// test access to the field subtable
void TableDataAccessTest::fieldTest()
{
  // because we're not accessing the buffers here, it shouldn't really
  // matter whether we open it with memory buffers or with disk buffers
  // and read-only table should be enough.
  itsTableInfoAccessor.reset(new TableInfoAccessor(
              casacore::Table(TableTestRunner::msName()),false));
  CPPUNIT_ASSERT(itsTableInfoAccessor);
  const IFieldSubtableHandler &fieldSubtable=itsTableInfoAccessor->
                      subtableInfo().getField();
  casacore::MEpoch time(casacore::MVEpoch(casacore::Quantity(50257.29,"d")),
                    casacore::MEpoch::Ref(casacore::MEpoch::UTC));
  casacore::MVDirection refDir(casacore::Quantity(0.,"deg"), casacore::Quantity(-50.,"deg"));
  CPPUNIT_ASSERT(fieldSubtable.getReferenceDir(time).getRef().getType() ==
                 casacore::MDirection::J2000);
  CPPUNIT_ASSERT(fieldSubtable.getReferenceDir(time).getValue().
                 separation(refDir)<1e-7);

  // test random access (for row=0)
  CPPUNIT_ASSERT(fieldSubtable.getReferenceDir(0).getRef().getType() ==
                 casacore::MDirection::J2000);
  CPPUNIT_ASSERT(fieldSubtable.getReferenceDir(0).getValue().
                 separation(refDir)<1e-7);
}


void TableDataAccessTest::doBufferTest() const
{
  const IBufferManager &bufferMgr=itsTableInfoAccessor->subtableInfo().
                                  getBufferManager();
  const casacore::uInt index=5;
  CPPUNIT_ASSERT(!bufferMgr.bufferExists("TEST",index));
  casacore::Cube<casacore::Complex> vis(5,10,2);
  vis.set(casacore::Complex(1.,-0.5));
  bufferMgr.writeBuffer(vis,"TEST",index);
  CPPUNIT_ASSERT(bufferMgr.bufferExists("TEST",index));
  casacore::Cube<casacore::Complex> vis2(5,1,2);
  vis2.set(casacore::Complex(-1.,0.5));
  CPPUNIT_ASSERT(!bufferMgr.bufferExists("TEST",index-1));
  bufferMgr.writeBuffer(vis2,"TEST",index-1);
  CPPUNIT_ASSERT(bufferMgr.bufferExists("TEST",index-1));
  bufferMgr.readBuffer(vis,"TEST",index-1);
  bufferMgr.readBuffer(vis2,"TEST",index);
  CPPUNIT_ASSERT(vis.shape()==casacore::IPosition(3,5,1,2));
  CPPUNIT_ASSERT(vis2.shape()==casacore::IPosition(3,5,10,2));
  ASKAPDEBUGASSERT(vis.shape()[0]>=0);
  ASKAPDEBUGASSERT(vis.shape()[1]>=0);
  ASKAPDEBUGASSERT(vis.shape()[2]>=0);
  for (casacore::uInt x=0;x<casacore::uInt(vis.shape()[0]);++x) {
       for (casacore::uInt y=0;y<casacore::uInt(vis.shape()[1]);++y) {
            for (casacore::uInt z=0;z<casacore::uInt(vis.shape()[2]);++z) {
	         CPPUNIT_ASSERT(abs(vis2(x,y,z)+vis(x,0,z))<1e-9);
	    }
       }
  }
}

/// test access to the antenna subtable
void TableDataAccessTest::antennaTest()
{
  // because we're not accessing the buffers here, it shouldn't really
  // matter whether we open it with memory buffers or with disk buffers
  // and read-only table should be enough.
  itsTableInfoAccessor.reset(new TableInfoAccessor(
              casacore::Table(TableTestRunner::msName()),false));
  CPPUNIT_ASSERT(itsTableInfoAccessor);
  const IAntennaSubtableHandler &antennaSubtable=itsTableInfoAccessor->
                      subtableInfo().getAntenna();
  for (casacore::uInt ant=0;ant<6;++ant) {
       CPPUNIT_ASSERT(antennaSubtable.getMount(ant) == "ALT-AZ");
      for (casacore::uInt ant2=0; ant2<ant; ++ant2) {
           CPPUNIT_ASSERT(antennaSubtable.getPosition(ant).getValue().
              separation(antennaSubtable.getPosition(ant2).getValue(),"deg").
                         getValue()<0.1);
      }
  }
}

/// test access to antenna positions via a short cut method, specific to
/// table-based implementation
void TableDataAccessTest::antennaPositionShortcutTest() {
  TableConstDataSource ds(TableTestRunner::msName());
  // this depends on the content of the test measurement set
  std::vector<casacore::MVPosition> expectation;
  expectation.push_back(casacore::MVPosition(-4.7522e+06, 2.79072e+06, -3.20048e+06));
  expectation.push_back(casacore::MVPosition(-4.75193e+06, 2.79118e+06, -3.20048e+06));
  expectation.push_back(casacore::MVPosition(-4.75155e+06, 2.79183e+06, -3.20048e+06));
  expectation.push_back(casacore::MVPosition(-4.75107e+06, 2.79264e+06, -3.20048e+06));
  expectation.push_back(casacore::MVPosition(-4.75092e+06, 2.79291e+06, -3.20048e+06));
  expectation.push_back(casacore::MVPosition(-4.7496e+06, 2.79514e+06, -3.20048e+06));
  // try non-short cut approach to test it too + test the number of antennas
  IConstDataSharedIter it=ds.createConstIterator();
  boost::shared_ptr<TableConstDataIterator> tabIt = it.dynamicCast<TableConstDataIterator>();
  CPPUNIT_ASSERT(tabIt);
  boost::shared_ptr<const ITableManager> mgr = tabIt->getTableManager();
  CPPUNIT_ASSERT(mgr);
  const casacore::uInt nAnt = mgr->getAntenna().getNumberOfAntennas();
  CPPUNIT_ASSERT_EQUAL(static_cast<casacore::uInt>(expectation.size()), nAnt);
  CPPUNIT_ASSERT_EQUAL(nAnt, ds.getNumberOfAntennas());
  //

  for (casacore::uInt ant=0;ant<expectation.size();++ant) {
       const casacore::MPosition pos = ds.getAntennaPosition(ant);
       CPPUNIT_ASSERT_EQUAL(casacore::MPosition::Ref(casacore::MPosition::ITRF).getType(), pos.getRef().getType());
       const casacore::MVPosition diff = pos.getValue() - expectation[ant];
       CPPUNIT_ASSERT_DOUBLES_EQUAL(0., diff.getLength().getValue("m") / expectation[ant].getLength().getValue("m"), 1e-5);
  }
}

/// test read/write with channel selection
void TableDataAccessTest::channelSelectionTest()
{
  TableDataSource tds(TableTestRunner::msName(), TableDataSource::WRITE_PERMITTED);
  IDataSource &ds=tds; // to have all interface methods available without
                       // ambiguity (otherwise methods overridden in
                       // TableDataSource would get a priority)
  for (IDataSharedIter it=ds.createIterator(); it!=it.end(); ++it) {
       // store original visibilities in a buffer
       it.buffer("BACKUP").rwVisibility() = it->visibility();
       // set new values for all spectral channels, rows and pols
       it->rwVisibility().set(casacore::Complex(1.,0.5));
  }


  IDataSelectorPtr sel = ds.createSelector();
  ASKAPASSERT(sel);
  sel->chooseChannels(2, 3);
  for (IDataSharedIter it=ds.createIterator(sel); it!=it.end(); ++it) {
       // different value corresponding to selected channels
       it->rwVisibility().set(casacore::Complex(-0.5,1.0));
  }

  // check that the visibilities are set to a required constant for the selected subset of channels
  for (IConstDataSharedIter cit = ds.createConstIterator(sel);
                                        cit != cit.end(); ++cit) {
       const casacore::Cube<casacore::Complex> &vis = cit->visibility();
       // selected just two channels
       ASKAPASSERT(vis.ncolumn() == 2);
       for (casacore::uInt row = 0; row < vis.nrow(); ++row) {
            for (casacore::uInt column = 0; column < vis.ncolumn(); ++column) {
                 for (casacore::uInt plane = 0; plane < vis.nplane(); ++plane) {
                      CPPUNIT_ASSERT(abs(vis(row,column,plane)-
                                         casacore::Complex(-0.5,1.0))<1e-7);
                 }
            }
       }
  }

  // check that the visibilities are set to a required constant in the whole cube
  for (IConstDataSharedIter cit = ds.createConstIterator();
                                        cit != cit.end(); ++cit) {
       const casacore::Cube<casacore::Complex> &vis = cit->visibility();
       // selected just two channels
       ASKAPASSERT(vis.ncolumn() == 13);
       for (casacore::uInt row = 0; row < vis.nrow(); ++row) {
            for (casacore::uInt column = 0; column < vis.ncolumn(); ++column) {
                 for (casacore::uInt plane = 0; plane < vis.nplane(); ++plane) {
                      const casacore::Complex result = (column==3) || (column==4) ?
                                  casacore::Complex(-0.5,1.0) : casacore::Complex(1.0,0.5);
                      CPPUNIT_ASSERT(abs(vis(row,column,plane) - result)<1e-7);
                 }
            }
       }
  }



  // set visibilities back to the original values
  for (IDataSharedIter it=ds.createIterator(); it!=it.end(); ++it) {
       // store original visibilities in a buffer
       it->rwVisibility() = it.buffer("BACKUP").visibility();
  }
}

/// test read/write with frequency selection
void TableDataAccessTest::freqSelectionTest()
{
  TableDataSource tds(TableTestRunner::msName(), TableDataSource::WRITE_PERMITTED);
  IDataSource &ds=tds; // to have all interface methods available without
                       // ambiguity (otherwise methods overridden in
                       // TableDataSource would get a priority)
  casacore::Vector<casacore::Double> freqs;
  for (IDataSharedIter it=ds.createIterator(); it!=it.end(); ++it) {
       // store original visibilities in a buffer
       it.buffer("BACKUP").rwVisibility() = it->visibility();
       // set new values for all spectral channels, rows and pols
       it->rwVisibility().set(casacore::Complex(1.,0.5));
       if (freqs.nelements()==0) {
           freqs = it->frequency();
       }
  }


  IDataSelectorPtr sel = ds.createSelector();
  ASKAPASSERT(sel);
  // choose 3rd freq (i.e. freq(2)), note only 0 width is supported at present
  sel->chooseFrequencies(1, casacore::MVFrequency(freqs(2)), 0.);
  for (IDataSharedIter it=ds.createIterator(sel); it!=it.end(); ++it) {
       // different value corresponding to selected channels
       it->rwVisibility().set(casacore::Complex(-0.5,1.0));
  }

  // check that the visibilities are set to a required constant for the selected subset of channels
  for (IConstDataSharedIter cit = ds.createConstIterator(sel);
                                        cit != cit.end(); ++cit) {
       const casacore::Cube<casacore::Complex> &vis = cit->visibility();
       // selected just two channels
       ASKAPASSERT(vis.ncolumn() == 1);
       for (casacore::uInt row = 0; row < vis.nrow(); ++row) {
            for (casacore::uInt column = 0; column < vis.ncolumn(); ++column) {
                 for (casacore::uInt plane = 0; plane < vis.nplane(); ++plane) {
                      CPPUNIT_ASSERT(abs(vis(row,column,plane)-
                                         casacore::Complex(-0.5,1.0))<1e-7);
                 }
            }
       }
  }

  // check that the visibilities are set to a required constant in the whole cube
  for (IConstDataSharedIter cit = ds.createConstIterator();
                                        cit != cit.end(); ++cit) {
       const casacore::Cube<casacore::Complex> &vis = cit->visibility();
       // selected just two channels
       ASKAPASSERT(vis.ncolumn() == 13);
       for (casacore::uInt row = 0; row < vis.nrow(); ++row) {
            for (casacore::uInt column = 0; column < vis.ncolumn(); ++column) {
                 for (casacore::uInt plane = 0; plane < vis.nplane(); ++plane) {
                      const casacore::Complex result = (column==2 ?
                                  casacore::Complex(-0.5,1.0) : casacore::Complex(1.0,0.5));
                      CPPUNIT_ASSERT(abs(vis(row,column,plane) - result)<1e-7);
                 }
            }
       }
  }



  // set visibilities back to the original values
  for (IDataSharedIter it=ds.createIterator(); it!=it.end(); ++it) {
       // store original visibilities in a buffer
       it->rwVisibility() = it.buffer("BACKUP").visibility();
  }
}

void TableDataAccessTest::originalFlagRewriteTest()
{
   TableDataSource tds(TableTestRunner::msName(), TableDataSource::WRITE_PERMITTED);
   IDataSource &ds=tds; // to have all interface methods available without
                        // ambiguity (otherwise methods overridden in
                        // TableDataSource would get a priority)
   casacore::uInt iterCntr=0;
   for (IDataSharedIter it=ds.createIterator(); it!=it.end(); ++it,++iterCntr) {
        // first check that read-only and read-write access return the same data
        const casacore::Cube<casacore::Bool>& roFlags = it->flag();
        IFlagDataAccessor& acc = dynamic_cast<IFlagDataAccessor&>(*it);
        const casacore::Cube<casacore::Bool>& rwFlags = acc.rwFlag();
        CPPUNIT_ASSERT(roFlags.shape() == rwFlags.shape());
        CPPUNIT_ASSERT(roFlags.shape() == casacore::IPosition(3,it->nRow(), it->nChannel(), it->nPol()));
        for (casacore::uInt row=0; row < it->nRow(); ++row) {
             for (casacore::uInt chan=0; chan < it->nChannel(); ++chan) {
                  for (casacore::uInt pol =0; pol < it->nPol(); ++pol) {
                       CPPUNIT_ASSERT_EQUAL(roFlags(row,chan,pol), rwFlags(row,chan,pol));
                  }
             }
        }
   }
   CPPUNIT_ASSERT_EQUAL(420u, iterCntr);

   casacore::Vector<casacore::Cube<casacore::Bool> > memoryBuffer(iterCntr);
   iterCntr=0;
   for (IDataSharedIter it=ds.createIterator(); it!=it.end(); ++it,++iterCntr) {
        // copy the original flags
        const casacore::Cube<casacore::Bool>& roFlags = it->flag();
        memoryBuffer[iterCntr] = roFlags;
        IFlagDataAccessor& acc = dynamic_cast<IFlagDataAccessor&>(*it);
        casacore::Cube<casacore::Bool>& rwFlags = acc.rwFlag();
        for (casacore::uInt row=0; row < it->nRow(); ++row) {
             for (casacore::uInt chan=0; chan < it->nChannel(); ++chan) {
                  for (casacore::uInt pol =0; pol < it->nPol(); ++pol) {
                       // the test dataset uses row-based flagging mechanism, so can't just flip
                       // the flag to opposite - just flag all samples for the test
                       rwFlags(row,chan,pol) = true;
                  }
             }
        }
   }

   iterCntr=0;
   for (IDataSharedIter it=ds.createIterator(); it!=it.end(); ++it,++iterCntr) {
        const casacore::Cube<casacore::Bool>& roFlags = it->flag();
        IFlagDataAccessor& acc = dynamic_cast<IFlagDataAccessor&>(*it);
        casacore::Cube<casacore::Bool>& rwFlags = acc.rwFlag();
        for (casacore::uInt row=0; row < it->nRow(); ++row) {
             for (casacore::uInt chan=0; chan < it->nChannel(); ++chan) {
                  for (casacore::uInt pol =0; pol < it->nPol(); ++pol) {
                       // check that the flag is now always set,  then reset if it was unflagged originally
                       CPPUNIT_ASSERT_EQUAL(true, roFlags(row,chan,pol));
                       if (!memoryBuffer[iterCntr](row,chan,pol)) {
                           rwFlags(row,chan,pol) = false;
                       }
                  }
             }
        }
   }

}

/// test to rewrite original visibilities
void TableDataAccessTest::originalVisRewriteTest()
{
  TableDataSource tds(TableTestRunner::msName(), TableDataSource::WRITE_PERMITTED);
  IDataSource &ds=tds; // to have all interface methods available without
                       // ambiguity (otherwise methods overridden in
                       // TableDataSource would get a priority)
  casacore::uInt iterCntr=0;
  for (IDataSharedIter it=ds.createIterator(); it!=it.end(); ++it,++iterCntr) {
       // store original visibilities in a buffer
       it.buffer("BACKUP").rwVisibility() = it->visibility();
  }
  casacore::Vector<casacore::Cube<casacore::Complex> > memoryBuffer(iterCntr);
  iterCntr=0;
  for (IDataSharedIter it = ds.createIterator(); it!=it.end(); ++it,++iterCntr) {
       // save original values in memory to check buffers as well
       memoryBuffer[iterCntr] = it->visibility();
       // reset visibilities to a constant
       it->rwVisibility().set(casacore::Complex(1.,0.5));
  }
  // check that the visibilities are set to a required constant
  for (IConstDataSharedIter cit = ds.createConstIterator();
                                        cit != cit.end(); ++cit) {
       const casacore::Cube<casacore::Complex> &vis = cit->visibility();
       for (casacore::uInt row = 0; row < vis.nrow(); ++row) {
            for (casacore::uInt column = 0; column < vis.ncolumn(); ++column) {
                 for (casacore::uInt plane = 0; plane < vis.nplane(); ++plane) {
                      CPPUNIT_ASSERT(abs(vis(row,column,plane)-
                                         casacore::Complex(1.,0.5))<1e-7);
                 }
            }
       }
  }
  // set visibilities back to the original values
  for (IDataSharedIter it=ds.createIterator(); it!=it.end(); ++it) {
       // store original visibilities in a buffer
       it->rwVisibility() = it.buffer("BACKUP").visibility();
  }

  // compare with the values stored in the memory
  iterCntr=0;
  for (IConstDataSharedIter cit = ds.createConstIterator();
                                  cit != cit.end(); ++cit,++iterCntr) {
       const casacore::Cube<casacore::Complex> &vis = cit->visibility();
       for (casacore::uInt row = 0; row < vis.nrow(); ++row) {
            for (casacore::uInt column = 0; column < vis.ncolumn(); ++column) {
                 for (casacore::uInt plane = 0; plane < vis.nplane(); ++plane) {
                      CPPUNIT_ASSERT(abs(vis(row,column,plane)-
                             memoryBuffer[iterCntr](row,column,plane))<1e-7);
                 }
            }
       }
  }
}

} // namespace accessors

} // namespace askap

#endif // #ifndef TABLE_DATA_ACCESS_TEST_H
