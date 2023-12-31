/// @file
///
/// Unit test for the table-based implementation of the interface to access
/// calibration solutions
///
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

#include <casacore/casa/aipstype.h>
#include <cppunit/extensions/HelperMacros.h>
#include <askap/calibaccess/TableCalSolutionSource.h>
#include <askap/calibaccess/ParsetCalSolutionAccessor.h>
#include <askap/calibaccess/JonesIndex.h>
#include <askap/calibaccess/ChanAdapterCalSolutionConstSource.h>
#include <casacore/tables/Tables/Table.h>



#include <boost/shared_ptr.hpp>
#include <string>

namespace askap {

namespace accessors {

class TableCalSolutionTest : public CppUnit::TestFixture
{
   CPPUNIT_TEST_SUITE(TableCalSolutionTest);
   CPPUNIT_TEST(testCreate);
   CPPUNIT_TEST(testRead);
   CPPUNIT_TEST(testBlankEntries);
   CPPUNIT_TEST(testTrailingBlankEntry);
   CPPUNIT_TEST(testChanAdapterRead);
   CPPUNIT_TEST(testDelayedWrite);
//   CPPUNIT_TEST_EXCEPTION(testChanAdapterUndefinedBandpass, AskapError);
   CPPUNIT_TEST_EXCEPTION(testUndefinedGains, AskapError);
   CPPUNIT_TEST_EXCEPTION(testUndefinedLeakages, AskapError);
//   CPPUNIT_TEST_EXCEPTION(testUndefinedBandpasses, AskapError);
   CPPUNIT_TEST_EXCEPTION(testUndefinedSolution, AskapError);
   CPPUNIT_TEST_EXCEPTION(testTooFarIntoThePast, AskapError);
   CPPUNIT_TEST(testCreateManyRows);
   CPPUNIT_TEST(testCreateManyRows1);
   CPPUNIT_TEST_SUITE_END();
protected:

   static boost::shared_ptr<ICalSolutionSource> rwSource(bool doRemove) {
       const std::string fname("calibdata.tab");
       if (doRemove) {
           TableCalSolutionSource::removeOldTable(fname);
       }
       boost::shared_ptr<TableCalSolutionSource> css(new TableCalSolutionSource(fname,6,3,8));
       CPPUNIT_ASSERT(css);
       return css;
   }

   static boost::shared_ptr<ICalSolutionConstSource> roSource() {
       const std::string fname("calibdata.tab");
       boost::shared_ptr<TableCalSolutionConstSource> css(new TableCalSolutionConstSource(fname));
       CPPUNIT_ASSERT(css);
       return css;
   }

   static void testComplex(const casacore::Complex &expected, const casacore::Complex &obtained, const float tol = 1e-5) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(real(expected),real(obtained),tol);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(imag(expected),imag(obtained),tol);
   }

   static boost::shared_ptr<ICalSolutionConstAccessor> accessorForExistingTable() {
       const boost::shared_ptr<ICalSolutionConstSource> css = roSource();
       CPPUNIT_ASSERT(css);
       const long sID = css->mostRecentSolution();
       CPPUNIT_ASSERT_EQUAL(3l, sID);
       const boost::shared_ptr<ICalSolutionConstAccessor> acc = css->roSolution(sID);
       CPPUNIT_ASSERT(acc);
       return acc;
   }

public:

   void testCreate() {
       boost::shared_ptr<ICalSolutionSource> css = rwSource(true);
       long newID = css->newSolutionID(0.);
       CPPUNIT_ASSERT_EQUAL(0l, newID);
       boost::shared_ptr<ICalSolutionAccessor> acc = css->rwSolution(newID);
       acc->setGain(JonesIndex(0u,0u),JonesJTerm(casacore::Complex(1.0,-1.0),true,casacore::Complex(-1.0,1.0),true));
       // reuse existing table
       acc.reset();
       css.reset();
       css = rwSource(false);
       newID = css->newSolutionID(60.);
       CPPUNIT_ASSERT_EQUAL(1l, newID);
       acc = css->rwSolution(newID);
       acc->setLeakage(JonesIndex(2u,1u),JonesDTerm(casacore::Complex(0.1,-0.1),true,casacore::Complex(-0.1,0.4),false));
       // once again reuse the table
       acc.reset();
       css.reset();
       css = rwSource(false);
       newID = css->newSolutionID(120.);
       CPPUNIT_ASSERT_EQUAL(2l, newID);
       acc = css->rwSolution(newID);
       acc->setBandpass(JonesIndex(1u,1u),JonesJTerm(casacore::Complex(1.0,-0.2),true,casacore::Complex(0.9,-0.1),true),1u);
       // once again reuse the table
       acc.reset();
       css.reset();
       css = rwSource(false);
       newID = css->newSolutionID(180.);
       CPPUNIT_ASSERT_EQUAL(3l, newID);
       acc = css->rwSolution(newID);
       acc->setBPLeakage(JonesIndex(1u,1u),JonesDTerm(casacore::Complex(0.1,-0.2),true,casacore::Complex(-0.1,-0.1),true),1u);
       CPPUNIT_ASSERT_EQUAL(css->solutionIDBefore(180.0).first, newID);
       CPPUNIT_ASSERT_EQUAL(css->solutionIDAfter(120.0).first, newID);
       CPPUNIT_ASSERT_EQUAL(css->solutionIDBefore(180.0).second, 180.);
       CPPUNIT_ASSERT_EQUAL(css->solutionIDAfter(120.0).second, 180.);
   }

   // common code testing leakages and gains in the test table
   // it is used in both normal source/accessor test
   void doGainAndLeakageTest(const boost::shared_ptr<ICalSolutionConstAccessor> &acc) {
       CPPUNIT_ASSERT(acc);
       // test gains
       for (casacore::uInt ant = 0; ant<6; ++ant) {
            for (casacore::uInt beam = 0; beam<3; ++beam) {
                 const JonesJTerm gain = acc->gain(JonesIndex(ant,beam));
                 if ((ant == 0) && (beam == 0)) {
                     testComplex(casacore::Complex(1.,-1.), gain.g1());
                     testComplex(casacore::Complex(-1.,1.), gain.g2());
                     CPPUNIT_ASSERT(gain.g1IsValid());
                     CPPUNIT_ASSERT(gain.g2IsValid());
                 } else {
                     // default gain is 1.0
                     testComplex(casacore::Complex(1.0,0.), gain.g1());
                     testComplex(casacore::Complex(1.0,0.), gain.g2());
                     CPPUNIT_ASSERT(!gain.g1IsValid());
                     CPPUNIT_ASSERT(!gain.g2IsValid());
                 }
            }
       }
       // test leakages
       for (casacore::uInt ant = 0; ant<6; ++ant) {
            for (casacore::uInt beam = 0; beam<3; ++beam) {
                 const JonesDTerm leakage = acc->leakage(JonesIndex(ant,beam));
                 if ((ant == 2) && (beam == 1)) {
                     testComplex(casacore::Complex(0.1,-0.1), leakage.d12());
                     testComplex(casacore::Complex(-0.1,0.4), leakage.d21());
                     CPPUNIT_ASSERT(leakage.d12IsValid());
                     CPPUNIT_ASSERT(!leakage.d21IsValid());
                 } else {
                     // default leakage is 0.0
                     testComplex(casacore::Complex(0.), leakage.d12());
                     testComplex(casacore::Complex(0.), leakage.d21());
                     CPPUNIT_ASSERT(!leakage.d12IsValid());
                     CPPUNIT_ASSERT(!leakage.d21IsValid());
                 }
            }
       }
   }
   // common code testing bandpass calibration data
   // it is used in both normal source/accessor test
   void doBandpassTest(const boost::shared_ptr<ICalSolutionConstAccessor> &acc) {
       CPPUNIT_ASSERT(acc);
       for (casa::uInt ant = 0; ant<6; ++ant) {
            for (casa::uInt beam = 0; beam<3; ++beam) {
                 const JonesIndex index(ant,beam);
                 for (casa::uInt chan = 0; chan < 8; ++chan) {
                      const JonesJTerm bp = acc->bandpass(index,chan);
                      if ((ant == 1) && (beam == 1) && (chan == 1)) {
                          testComplex(casa::Complex(1.0,-0.2), bp.g1());
                          testComplex(casa::Complex(0.9,-0.1), bp.g2());
                          CPPUNIT_ASSERT(bp.g1IsValid());
                          CPPUNIT_ASSERT(bp.g2IsValid());
                      } else {
                          // default bandpass gain is 1.0
                          testComplex(casa::Complex(1.0,0.), bp.g1());
                          testComplex(casa::Complex(1.0,0.), bp.g2());
                          CPPUNIT_ASSERT(!bp.g1IsValid());
                          CPPUNIT_ASSERT(!bp.g2IsValid());
                      }
                 }
            }
       }
       // test leakages
       for (casa::uInt ant = 0; ant<6; ++ant) {
            for (casa::uInt beam = 0; beam<3; ++beam) {
                 const JonesIndex index(ant,beam);
                 for (casa::uInt chan = 0; chan < 8; ++chan) {
                      const JonesDTerm bpl = acc->bpleakage(index,chan);
                      if ((ant == 1) && (beam == 1) && (chan == 1)) {
                          testComplex(casa::Complex(0.1,-0.2), bpl.d12());
                          testComplex(casa::Complex(-0.1,-0.1), bpl.d21());
                          CPPUNIT_ASSERT(bpl.d12IsValid());
                          CPPUNIT_ASSERT(bpl.d21IsValid());
                      } else {
                          // default bandpass leakage is 0.0
                          testComplex(casa::Complex(0.), bpl.d12());
                          testComplex(casa::Complex(0.), bpl.d21());
                          CPPUNIT_ASSERT(!bpl.d12IsValid());
                          CPPUNIT_ASSERT(!bpl.d21IsValid());
                      }
                 }
            }
       }
   }

   void testTrailingBlankEntry() {
       // reuse generation code which initialises 4 entries with all available products between them
       testCreate();
       const boost::shared_ptr<ICalSolutionSource> css = rwSource(false);
       CPPUNIT_ASSERT(css);
       const long newID = css->newSolutionID(240.);
       CPPUNIT_ASSERT_EQUAL(4l, newID);

       // reading as most recent solution
       {
          const long sID = css->mostRecentSolution();
          CPPUNIT_ASSERT_EQUAL(newID, sID);
          const boost::shared_ptr<ICalSolutionConstAccessor> acc = css->roSolution(sID);
          CPPUNIT_ASSERT(acc);
          doGainAndLeakageTest(acc);
          doBandpassTest(acc);
       }
       // reading by giving solution ID directly
       {
          const boost::shared_ptr<ICalSolutionConstAccessor> acc = css->roSolution(newID);
          CPPUNIT_ASSERT(acc);
          doGainAndLeakageTest(acc);
          doBandpassTest(acc);
       }
   }

   void testBlankEntries() {
       // rerun the code creating a table to ensure we always get the same starting point in the spirit of unit tests
       testCreate();
       // although not strictly necessary, run the following code inside the block to ensure destructors are called
       {
          const boost::shared_ptr<ICalSolutionSource> css = rwSource(false);
          for (long id = 4; id < 10; ++id) {
               const long newID = css->newSolutionID(60. * id);
               CPPUNIT_ASSERT_EQUAL(id, newID);
               // deliberatly don't set any calibration information for this solution ID
          }

          const long newID = css->newSolutionID(600.);
          CPPUNIT_ASSERT_EQUAL(10l, newID);
          boost::shared_ptr<ICalSolutionAccessor> acc = css->rwSolution(newID);
          acc->setGain(JonesIndex(0u,0u),JonesJTerm(casa::Complex(1.0,-1.0),true,casa::Complex(-1.0,1.0),true));
          acc->setLeakage(JonesIndex(2u,1u),JonesDTerm(casa::Complex(0.1,-0.1),true,casa::Complex(-0.1,0.4),false));
          acc->setBandpass(JonesIndex(1u,1u),JonesJTerm(casa::Complex(1.0,-0.2),true,casa::Complex(0.9,-0.1),true),1u);
          acc->setBPLeakage(JonesIndex(1u,1u),JonesDTerm(casa::Complex(0.1,-0.2),true,casa::Complex(-0.1,-0.1),true),1u);

       }
       // reading
       const boost::shared_ptr<ICalSolutionConstSource> css = roSource();
       CPPUNIT_ASSERT(css);

       {
          const long sID = css->mostRecentSolution();
          CPPUNIT_ASSERT_EQUAL(10l, sID);
          const boost::shared_ptr<ICalSolutionConstAccessor> acc = css->roSolution(sID);
          CPPUNIT_ASSERT(acc);
          doGainAndLeakageTest(acc);
          doBandpassTest(acc);
       }

       // rows with empty cells
       for (long id = 9; id >= 4; --id) {
            const long sID = css->solutionID(60. * id);
            CPPUNIT_ASSERT_EQUAL(id, sID);
            const boost::shared_ptr<ICalSolutionConstAccessor> acc = css->roSolution(sID);
            CPPUNIT_ASSERT(acc);
            doGainAndLeakageTest(acc);
            doBandpassTest(acc);
       }
   }

   void testRead() {
       // rerun the code creating a table, although we could've just relied on the fact that testCreate() is executed
       // just before this test
       testCreate();
       const boost::shared_ptr<ICalSolutionConstSource> css = roSource();
       CPPUNIT_ASSERT(css);
       const long sID = css->mostRecentSolution();
       CPPUNIT_ASSERT_EQUAL(3l, sID);
       for (long id = 0; id<4; ++id) {
            CPPUNIT_ASSERT_EQUAL(id, css->solutionID(0.5+60.*id));
            CPPUNIT_ASSERT_EQUAL(id, css->solutionIDBefore(0.5+60.*id).first);
            CPPUNIT_ASSERT(abs(60.*id - css->solutionIDBefore(0.5+60.*id).second) < 0.01);
            if (id > 0) {
                CPPUNIT_ASSERT_EQUAL(id, css->solutionIDAfter(0.5+60.*(id-1)).first);
                CPPUNIT_ASSERT(abs(60.*id - css->solutionIDAfter(0.5+60.*(id-1)).second) < 0.01);
            }
       }
       const boost::shared_ptr<ICalSolutionConstAccessor> acc = css->roSolution(sID);
       CPPUNIT_ASSERT(acc);
       doGainAndLeakageTest(acc);

       // test bandpasses
       doBandpassTest(acc);
   }

   void testDelayedWrite() {
       // do essentially the same as testRead case but write the table in a different fashion - first request IDs for all entries then write.
       // this test would've caught read-write vs. read-only bug we lived with for a while (see ASKAPSDP-3731)

       const boost::shared_ptr<ICalSolutionSource> css = rwSource(true);
       CPPUNIT_ASSERT(css);
       long ids[4] = {-1l,-1l,-1l,-1l};
       for (int row = 0; row<4; ++row) {
            ids[row] = css->newSolutionID(60. * row);
            CPPUNIT_ASSERT_EQUAL(static_cast<long>(row), ids[row]);
       }
       boost::shared_ptr<ICalSolutionAccessor> acc = css->rwSolution(ids[0]);
       CPPUNIT_ASSERT(acc);
       acc->setGain(JonesIndex(0u,0u),JonesJTerm(casa::Complex(1.0,-1.0),true,casa::Complex(-1.0,1.0),true));
       acc = css->rwSolution(ids[1]);
       acc->setLeakage(JonesIndex(2u,1u),JonesDTerm(casa::Complex(0.1,-0.1),true,casa::Complex(-0.1,0.4),false));
       acc = css->rwSolution(ids[2]);
       acc->setBandpass(JonesIndex(1u,1u),JonesJTerm(casa::Complex(1.0,-0.2),true,casa::Complex(0.9,-0.1),true),1u);
       acc = css->rwSolution(ids[3]);
       acc->setBPLeakage(JonesIndex(1u,1u),JonesDTerm(casa::Complex(0.1,-0.2),true,casa::Complex(-0.1,-0.1),true),1u);
       acc.reset();

       // now test the content  first try read-write accesspr
       const long sID = css->mostRecentSolution();
       CPPUNIT_ASSERT_EQUAL(3l, sID);
       for (long id = 0; id<4; ++id) {
            CPPUNIT_ASSERT_EQUAL(id, css->solutionID(0.5+60.*id));
       }
       boost::shared_ptr<ICalSolutionConstAccessor> accRO = css->roSolution(sID);
       CPPUNIT_ASSERT(accRO);
       doGainAndLeakageTest(accRO);
       doBandpassTest(accRO);

       // now open the same table with read-only access and redo the test
       const boost::shared_ptr<ICalSolutionConstSource> cssRO = roSource();
       CPPUNIT_ASSERT_EQUAL(sID, cssRO->mostRecentSolution());
       for (long id = 0; id<4; ++id) {
            CPPUNIT_ASSERT_EQUAL(id, cssRO->solutionID(0.5+60.*id));
       }
       accRO = cssRO->roSolution(sID);
       CPPUNIT_ASSERT(accRO);
       doGainAndLeakageTest(accRO);
       doBandpassTest(accRO);
   }

   void testTooFarIntoThePast() {
       boost::shared_ptr<ICalSolutionSource> css = rwSource(true);
       CPPUNIT_ASSERT(css);
       long newID = css->newSolutionID(1000.);
       CPPUNIT_ASSERT_EQUAL(0l, newID);
       boost::shared_ptr<ICalSolutionAccessor> acc = css->rwSolution(newID);
       CPPUNIT_ASSERT(acc);
       acc->setGain(JonesIndex(0u,0u),JonesJTerm(casa::Complex(1.0,-1.0),true,casa::Complex(-1.0,1.0),true));
       newID = css->newSolutionID(1060.);
       CPPUNIT_ASSERT_EQUAL(1l, newID);
       acc = css->rwSolution(newID);
       acc->setLeakage(JonesIndex(2u,1u),JonesDTerm(casa::Complex(0.1,-0.1),true,casa::Complex(-0.1,0.4),false));
       acc.reset();
       // add empty row
       css->newSolutionID(1120.);
       // first read using time within the table
       const long id1 = css->solutionID(1000.5);
       CPPUNIT_ASSERT_EQUAL(0l, id1);
       const long id2 = css->solutionID(1060.5);
       CPPUNIT_ASSERT_EQUAL(1l, id2);
       // this should work too as we the solution remains valid until there is a new entry, i.e. forever at the end of the table
       const long id3 = css->solutionID(1200.5);
       CPPUNIT_ASSERT_EQUAL(2l, id3);
       // reading test for gains and leakages should be successful too
       boost::shared_ptr<ICalSolutionConstAccessor> accRO = css->roSolution(id3);
       doGainAndLeakageTest(accRO);
       // the following would cause an exception
       css->solutionID(990.);
   }

   void testChanAdapterRead() {
       // rerun the code creating a table, although we could've just relied on the fact that testCreate() is executed
       // just before this test
       testCreate();
       // adapter to offset everything by one channel
       const boost::shared_ptr<ICalSolutionConstSource> css(new ChanAdapterCalSolutionConstSource(roSource(),1u));
       CPPUNIT_ASSERT(css);
       const long sID = css->mostRecentSolution();
       CPPUNIT_ASSERT_EQUAL(3l, sID);
       for (long id = 0; id<4; ++id) {
            CPPUNIT_ASSERT_EQUAL(id, css->solutionID(0.5+60.*id));
       }
       const boost::shared_ptr<ICalSolutionConstAccessor> acc = css->roSolution(sID);
       CPPUNIT_ASSERT(acc);
       doGainAndLeakageTest(acc);

       // test bandpasses
       for (casacore::uInt ant = 0; ant<6; ++ant) {
            for (casacore::uInt beam = 0; beam<3; ++beam) {
                 const JonesIndex index(ant,beam);
                 for (casacore::uInt chan = 0; chan < 7; ++chan) {
                      const JonesJTerm bp = acc->bandpass(index,chan);
                      if ((ant == 1) && (beam == 1) && (chan == 0)) {
                          testComplex(casacore::Complex(1.0,-0.2), bp.g1());
                          testComplex(casacore::Complex(0.9,-0.1), bp.g2());
                          CPPUNIT_ASSERT(bp.g1IsValid());
                          CPPUNIT_ASSERT(bp.g2IsValid());
                      } else {
                          // default bandpass gain is 1.0
                          testComplex(casacore::Complex(1.0,0.), bp.g1());
                          testComplex(casacore::Complex(1.0,0.), bp.g2());
                          CPPUNIT_ASSERT(!bp.g1IsValid());
                          CPPUNIT_ASSERT(!bp.g2IsValid());
                      }
                 }
            }
       }
   }

   void testUndefinedGains() {
       const boost::shared_ptr<ICalSolutionConstAccessor> acc = accessorForExistingTable();
       CPPUNIT_ASSERT(acc);
       // only 6 antennas, 3 beams and 8 channels are defined
       acc->gain(JonesIndex(7u,0u));
   }

   void testUndefinedLeakages() {
       const boost::shared_ptr<ICalSolutionConstAccessor> acc = accessorForExistingTable();
       CPPUNIT_ASSERT(acc);
       // only 6 antennas, 3 beams and 8 channels are defined
       acc->leakage(JonesIndex(3u,3u));
   }

   void testUndefinedBandpasses() {
       const boost::shared_ptr<ICalSolutionConstAccessor> acc = accessorForExistingTable();
       CPPUNIT_ASSERT(acc);
       // only 6 antennas, 3 beams and 8 channels are defined
       acc->bandpass(JonesIndex(0u,0u),8);
   }

   void testChanAdapterUndefinedBandpass() {
       // rerun the code creating a table, although we could've just relied on the fact that testCreate() is executed
       // just before this test
       testCreate();
       // adapter to offset everything by one channel
       const boost::shared_ptr<ICalSolutionConstSource> css(new ChanAdapterCalSolutionConstSource(roSource(),1u));
       CPPUNIT_ASSERT(css);
       const long sID = css->mostRecentSolution();
       CPPUNIT_ASSERT_EQUAL(3l, sID);
       const boost::shared_ptr<ICalSolutionConstAccessor> acc = css->roSolution(sID);
       CPPUNIT_ASSERT(acc);

       // only 6 antennas, 3 beams and 8 channels are defined, after offseting by one channel we got 7 channels only
       acc->bandpass(JonesIndex(0u,0u),7);
   }

   void testUndefinedSolution() {
       const boost::shared_ptr<ICalSolutionConstSource> css = roSource();
       CPPUNIT_ASSERT(css);
       const long id = css->solutionID(0.5);
       CPPUNIT_ASSERT_EQUAL(0l, id);
       const boost::shared_ptr<ICalSolutionConstAccessor> acc = css->roSolution(id);
       CPPUNIT_ASSERT(acc);
       try {
          // the following should be successful because the first solution in the table
          // was the gain solution
          const JonesJTerm gain = acc->gain(JonesIndex(0u,0u));
          testComplex(casacore::Complex(1.,-1.), gain.g1());
          testComplex(casacore::Complex(-1.,1.), gain.g2());
          CPPUNIT_ASSERT(gain.g1IsValid());
          CPPUNIT_ASSERT(gain.g2IsValid());
       }
       catch (const AskapError &) {
          // this shouldn't have happened
          CPPUNIT_ASSERT(false);
       }
       // the following should cause an exception because no leakage is defined at or before row 0.
       acc->leakage(JonesIndex(0u,0u));
   }

   void testCreateManyRows() {
       boost::shared_ptr<ICalSolutionSource> css = rwSource(true);
       const size_t nSolutions = 10;
       for (size_t sol = 0; sol<nSolutions; ++sol) {
            const long newID = css->newSolutionID(55553.*86400+3600.*double(sol));
            CPPUNIT_ASSERT_EQUAL(long(sol), newID);
            boost::shared_ptr<ICalSolutionAccessor> acc = css->rwSolution(newID);
            for (casacore::Short ant = 0; ant<6; ++ant) {
                 for (casacore::Short beam=0; beam<3; ++beam) {
                      const float amp = float(ant)/6. + 0.5 + 0.1*float(sol);
                      const float phase = casacore::C::pi/3.*float(beam);
                      const casacore::Complex val = casacore::Complex(sin(phase),cos(phase))*amp;
                      acc->setGain(JonesIndex(ant,beam),JonesJTerm(val,true,-val,true));
                 }
            }
       }
       // reuse solution source
       css.reset();
       css = rwSource(false);
       for (size_t sol = 0; sol<nSolutions; ++sol) {
            const long id = css->solutionID(55553.*86400+3600.*double(sol)+0.5);
            CPPUNIT_ASSERT_EQUAL(long(sol), id);
            const boost::shared_ptr<ICalSolutionConstAccessor> acc = css->roSolution(id);
            CPPUNIT_ASSERT(acc);
            for (casacore::Short ant = 0; ant<6; ++ant) {
                 for (casacore::Short beam=0; beam<3; ++beam) {
                      const float amp = float(ant)/6. + 0.5 + 0.1*float(sol);
                      const float phase = casacore::C::pi/3.*float(beam);
                      const casacore::Complex val = casacore::Complex(sin(phase),cos(phase))*amp;
                      const JonesJTerm gain = acc->gain(JonesIndex(ant,beam));
                      testComplex(val, gain.g1());
                      testComplex(-val, gain.g2());
                      CPPUNIT_ASSERT(gain.g1IsValid());
                      CPPUNIT_ASSERT(gain.g2IsValid());
                 }
            }
       }
   }
   void testCreateManyRows1() {
       boost::shared_ptr<ICalSolutionSource> css = rwSource(true);
       const size_t nSolutions = 10;
       for (size_t sol = 0; sol<nSolutions; ++sol) {
            const long newID = css->newSolutionID(55553.*86400+3600.*double(sol));
            CPPUNIT_ASSERT_EQUAL(long(sol), newID);
            boost::shared_ptr<ICalSolutionAccessor> acc = css->rwSolution(newID);
            ASKAPDEBUGASSERT(acc);
            for (casacore::uInt ant = 0; ant<6; ++ant) {
                 for (casacore::uInt beam=0; beam<3; ++beam) {
                      const float amp = float(ant)/6. + 0.5 + 0.1*float(sol);
                      const float phase = casacore::C::pi/3.*float(beam);
                      const casacore::Complex val = casacore::Complex(sin(phase),cos(phase))*amp;
                      const accessors::JonesIndex index(ant,beam);
                      acc->setJonesElement(index, casacore::Stokes::XX, val);
                      acc->setJonesElement(index, casacore::Stokes::YY, -val);
                      acc->setJonesElement(index, casacore::Stokes::XY, val*0.1f);
                      acc->setJonesElement(index, casacore::Stokes::YX, -val*0.1f);
                 }
            }
       }
       // reuse solution source
       css.reset();
       css = rwSource(false);
       for (size_t sol = 0; sol<nSolutions; ++sol) {
            const long id = css->solutionID(55553.*86400+3600.*double(sol)+0.5);
            CPPUNIT_ASSERT_EQUAL(long(sol), id);
            const boost::shared_ptr<ICalSolutionConstAccessor> acc = css->roSolution(id);
            CPPUNIT_ASSERT(acc);
            for (casacore::Short ant = 0; ant<6; ++ant) {
                 for (casacore::Short beam=0; beam<3; ++beam) {
                      const float amp = float(ant)/6. + 0.5 + 0.1*float(sol);
                      const float phase = casacore::C::pi/3.*float(beam);
                      const casacore::Complex val = casacore::Complex(sin(phase),cos(phase))*amp;
                      const JonesJTerm gain = acc->gain(JonesIndex(ant,beam));
                      testComplex(val, gain.g1());
                      testComplex(-val, gain.g2());
                      CPPUNIT_ASSERT(gain.g1IsValid());
                      CPPUNIT_ASSERT(gain.g2IsValid());
                      const JonesDTerm leakage = acc->leakage(JonesIndex(ant,beam));
                      testComplex(val*0.1f, leakage.d12());
                      testComplex(-val*0.1f, leakage.d21());
                      CPPUNIT_ASSERT(leakage.d12IsValid());
                      CPPUNIT_ASSERT(leakage.d21IsValid());
                 }
            }
       }
   }

};

} // namespace accessors

} // namespace askap
