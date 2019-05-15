/// @file
///
/// Unit test for the parset-based implementation of the interface to access
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
#include <askap/calibaccess/ParsetCalSolutionSource.h>
#include <askap/calibaccess/ParsetCalSolutionAccessor.h>
#include <askap/calibaccess/JonesIndex.h>

#include <boost/shared_ptr.hpp>


namespace askap {

namespace accessors {

class ParsetCalSolutionTest : public CppUnit::TestFixture 
{
   CPPUNIT_TEST_SUITE(ParsetCalSolutionTest);
   CPPUNIT_TEST(testReadWrite);
   CPPUNIT_TEST(testOverwrite);
   CPPUNIT_TEST(testPartiallyUndefined);
   CPPUNIT_TEST(testSolutionSource);
   CPPUNIT_TEST_SUITE_END();
protected:
   static void createDummyParset(ICalSolutionAccessor &acc) {
       for (casacore::uInt ant=0; ant<5; ++ant) {
            for (casacore::uInt beam=0; beam<4; ++beam) { 
                 const float tag = float(ant)/100. + float(beam)/1000.;
                 acc.setJonesElement(ant,beam,casacore::Stokes::XX,casacore::Complex(1.1+tag,0.1));
                 acc.setJonesElement(ant,beam,casacore::Stokes::YY,casacore::Complex(1.1,-0.1-tag));
                 acc.setJonesElement(ant,beam,casacore::Stokes::XY,casacore::Complex(0.1+tag,-0.1));
                 acc.setJonesElement(ant,beam,casacore::Stokes::YX,casacore::Complex(-0.1,0.1+tag));

                 for (casacore::uInt chan=0; chan<20; ++chan) {
                     acc.setBandpassElement(ant,beam,casacore::Stokes::XX,chan,casacore::Complex(1.,0.));
                     acc.setBandpassElement(ant,beam,casacore::Stokes::YY,chan,casacore::Complex(1.,0.));
                 }
            }
       }
   }
   
   static void createDummyParset(const std::string &fname) {
       ParsetCalSolutionAccessor acc(fname);
       createDummyParset(acc);
   }

   static void testComplex(const casacore::Complex &expected, const casacore::Complex &obtained, const float tol = 1e-5) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(real(expected),real(obtained),tol);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(imag(expected),imag(obtained),tol);
   }

   static void testDummyParset(const ICalSolutionConstAccessor &acc) {
        for (casacore::uInt ant=0; ant<5; ++ant) {
            for (casacore::uInt beam=0; beam<4; ++beam) {
                 CPPUNIT_ASSERT(acc.jonesValid(ant,beam,0));
                 const casacore::SquareMatrix<casacore::Complex, 2> jones = acc.jones(ant,beam,0);
                 const float tag = float(ant)/100. + float(beam)/1000.;
                 testComplex(casacore::Complex(1.1+tag,0.1), jones(0,0));
                 testComplex(casacore::Complex(1.1,-0.1-tag), jones(1,1));
                 testComplex(casacore::Complex(0.1+tag,-0.1) * casacore::Complex(1.1+tag,0.1), jones(0,1));
                 testComplex(casacore::Complex(-0.1,0.1+tag) * casacore::Complex(1.1,-0.1-tag), -jones(1,0));

                 const JonesIndex index(ant,beam);
                 CPPUNIT_ASSERT(index.antenna() == casacore::Short(ant));
                 CPPUNIT_ASSERT(index.beam() == casacore::Short(beam));

                 const casacore::SquareMatrix<casacore::Complex, 2> jones2 = acc.jones(index,10);
                 testComplex(casacore::Complex(1.1+tag,0.1), jones2(0,0));
                 testComplex(casacore::Complex(1.1,-0.1-tag), jones2(1,1));
                 testComplex(casacore::Complex(0.1+tag,-0.1) * casacore::Complex(1.1+tag,0.1), jones2(0,1));
                 testComplex(casacore::Complex(-0.1,0.1+tag) * casacore::Complex(1.1,-0.1-tag), -jones2(1,0));

                 const JonesJTerm jTerm = acc.gain(index);
                 CPPUNIT_ASSERT(jTerm.g1IsValid() && jTerm.g2IsValid());
                 testComplex(casacore::Complex(1.1+tag,0.1), jTerm.g1());
                 testComplex(casacore::Complex(1.1,-0.1-tag), jTerm.g2());

                 const JonesDTerm dTerm = acc.leakage(index);
                 CPPUNIT_ASSERT(dTerm.d12IsValid() && dTerm.d21IsValid());
                 testComplex(casacore::Complex(0.1+tag,-0.1), dTerm.d12());
                 testComplex(casacore::Complex(-0.1,0.1+tag), dTerm.d21()); 

                 for (casacore::uInt chan=0; chan<20; ++chan) {
                      const JonesJTerm bpTerm = acc.bandpass(index, chan);
                      CPPUNIT_ASSERT(bpTerm.g1IsValid() && bpTerm.g2IsValid());
                      testComplex(casacore::Complex(1.,0.), bpTerm.g1());
                      testComplex(casacore::Complex(1.,0.), bpTerm.g2());
                 }
            }
        }
   }

public:
   void testReadWrite() {
        const std::string fname = "tmp.testparset";
        createDummyParset(fname);
        ParsetCalSolutionAccessor acc(fname);
        testDummyParset(acc);
   }

   void testOverwrite() {
        const std::string fname = "tmp.testparset";
        createDummyParset(fname);
        {
          // now write again and overwrite the first antenna/beam only
          // actual write happens in destructor, hence the curly brackets
          ParsetCalSolutionAccessor acc(fname);
          acc.setJonesElement(0,0,casacore::Stokes::XX,casacore::Complex(1.1,0.1));
          acc.setJonesElement(0,0,casacore::Stokes::YY,casacore::Complex(1.05,-0.1));
          acc.setJonesElement(0,0,casacore::Stokes::XY,casacore::Complex(0.13,-0.12));
          acc.setJonesElement(0,0,casacore::Stokes::YX,casacore::Complex(-0.14,0.11));

          // Write bandpass for the first channel/antenna/beam.
          acc.setBandpassElement(0,0,casacore::Stokes::XX,0,casacore::Complex(1.,0.));
          acc.setBandpassElement(0,0,casacore::Stokes::YY,0,casacore::Complex(1.,0.));
        }
        // now read
        ParsetCalSolutionAccessor acc(fname);
        for (casacore::uInt ant=0; ant<10; ++ant) {
            for (casacore::uInt beam=0; beam<6; ++beam) {
                 CPPUNIT_ASSERT_EQUAL(!ant && !beam, acc.jonesValid(ant,beam,0));
                 const JonesIndex index(ant,beam);
                 CPPUNIT_ASSERT(index.antenna() == casacore::Short(ant));
                 CPPUNIT_ASSERT(index.beam() == casacore::Short(beam));
                 const casacore::SquareMatrix<casacore::Complex, 2> jones = acc.jones(index,0);
                 if (!ant && !beam) {
                     testComplex(casacore::Complex(1.1,0.1), jones(0,0));
                     testComplex(casacore::Complex(1.05,-0.1), jones(1,1));
                     testComplex(casacore::Complex(0.13,-0.12) * casacore::Complex(1.1,0.1), jones(0,1));
                     testComplex(casacore::Complex(-0.14,0.11) * casacore::Complex(1.05,-0.1), -jones(1,0));
                 } else {
                     // expect default values for undefined gains/leakages
                     testComplex(casacore::Complex(1.,0.), jones(0,0));
                     testComplex(casacore::Complex(1.,0.), jones(1,1));
                     testComplex(casacore::Complex(0.,0.), jones(0,1));
                     testComplex(casacore::Complex(0.,0.), -jones(1,0));
                 }
            }
        }
   }

   void testPartiallyUndefined() {
        const std::string fname = "tmp.testparset";
        const JonesIndex index(0u,0u); 
        {
          // actual write happens in destructor, hence the curly brackets
          ParsetCalSolutionAccessor acc(fname);
          const JonesJTerm gains(casacore::Complex(1.1,0.1), true, casacore::Complex(1.05,-0.1), false);
          acc.setGain(index, gains);
          const JonesDTerm leakages(casacore::Complex(0.13,-0.12),false, casacore::Complex(-0.14,0.11), true);
          acc.setLeakage(index, leakages);
        }
        // now read and check
        ParsetCalSolutionAccessor acc(fname);
        CPPUNIT_ASSERT_EQUAL(false, acc.jonesValid(index,0));
        const casacore::SquareMatrix<casacore::Complex, 2> jones = acc.jones(index,0);

        testComplex(casacore::Complex(1.1,0.1), jones(0,0));
        // undefined gain is one
        testComplex(casacore::Complex(1.0,0.), jones(1,1));
        // undefined leakage is zero
        testComplex(casacore::Complex(0.,0.), jones(0,1));
        testComplex(casacore::Complex(-0.14,0.11), -jones(1,0));
   }

   void testSolutionSource() {
        const std::string fname = "tmp.testparset";
        ParsetCalSolutionSource ss(fname);
        long id = ss.newSolutionID(0.);
        boost::shared_ptr<ICalSolutionAccessor> rwAcc = ss.rwSolution(id);
        CPPUNIT_ASSERT(rwAcc);
        createDummyParset(*rwAcc);
        CPPUNIT_ASSERT_EQUAL(id, ss.mostRecentSolution());
        CPPUNIT_ASSERT_EQUAL(id, ss.solutionID(1e-6));
        boost::shared_ptr<ICalSolutionConstAccessor> roAcc = ss.roSolution(id);
        CPPUNIT_ASSERT(roAcc);
        testDummyParset(*roAcc);
   }
};

} // namespace accessors

} // namespace askap

