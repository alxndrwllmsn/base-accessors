/// @file GainFilterAdapterCalSolutionTest.h
///
/// @brief Unit tests for gain calibration normalisation
/// @details The tests gathered in this file normalise gains and
/// check the results
///
/// @copyright (c) 2025 CSIRO
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
/// @author Mark Wieringa <mark.wieringa@csiro.au>


#include <askap/calibaccess/GainFilterAdapterCalSolutionSource.h>
#include <askap/calibaccess/TableCalSolutionSource.h>


#include <Common/ParameterSet.h>
#include <boost/shared_ptr.hpp>
#include <cppunit/extensions/HelperMacros.h>

namespace askap {
namespace accessors {
    using namespace casacore;

    class GainFilterTest : public CppUnit::TestFixture
    {
        CPPUNIT_TEST_SUITE(GainFilterTest);
        CPPUNIT_TEST(testGainNormalisationNone);
        CPPUNIT_TEST(testGainNormalisationMean);
        CPPUNIT_TEST(testGainNormalisationMeanAllBeams);
        CPPUNIT_TEST(testGainNormalisationMedian);
        CPPUNIT_TEST_SUITE_END();

        public:
        void setUp() {
            itsCalTable = "calibdata.tab";
        };

        void tearDown() {
        }

        boost::shared_ptr<ICalSolutionSource> rwSource(LOFAR::ParameterSet& parset, bool doRemove, uInt nAnt, uInt nBeam, uInt nChan) {
            if (doRemove) {
                TableCalSolutionSource::removeOldTable(itsCalTable);
            }
            boost::shared_ptr<TableCalSolutionSource> tcss(new TableCalSolutionSource(itsCalTable, nAnt, nBeam, nChan));
            CPPUNIT_ASSERT(tcss);
            boost::shared_ptr<ICalSolutionSource> css(new GainFilterAdapterCalSolutionSource(tcss,parset));
            return css;
        }
        boost::shared_ptr<ICalSolutionSource> roSource(uInt nAnt, uInt nBeam, uInt nChan) {
            boost::shared_ptr<TableCalSolutionSource> css(new TableCalSolutionSource(itsCalTable, nAnt, nBeam, nChan));
            CPPUNIT_ASSERT(css);
            return css;
        }

        void fillGains(LOFAR::ParameterSet& parset, uInt nAnt, uInt nBeam, uInt nTime) {
            // set up gains for solutions at t=0, t=60, ...
            const uInt nChan = 1;
            boost::shared_ptr<ICalSolutionSource> css = rwSource(parset, true, nAnt, nBeam, nChan);
            // store in Cube for easy access, 2 pols and time share first index
            itsGains.resize(2*nTime,nAnt,nBeam);
            for (long t = 0; t < nTime; t++) {
                long newID = css->newSolutionID(t*60.);
                CPPUNIT_ASSERT_EQUAL(t, newID);
                boost::shared_ptr<ICalSolutionAccessor> acc = css->rwSolution(newID);
                const float pi4 = C::pi/4;
                for (uInt j = 0; j< nBeam; j++) {
                    for (uInt i = 0; i < nAnt; i++) {
                        long index = 2 * newID;
                        itsGains(index + 0,i,j) = (1.1f+j/10.f-t/5.f)*exp(casacore::Complex(0.,pi4*(-1.+2*(i%2))));
                        itsGains(index + 1,i,j) = (1.2f-t/5.f)*exp(casacore::Complex(0.,pi4/2));
                        acc->setGain(JonesIndex(i,j),JonesJTerm(itsGains(index,i,j),true,itsGains(index+1,i,j),true));
                    }
                }
            }
        }

        Cube<Complex> readGains(uInt nAnt, uInt nBeam, uInt nTime) {
            const uInt nChan = 1;
            boost::shared_ptr<ICalSolutionConstSource> css = roSource(nAnt, nBeam, nChan);
            // store in Cube for easy access, 2 pols and time share first index
            Cube<Complex> gains(2*nTime,nAnt,nBeam);
            for (long t = 0; t < nTime; t++) {
                long newID = css->solutionID(t*60.);
                CPPUNIT_ASSERT_EQUAL(t, newID);
                boost::shared_ptr<ICalSolutionConstAccessor> acc = css->roSolution(newID);
                const float pi4 = 3.14159265358/4;
                for (uInt j = 0; j< nBeam; j++) {
                    for (uInt i = 0; i < nAnt; i++) {
                        long index = 2 * newID;
                        JonesJTerm jj = acc->gain(JonesIndex(i,j));
                        gains(index,i,j) = jj.g1();
                        gains(index+1,i,j) = jj.g2();
                    }
                }
            }
            return gains;
        }


        void testGainNormalisationNone() {
            testGainNormalisation("none");
        }
        void testGainNormalisationMean() {
            testGainNormalisation("mean");
        }
        void testGainNormalisationMeanAllBeams() {
            testGainNormalisation("mean",true);
        }
        void testGainNormalisationMedian() {
            testGainNormalisation("median");
        }

        void testGainNormalisation(const std::string& fluxnorm, bool allbeams=false) {
            const int nAnt = 36;
            const int nBeam = 2;
            const int nTime = 2;
            const int nChan = 1;

            // set up ParameterSet
            LOFAR::ParameterSet parset;
            parset.add("fluxnorm",fluxnorm);
            parset.replace(LOFAR::KVpair("fluxnorm.allbeams",allbeams));
            parset.add("calibaccess","table");
            parset.add("calibaccess.table",itsCalTable);
            parset.replace(LOFAR::KVpair("calibaccess.table.maxant",nAnt));
            parset.replace(LOFAR::KVpair("calibaccess.table.maxbeam",nBeam));
            parset.replace(LOFAR::KVpair("calibaccess.table.maxchan",nChan));
            parset.replace(LOFAR::KVpair("calibaccess.table.reuse",true));

            fillGains(parset, nAnt,nBeam,nTime);

            Cube<Complex> gains = readGains(nAnt,nBeam,nTime);
            for (int beam = 0; beam < nBeam; beam++) {
                if (fluxnorm == "" || fluxnorm == "none") {
                    // check gains are unchanged
                    CPPUNIT_ASSERT(allEQ(itsGains(Slice(),Slice(),Slice(beam)),gains(Slice(),Slice(),Slice(beam))));
                } else if (fluxnorm == "mean") {
                    // check mean gains are 1
                    if (beam == 0 || allbeams) {
                        CPPUNIT_ASSERT_DOUBLES_EQUAL(sqrt(mean(square(amplitude(gains(Slice(0,nTime,2),Slice(),Slice(beam)))))),1,1e-5);
                        CPPUNIT_ASSERT_DOUBLES_EQUAL(sqrt(mean(square(amplitude(gains(Slice(1,nTime,2),Slice(),Slice(beam)))))),1,1e-5);
                    } else {
                        // check gains are unchanged
                        CPPUNIT_ASSERT(allEQ(itsGains(Slice(),Slice(),Slice(beam)),gains(Slice(),Slice(),Slice(beam))));
                    }
                } else if (fluxnorm == "median") {
                    if (beam == 0 || allbeams) {
                        // check median gains are 1
                        CPPUNIT_ASSERT_DOUBLES_EQUAL(sqrt(median(square(amplitude(gains(Slice(0,nTime,2),Slice(),Slice(beam)))))),1,1e-5);
                        CPPUNIT_ASSERT_DOUBLES_EQUAL(sqrt(median(square(amplitude(gains(Slice(1,nTime,2),Slice(),Slice(beam)))))),1,1e-5);
                    } else {
                        // check gains are unchanged
                        CPPUNIT_ASSERT(allEQ(itsGains(Slice(),Slice(),Slice(beam)),gains(Slice(),Slice(),Slice(beam))));
                    }
                }
            }
        }

     private:
        Cube<Complex> itsGains;
        string itsCalTable;

   };
  } // namespace synthesis

} // namespace askap
