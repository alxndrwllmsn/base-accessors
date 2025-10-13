//
// @file tVerifyUVW.cc : read UVWs from the given measurement set via the 
//                       standard accessor interface, check them vs. predicted
//                       values given times, array layout and phasing info
//
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


#include <askap/dataaccess/TableDataSource.h>
#include <askap_accessors.h>
#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".tVerifyUVW");

#include <askap/askap/AskapError.h>
#include <askap/dataaccess/SharedIter.h>
#include <askap/dataaccess/ParsetInterface.h>

#include <askap/dataaccess/TableManager.h>
#include <askap/dataaccess/IDataConverterImpl.h>

// casa
#include <casacore/measures/Measures/MFrequency.h>
#include <casacore/tables/Tables/Table.h>
#include <casacore/casa/OS/Timer.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>
#include <casacore/measures/Measures/UVWMachine.h>

// boost
#include <boost/noncopyable.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>

// std
#include <stdexcept>
#include <iostream>

using namespace askap;
using namespace accessors;
using namespace boost::accumulators;

class UVWChecker : public boost::noncopyable {
public:
    /// @brief constructor
    /// @param[in] ds table data source to work with
    /// @param[in] verbose if true, more details are printed in the log
    /// @note other constructors can be written later (i.e. for some streaming applications which cannot access table-specific info)
    explicit UVWChecker(const TableConstDataSource &ds, bool verbose = false);

    /// @brief iterate over the data - main entry point
    void run();

    /// @brief method predicting UVWs for the given accessor
    /// @details For each row a UVW vector is calculated
    /// @param[in] acc const accessor to work with (only metadata are used)
    /// @return a vector with simulated UVWs, one item per row
    casacore::Vector<casacore::RigidVector<casacore::Double, 3> > simulateUVW(const IConstDataAccessor &acc) const;

private:
    /// @brief const reference data source - it is set in the constructor only
    const IConstDataSource& itsDataSource;

    /// @brief buffer with antenna positions
    std::vector<casacore::MPosition> itsLayout;

    /// @brief reference MJD for time in the accessor
    const double itsRefMJD;

    /// @brief if true, more details are printed
    const bool itsVerbose;
};

/// @brief constructor
/// @param[in] ds table data source to work with
/// @param[in] verbose if true, more details are printed in the log
UVWChecker::UVWChecker(const TableConstDataSource &ds, bool verbose) : itsDataSource(ds), itsRefMJD(59000.), itsVerbose(verbose)
{
   // all operations specific to the table-based accessor are confined to this constructor
   // the rest of the work can proceed through the general interface (hence the type of 
   // itsDataSource in this class w.r.t. the type of ds)
   const casacore::uInt nAnt = ds.getNumberOfAntennas();
   itsLayout.reserve(nAnt);
   for (casacore::uInt ant = 0; ant < nAnt; ++ant) {
        itsLayout.push_back(ds.getAntennaPosition(ant));
   }
}

void UVWChecker::run() {
  IDataSelectorPtr sel=itsDataSource.createSelector();
  IDataConverterPtr conv=itsDataSource.createConverter();  
  conv->setFrequencyFrame(casacore::MFrequency::Ref(casacore::MFrequency::TOPO),"MHz");
  conv->setEpochFrame(casacore::MEpoch(casacore::Quantity(itsRefMJD,"d"),
                      casacore::MEpoch::Ref(casacore::MEpoch::UTC)),"s");
  conv->setDirectionFrame(casacore::MDirection::Ref(casacore::MDirection::J2000)); 

  sel->chooseCrossCorrelations();
    
  for (IConstDataSharedIter it=itsDataSource.createConstIterator(sel,conv);it!=it.end();++it) {  
       const casacore::Vector<casacore::RigidVector<casacore::Double, 3> > testUVWs = simulateUVW(*it);
       const casacore::Vector<casacore::RigidVector<casacore::Double, 3> > measUVWs = it->uvw();
       ASKAPDEBUGASSERT(testUVWs.nelements() == measUVWs.nelements());

       const casacore::MEpoch epoch(casacore::Quantity(it->time()/86400. + itsRefMJD,"d"), casacore::MEpoch::Ref(casacore::MEpoch::UTC));

       bool epochPrinted = false;
       accumulator_set<double, features<tag::min, tag::max> > angleStats;
       accumulator_set<double, features<tag::min, tag::max> > stretchStats;

       const casacore::Cube<casacore::Bool>& flags = it->flag();
       bool hasUnflaggedData = false;
       for (casacore::uInt row = 0; row < it->nRow(); ++row) {
            if (!casacore::allTrue(flags.xyPlane(row))) {
                // this row has unflagged data, UVWs should be good
                hasUnflaggedData = true;
                const casacore::RigidVector<casacore::Double, 3> testUVW = testUVWs[row];
                const casacore::RigidVector<casacore::Double, 3> measUVW = measUVWs[row];
                const double simBslnLength = casacore::sqrt(testUVW * testUVW);
                if (simBslnLength < 1e-6) {
                    ASKAPLOG_WARN_STR(logger, "Encountered zero baseline length for a cross-correlation with unflagged data, antenna ids: "<<
                                      it->antenna1()[row]<<" "<<it->antenna2()[row]<<" beam: "<<it->feed1()[row]<<" epoch: "<<epoch);
                    continue;
                }
                const double measBslnLength = casacore::sqrt(measUVW * measUVW);
                const double cosAngle = simBslnLength > 0. && measBslnLength > 0. ? testUVW * measUVW / (simBslnLength * measBslnLength): 0.;
                ASKAPASSERT(cosAngle <= 1. && cosAngle >= -1.);
                const casacore::RigidVector<casacore::Double, 3> diffUVW = measUVW - testUVW;
                const double diffLength = casacore::sqrt(diffUVW * diffUVW);
                // relative stretch minus 1. (i.e. 0 is the perfect match)
                const double stretch = simBslnLength > 0. ? measBslnLength / simBslnLength - 1. : 0.;
                // angle between two UVW vectors in degrees
                const double angle = casacore::acos(cosAngle) / casacore::C::pi * 180.;
                if (itsVerbose) {
                    if (!epochPrinted) {
                        ASKAPLOG_INFO_STR(logger, "UVW comparison for "<<epoch);
                        epochPrinted = true;
                    }
                    ASKAPLOG_INFO_STR(logger, " "<<it->antenna1()[row]<<" "<<it->antenna2()[row]<<" "<<it->feed1()[row]<<" "<<measUVW<<" "<<testUVW<<" "<<diffUVW<<" "<<diffLength<<" "<<stretch<<" "<<angle);
                }
                angleStats(angle);
                stretchStats(stretch);
            }
       }
       if (hasUnflaggedData) {
           ASKAPLOG_INFO_STR(logger, "For "<<epoch<<" UVW min/max stretch values are "<<min(stretchStats)<<" "<<max(stretchStats)<<
                " min/max angles (deg) are "<<min(angleStats)<<" "<<max(angleStats));
       }
  }
}

/// @brief method predicting UVWs for the given accessor
/// @details For each row a UVW vector is calculated
/// @param[in] acc const accessor to work with (only metadata are used)
/// @return a vector with simulated UVWs, one item per row
casacore::Vector<casacore::RigidVector<casacore::Double, 3> > UVWChecker::simulateUVW(const IConstDataAccessor &acc) const
{
   ASKAPDEBUGASSERT(itsLayout.size() > 0);
   // first, find the number of beams and their phase centres, it looks like we have to handle potential sparse nature of the measurement set 
   // we could've cached the map/phaseCentres because they're changing not that often - possible optimisation for the future
   std::vector<casacore::MVDirection> phaseCentres;
   phaseCentres.reserve(acc.nRow());
   // map to covert real beam indices used in the accessor to indices in the phaseCentres vector or UVW arrays
   std::map<casacore::uInt, size_t> beamIndices;
   for (casacore::uInt row = 0; row < acc.nRow(); ++row) {
        const casa::uInt beam = acc.feed1()[row];
        ASKAPCHECK(beam == acc.feed2()[row], "Cross-beam products are not supported!");
        const casacore::MVDirection phc = acc.pointingDir1()[row];
        ASKAPCHECK(phc.separation(acc.pointingDir2()[row]) < 2e-5, "Phase centres are different for antenna 1 and 2 of the baseline - this is not supported: "<<
                   phc.separation(acc.pointingDir2()[row])<<" "<<printDirection(phc)<<" "<<printDirection(acc.pointingDir2()[row]));
        const std::map<casacore::uInt, size_t>::const_iterator ci = beamIndices.find(beam);
        if (ci != beamIndices.end()) {
            // check that the phase centre is the same as before
            ASKAPDEBUGASSERT(ci->second < phaseCentres.size());
            ASKAPCHECK(phaseCentres[ci->second].separation(phc) < 2e-5, "Phase centres for beam "<<beam + 1<<" (1-based) have changed within one accessor - this is not supported: "<<
                       phaseCentres[ci->second].separation(phc)<<" "<<printDirection(phc)<<" "<<printDirection(phaseCentres[ci->second]));
        } else {
            // this is the new beam
            const size_t newIndex = phaseCentres.size();
            phaseCentres.push_back(phc);
            ASKAPASSERT(phaseCentres.back().separation(phc) < 3e-6);
            beamIndices[beam] = newIndex;
        }
   }

   // now we're ready to compute per-antenna per-beam UVWs
   const casacore::MEpoch epoch(casacore::Quantity(acc.time()/86400. + itsRefMJD,"d"), casacore::MEpoch::Ref(casacore::MEpoch::UTC));

   const casacore::uInt nBeams = static_cast<casacore::uInt>(phaseCentres.size());
   // geocentric U and V per antenna/beam
   casa::Matrix<double> antUs(itsLayout.size(), nBeams, 0.);
   casa::Matrix<double> antVs(itsLayout.size(), nBeams, 0.);
   casa::Matrix<double> antWs(itsLayout.size(), nBeams, 0.);
   std::vector<boost::shared_ptr<casacore::UVWMachine> > uvwMachines(nBeams);
   for (size_t ant = 0; ant < itsLayout.size(); ++ant) {
        const casacore::MPosition antPos = itsLayout[ant];
        const casa::MeasFrame frame(epoch, antPos);
        // antenna position in metres
        const casacore::Vector<casacore::Double> xyz = antPos.getValue().getValue();
        for (casa::uInt beam =0; beam < nBeams; ++beam) {
             // Current APP phase center
             const casa::MDirection fpc = casa::MDirection::Convert(phaseCentres[beam],
                                   casa::MDirection::Ref(casa::MDirection::TOPO, frame))();
             const casa::MDirection hadec = casa::MDirection::Convert(phaseCentres[beam],
                                   casa::MDirection::Ref(casa::MDirection::HADEC, frame))();
             if (ant == 0) {
                 // for uvw rotation
                 // HADEC frame doesn't seem to work correctly with UVW machine, even apart from inversion of the first coordinate
                 // However, it is required for phasing model/UVW itself 
                 // For details see ADESCOM-342.
                 uvwMachines[beam].reset(new casa::UVWMachine(casa::MDirection::Ref(casa::MDirection::J2000), fpc, frame));
                 /*
                 if (beam == 0) {
                     std::cout<<ant<<" "<<epoch<<" "<<antPos<<" "<<printDirection(fpc.getValue())<<" "<<printDirection(hadec.getValue())<<std::endl;
                 }
                 */
             }
             const double dec = hadec.getValue().getLat(); //fpc.getAngle().getValue()(1);
             // hour angle at latitude zero
             const double H0 = hadec.getValue().getLong() - antPos.getValue().getLong();

             // Transformation from antenna position to the geocentric delay
             const double sH0 = sin(H0);
             const double cH0 = cos(H0);
             const double cd = cos(dec);
             const double sd = sin(dec);
             antUs(ant, beam) = -sH0 * xyz(0) - cH0 * xyz(1);
             antVs(ant, beam) = sd * cH0 * xyz(0) - sd * sH0 * xyz(1) - cd * xyz(2);
             antWs(ant, beam) = -cd * cH0 * xyz(0) + cd * sH0 * xyz(1) - sd * xyz(2);
        }
   }
 
   // now everything is ready to compute the result for each baseline
   // note, in principle, the code above can be moved to a separate method or methods
   // something can be cached too, although when time changes recalculation will still need to happen, so probably
   // not a huge overhead to do it the current way
   casacore::Vector<casacore::RigidVector<casacore::Double, 3> > result(acc.nRow());
   casacore::Vector<double> uvwvec(3);
   for (casa::uInt row = 0; row < acc.nRow(); ++row) {
        const casa::uInt ant1 = acc.antenna1()[row];
        const casa::uInt ant2 = acc.antenna2()[row];
        const casa::uInt beamIndex = acc.feed1()[row];
        std::map<casacore::uInt, size_t>::const_iterator ci = beamIndices.find(beamIndex);
        ASKAPDEBUGASSERT(ci != beamIndices.end());
        const casa::uInt beam = static_cast<casa::uInt>(ci->second);
        
        ASKAPASSERT(ant1 < itsLayout.size());
        ASKAPASSERT(ant2 < itsLayout.size());
        ASKAPDEBUGASSERT(beam < nBeams);

        uvwvec(0) = antUs(ant2, beam) - antUs(ant1, beam);
        uvwvec(1) = antVs(ant2, beam) - antVs(ant1, beam);
        uvwvec(2) = antWs(ant2, beam) - antWs(ant1, beam);

        ASKAPDEBUGASSERT(beam < uvwMachines.size());
        ASKAPDEBUGASSERT(uvwMachines[beam]);
        // uvw rotation into J2000
        uvwMachines[beam]->convertUVW(uvwvec);
        ASKAPDEBUGASSERT(uvwvec.nelements() == 3);
        result[row] = uvwvec;
 
   }
   return result;
}

// don't use the whole application harness for now, we don't need any parallelism or passing a parset
// like it is for the standard application
int main(int argc, char **argv) {
  try {
     // borrowed this code from the parallel application class. Perhaps, needs to be abstracted to avoid
     // duplication
     if (!ASKAPLOG_ISCONFIGURED) {
        // Now we have to initialize the logger before we use it
        // If a log configuration exists in the current directory then
        // use it, otherwise try to use the programs default one
        const std::ifstream config("askap.log_cfg", std::ifstream::in);

        if (config) {
            ASKAPLOG_INIT("askap.log_cfg");
        } else {
            std::ostringstream ss;
            ss << argv[0] << ".log_cfg";
            ASKAPLOG_INIT(ss.str().c_str());
            std::cerr<<"initialised for "<<ss.str().c_str()<<std::endl;
        }
     }
     if (argc != 2 && argc !=3) {
         std::cerr<<"Usage "<<argv[0]<<" [-v] measurement_set"<<std::endl;
	 return -2;
     }
     const bool verbose = (argc == 3);
     if (verbose) {
         ASKAPCHECK(std::string("-v") == argv[1], "Incorrect first parameter, only -v is allowed");
     }

     casacore::Timer timer;

     timer.mark();
     TableDataSource ds(verbose ? argv[2] : argv[1],TableDataSource::MEMORY_BUFFERS); 
     UVWChecker checker(ds, verbose);
     ASKAPLOG_DEBUG_STR(logger, "Initialization: "<<timer.real());
     timer.mark();
     checker.run();
     ASKAPLOG_DEBUG_STR(logger,"Job: "<<timer.real());
  }
  catch(const AskapError &ce) {
     std::cerr<<"AskapError has been caught. "<<ce.what()<<std::endl;
     return -1;
  }
  catch(const std::exception &ex) {
     std::cerr<<"std::exception has been caught. "<<ex.what()<<std::endl;
     return -1;
  }
  catch(...) {
     std::cerr<<"An unexpected exception has been caught"<<std::endl;
     return -1;
  }
  return 0;
}
