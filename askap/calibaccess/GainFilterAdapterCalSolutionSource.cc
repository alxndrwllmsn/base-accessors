/// @file GainFilterAdapterCalSolutionSource.cc
/// @brief implementation of the calibration solution source with gain filtering
/// @details This implementation can apply filtering (e.g., flux normalisation, smoothing)
/// to the gain solutions before they are written out.
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


// own includes
#include <askap/calibaccess/GainFilterAdapterCalSolutionSource.h>
#include <askap/calibaccess/TableCalSolutionSource.h>
#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(logger, "calibaccess.GainFilter");


using namespace std;
using namespace casacore;

namespace askap {

namespace accessors {

/// @brief constructor from cal solution source and parset
/// @details a calibration solution source that can filter the gain solutions
/// @param[in] src a cal solution source
/// @param[in] parset a parameterset with filter parameters
GainFilterAdapterCalSolutionSource::GainFilterAdapterCalSolutionSource(const boost::shared_ptr<TableCalSolutionSource> &src, const LOFAR::ParameterSet &parset):
itsCalSolutionSource(src),itsFluxNorm(parset.getString("fluxnorm","none")),itsAllBeams(parset.getBool("fluxnorm.allbeams",false))
{}


/// @brief destructor
/// @details We need it to modify and flush the solutions at the end
GainFilterAdapterCalSolutionSource::~GainFilterAdapterCalSolutionSource() {
    filter();
}

// const methods of interface

/// @brief obtain ID for the most recent solution
/// @return ID for the most recent solution
long GainFilterAdapterCalSolutionSource::mostRecentSolution() const
{
    return itsCalSolutionSource->mostRecentSolution();
}

/// @brief obtain solution ID for a given time
/// @details This method looks for a solution valid at the given time
/// and returns its ID. It is equivalent to mostRecentSolution() if
/// called with a time sufficiently into the future.
/// @param[in] time time stamp in seconds since MJD of 0.
/// @return solution ID
long GainFilterAdapterCalSolutionSource::solutionID(const double time) const
{
    return itsCalSolutionSource->solutionID(time);
}

/// @brief obtain solution ID for a given time
/// @details This method looks for a solution valid at the given time
/// and returns its ID. This is another name for solutionID indicating
/// it is the last solution before the given time.
/// @param[in] time time stamp in seconds since MJD of 0.
/// @return solution ID, time of solution
std::pair<long, double> GainFilterAdapterCalSolutionSource::solutionIDBefore(const double time) const
{
    return itsCalSolutionSource->solutionIDBefore(time);
}

/// @brief obtain closest solution ID after a given time
/// @details This method looks for the first solution valid after
/// the given time and returns its ID.
/// @param[in] time time stamp in seconds since MJD of 0.
/// @return solution ID, time of solution
std::pair<long, double>  GainFilterAdapterCalSolutionSource::solutionIDAfter(const double time) const
{
    return itsCalSolutionSource->solutionIDAfter(time);
}

/// @brief obtain read-only accessor for a given solution ID
/// @details This method returns a shared pointer to the solution accessor, which
/// can be used to read the parameters. If a solution with the given ID doesn't
/// exist, an exception is thrown. Existing solutions with undefined parameters
/// are managed via validity flags of gains, leakages and bandpasses
/// @param[in] id solution ID to read
/// @return shared pointer to an accessor object
boost::shared_ptr<ICalSolutionConstAccessor> GainFilterAdapterCalSolutionSource::roSolution(const long id) const
{
    ASKAPASSERT(id < itsSolAcc.size());
    return itsSolAcc[id];
}

// read/write virtual methods of the interface

/// @brief obtain a solution ID to store new solution
/// @details This method provides a solution ID for a new solution. It must
/// be called before any write operation (one needs a writable accessor to
/// write the actual solution and to get this accessor one needs an ID).
/// @param[in] time time stamp of the new solution in seconds since MJD of 0.
/// @return solution ID
long GainFilterAdapterCalSolutionSource::newSolutionID(const double time)
{
    const long id = itsCalSolutionSource->newSolutionID(time);
    itsSolAcc.push_back(itsCalSolutionSource->rwSolution(id));
    ASKAPASSERT(id == itsSolAcc.size()-1);
    return id;
}

/// @brief obtain a writeable accessor for a given solution ID
/// @details This method returns a shared pointer to the solution accessor, which
/// can be used to both read the parameters and write them back. If a solution with
/// the given ID doesn't exist, an exception is thrown. Existing solutions with undefined
/// parameters are managed via validity flags of gains, leakages and bandpasses
/// @param[in] id solution ID to access
/// @return shared pointer to an accessor object
boost::shared_ptr<ICalSolutionAccessor> GainFilterAdapterCalSolutionSource::rwSolution(const long id) const
{
    ASKAPASSERT(id < itsSolAcc.size());
    ASKAPASSERT(itsSolAcc[id]);
    return itsSolAcc[id];
}

/// @brief function to do the gain filtering
/// @details filters that needs a timeseries of gains can be added here
void GainFilterAdapterCalSolutionSource::filter()
{
    normaliseFlux();
    // smoothGains();
    // smoothGainAmplitudes();
}

/// @brief function to do flux normalisation of selfcal gains
/// @details when selfcalibrating on still incomplete models, the flux scale
/// can drift - this function takes the mean or median of the squared gain amplitudes
/// and uses that to normalise the gains to avoid changing the fluxscale.
/// The gains will be modified in place, i.e., existing files will be replaced.
/// The X and Y gains are normalised separately to avoid introducing leakage due to e.g.,
/// a single strong source with large (instrumental) polarisation.
void GainFilterAdapterCalSolutionSource::normaliseFlux()
{
    // the nAnt and nBeam parameters will take default values for non table based calibration
    const casacore::uInt nBeam =  itsAllBeams ? itsCalSolutionSource->nBeam() : 1;
    const casacore::uInt nAnt = itsCalSolutionSource->nAnt();
    const long nSol = itsSolAcc.size();

    if (itsFluxNorm == "mean" || itsFluxNorm == "median") {
        ASKAPLOG_INFO_STR(logger, "Rescaling gains to keep the "<<itsFluxNorm<<" flux scale constant");
        // read the valid gains for all timeslots, square them and take mean or median
        for (casacore::uInt beam = 0; beam < nBeam; beam++) {
            vector<float> xGains;
            vector<float> yGains;
            float xMean = 0;
            float yMean = 0;
            for (long id = 0; id < nSol; id++) {
                boost::shared_ptr<ICalSolutionConstAccessor> sol = roSolution(id);
                for (casacore::uInt ant = 0; ant < nAnt; ant++) {
                    JonesIndex i = JonesIndex(ant,beam);
                    JonesJTerm jj = sol->gain(i);
                    if (jj.g1IsValid()) xGains.push_back(std::norm(jj.g1()));
                    if (jj.g2IsValid()) yGains.push_back(std::norm(jj.g2()));
                }
            }
            if (xGains.size()>0 && yGains.size()>0) {
                // scale gains with sqrt of mean or median to keep flux scale the same
                if (itsFluxNorm == "mean") {
                    xMean = sqrt(mean((Vector<float>(xGains))));
                    yMean = sqrt(mean((Vector<float>(yGains))));
                } else {
                    xMean = sqrt(median((Vector<float>(xGains))));
                    yMean = sqrt(median((Vector<float>(yGains))));
                }
            } else {
                ASKAPLOG_WARN_STR(logger,"No valid gains found for X and/or Y pol in normaliseFlux for beam = "<<beam);
            }

            // update gain table
            for (long id = 0; id < nSol; id++) {
                boost::shared_ptr<ICalSolutionAccessor> sol = rwSolution(id);
                for (casacore::uInt ant = 0; ant < nAnt; ant++) {
                    const JonesIndex i = JonesIndex(ant,beam);
                    const JonesJTerm jj = sol->gain(i);
                    const bool v1 = jj.g1IsValid() && (abs(xMean) > 0);
                    Complex g1 = jj.g1();
                    if (v1) {
                        g1 /= xMean;
                    }
                    const bool v2 = jj.g2IsValid() && (abs(yMean) > 0);
                    Complex g2 = jj.g2();
                    if (v2) {
                        g2 /= yMean;
                    }
                    const JonesJTerm newJJ(g1, v1, g2, v2);
                    sol->setGain(i, newJJ);
                }
            }
        }
    } else {
        ASKAPLOG_INFO_STR(logger,"Flux normalisation skipped because fluxnorm="<<itsFluxNorm);
    }
}




} // namespace accessors

} // namespace askap
