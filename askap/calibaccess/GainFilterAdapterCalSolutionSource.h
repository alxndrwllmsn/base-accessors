/// @file GainFilterAdapterCalSolutionSource.h
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

#ifndef ASKAP_ACCESSORS_GAIN_FILTER_ADAPTER_CAL_SOLUTION_SOURCE_H
#define ASKAP_ACCESSORS_GAIN_FILTER_ADAPTER_CAL_SOLUTION_SOURCE_H

#include <askap/calibaccess/TableCalSolutionSource.h>
#include <Common/ParameterSet.h>

namespace askap {

namespace accessors {

/// @brief implementation of the calibration solution source (sink actually) that can filter the gain solutions
/// @details This implementation can apply filtering (e.g., flux normalisation, smoothing)
/// to the gain solutions before they are written out.
/// @ingroup calibaccess
class GainFilterAdapterCalSolutionSource : virtual public ICalSolutionSource {
public:

  /// @brief constructor from cal solution source and parset
  /// @details a calibration solution source that can filter the gain solutions
  /// @param[in] src a cal solution source
  /// @param[in] parset a parameterset with filter parameters
  explicit GainFilterAdapterCalSolutionSource(const boost::shared_ptr<TableCalSolutionSource> &src, const LOFAR::ParameterSet &parset);

   /// @brief destructor
   /// @details We need it to modify and flush the solutions at the end
   ~GainFilterAdapterCalSolutionSource();

  // const methods of interface

  /// @brief obtain ID for the most recent solution
  /// @return ID for the most recent solution
  virtual long mostRecentSolution() const override;

  /// @brief obtain solution ID for a given time
  /// @details This method looks for a solution valid at the given time
  /// and returns its ID. It is equivalent to mostRecentSolution() if
  /// called with a time sufficiently into the future.
  /// @param[in] time time stamp in seconds since MJD of 0.
  /// @return solution ID
  virtual long solutionID(const double time) const override;

  /// @brief obtain solution ID for a given time
  /// @details This method looks for a solution valid at the given time
  /// and returns its ID. This is another name for solutionID indicating
  /// it is the last solution before the given time.
  /// @param[in] time time stamp in seconds since MJD of 0.
  /// @return solution ID, time of solution
  virtual std::pair<long, double> solutionIDBefore(const double time) const override;

  /// @brief obtain closest solution ID after a given time
  /// @details This method looks for the first solution valid after
  /// the given time and returns its ID.
  /// @param[in] time time stamp in seconds since MJD of 0.
  /// @return solution ID, time of solution
  virtual std::pair<long, double>  solutionIDAfter(const double time) const override;

  /// @brief obtain read-only accessor for a given solution ID
  /// @details This method returns a shared pointer to the solution accessor, which
  /// can be used to read the parameters. If a solution with the given ID doesn't
  /// exist, an exception is thrown. Existing solutions with undefined parameters
  /// are managed via validity flags of gains, leakages and bandpasses
  /// @param[in] id solution ID to read
  /// @return shared pointer to an accessor object
  virtual boost::shared_ptr<ICalSolutionConstAccessor> roSolution(const long id) const override;

  // read/write virtual methods of the interface

  /// @brief obtain a solution ID to store new solution
  /// @details This method provides a solution ID for a new solution. It must
  /// be called before any write operation (one needs a writable accessor to
  /// write the actual solution and to get this accessor one needs an ID).
  /// @param[in] time time stamp of the new solution in seconds since MJD of 0.
  /// @return solution ID
  virtual long newSolutionID(const double time) override;

  /// @brief obtain a writeable accessor for a given solution ID
  /// @details This method returns a shared pointer to the solution accessor, which
  /// can be used to both read the parameters and write them back. If a solution with
  /// the given ID doesn't exist, an exception is thrown. Existing solutions with undefined
  /// parameters are managed via validity flags of gains, leakages and bandpasses
  /// @param[in] id solution ID to access
  /// @return shared pointer to an accessor object
  virtual boost::shared_ptr<ICalSolutionAccessor> rwSolution(const long id) const override;

    /// @brief shared pointer definition
  typedef boost::shared_ptr<GainFilterAdapterCalSolutionSource> ShPtr;

protected:

/// @brief function to do the gain filtering
/// @details filters that needs a timeseries of gains can be added here
void filter();

/// @brief function to do flux normalisation of selfcal gains
/// @details when selfcalibrating on still incomplete models, the flux scale
/// can drift - this function takes the mean or median of the squared gain amplitudes
/// and uses that to normalise the gains to avoid changing the fluxscale.
/// The gains will be modified in place, i.e., existing files will be replaced.
/// The X and Y gains are normalised separately to avoid introducing leakage due to e.g.,
/// a single strong source with large (instrumental) polarisation.
void normaliseFlux();

private:
   const std::string itsFluxNorm;
   const bool itsAllBeams;
   const boost::shared_ptr<TableCalSolutionSource> itsCalSolutionSource;
   std::vector<boost::shared_ptr<ICalSolutionAccessor>> itsSolAcc;

}; // class GainFilterAdapterCalSolutionSource

} // namespace accessors

} // namespace askap

#endif
