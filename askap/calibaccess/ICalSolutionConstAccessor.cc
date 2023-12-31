/// @file
///
/// @brief An interface for accessing calibration solutions for reading.
/// @details This interface is used to access calibration parameters
/// read-only. A writable version of the interface is derived from this
/// class. Various implementations are possible, i.e. parset-based,
/// table-based and working via database ice service.
/// @copyright (c) 2011 CSIRO
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
/// @author Max Voronkov <Maxim.Voronkov@csiro.au>
/// Based on the original version of this interface by Ben Humphreys <ben.humphreys@csiro.au>

#include <askap/calibaccess/ICalSolutionConstAccessor.h>
#include <askap/askap/AskapError.h>

namespace askap {

namespace accessors {

// just to keep compiler happy we define an empty virtual destructor
ICalSolutionConstAccessor::~ICalSolutionConstAccessor()
{
}

/// @brief obtain full 2x2 Jones Matrix taking all effects into account
/// @details This method returns resulting 2x2 matrix taking gain, leakage and
/// bandpass effects (for a given channel) into account. Invalid gains (and bandpass
/// values) are replaced by 1., invalid leakages are replaced by zeros. This method
/// calls gain, bandpass and leakage virtual methods
/// @param[in] index ant/beam index
/// @param[in] chan spectral channel of interest
/// @return 2x2 Jones matrix
/// @note The relation between leakage terms and Jones matrices matches
/// the definition of Hamaker, Bregman & Sault. See their equation
/// (14) for details. Our parameters d12 (corresponding to Stokes:XY) and
/// d21 (corresponding to Stokes::YX) correspond to d_{Ap} and d_{Aq} from
/// Hamaker, Bregman & Sault, respectively. It is assumed that the gain errors
/// are applied after leakages (i.e. R=GD).
casacore::SquareMatrix<casacore::Complex, 2> ICalSolutionConstAccessor::jones(const JonesIndex &index, const casacore::uInt chan) const
{
  return jonesAndValidity(index, chan).first;
}

/// @brief obtain full 2x2 Jones Matrix taking all effects into account
/// @details This version of the method accepts antenna and beam indices explicitly and
/// does extra checks before calling the main method expressed via JonesIndex.
/// @param[in] ant antenna index
/// @param[in] beam beam index
/// @param[in] chan spectral channel of interest
/// @return 2x2 Jones matrix
casacore::SquareMatrix<casacore::Complex, 2> ICalSolutionConstAccessor::jones(const casacore::uInt ant,
                                     const casacore::uInt beam, const casacore::uInt chan) const
{
  ASKAPCHECK(chan < 20736, "Channel number is supposed to be less than 20736");
  return jones(JonesIndex(ant, beam), chan);
}

/// @brief obtain validity flag for the full 2x2 Jones Matrix
/// @details This method combines all validity flags for parameters used to compose Jones
/// matrix and returns true if at least one component is defined and false if all constituents
/// are not valid
/// @param[in] index ant/beam index
/// @param[in] chan spectral channel of interest
/// @return true, if the matrix returned by jones(...) method called with the same parameters is
/// valid, false otherwise
bool ICalSolutionConstAccessor::jonesValid(const JonesIndex &index, const casacore::uInt chan) const
{
  return jonesAndValidity(index, chan).second;
}

/// @brief obtain validity flag for the full 2x2 Jones Matrix
/// @details This version of the method accepts antenna and beam indices explicitly and
/// does extra checks before calling the main method expressed via JonesIndex.
/// @param[in] ant antenna index
/// @param[in] beam beam index
/// @param[in] chan spectral channel of interest
/// @return true, if the matrix returned by jones(...) method called with the same parameters is
/// valid, false otherwise
bool ICalSolutionConstAccessor::jonesValid(const casacore::uInt ant, const casacore::uInt beam, const casacore::uInt chan) const
{
  ASKAPCHECK(chan < 20736, "Channel number is supposed to be less than 20736");
  return jonesValid(JonesIndex(ant, beam), chan);
}

/// @brief obtain validity flag for the full 2x2 Jones Matrix
/// @details This version of the method accepts antenna and beam indices explicitly and
/// does extra checks before calling the main method expressed via JonesIndex.
/// @param[in] ant antenna index
/// @param[in] beam beam index
/// @param[in] chan spectral channel of interest
/// @return true, if the matrix returned by jones(...) method called with the same parameters is
/// valid, false otherwise
bool ICalSolutionConstAccessor::jonesAllValid(const casacore::uInt ant, const casacore::uInt beam, const casacore::uInt chan) const
{
  ASKAPCHECK(chan < 20736, "Channel number is supposed to be less than 20736");
  return jonesAllValid(JonesIndex(ant, beam), chan);
}


/// @brief obtain validity flag for the full 2x2 Jones Matrix
/// @details This method combines all validity flags for parameters used to compose Jones
/// matrix and returns true if all elements are valid and false if at least one constituent
/// is not valid
/// @param[in] index ant/beam index
/// @param[in] chan spectral channel of interest
/// @return true, if the matrix returned by jones(...) method called with the same parameters is
/// valid, false otherwise
bool ICalSolutionConstAccessor::jonesAllValid(const JonesIndex &index, const casacore::uInt chan) const
{
  const JonesJTerm gTerm = gain(index);
  const JonesJTerm bpTerm = bandpass(index,chan);
  const JonesDTerm dTerm = leakage(index);
  const JonesDTerm bpdTerm = bpleakage(index,chan);

  return gTerm.g1IsValid() && gTerm.g2IsValid() && bpTerm.g1IsValid() && bpTerm.g2IsValid() &&
         dTerm.d12IsValid() && dTerm.d21IsValid() && bpdTerm.d12IsValid() && bpdTerm.d21IsValid();
}

std::pair<casa::SquareMatrix<casa::Complex, 2>, bool> ICalSolutionConstAccessor::jonesAndValidity(const JonesIndex &index, const casa::uInt chan) const
{
  const JonesJTerm gTerm = gain(index);
  const JonesJTerm bpTerm = bandpass(index,chan);
  const JonesDTerm dTerm = leakage(index);
  const JonesDTerm bpdTerm = bpleakage(index,chan);
  const bool leakageValid = dTerm.d12IsValid() && dTerm.d21IsValid();
  const bool bpleakageValid = bpdTerm.d12IsValid() && bpdTerm.d21IsValid();
  const bool anyLeakageValid = leakageValid || bpleakageValid;

  // MV: there is a serious issue/contradiction in how we use the system now vs. how it was designed
  // the following change is a hack changing the logic w.r.t. documentation (OR instead of AND),
  // we probably have to think how we approach it in the future when we have just bandpass or just gains, etc
  bool valid = (gTerm.g1IsValid() && gTerm.g2IsValid()) || anyLeakageValid ||
              (bpTerm.g1IsValid() && bpTerm.g2IsValid());

  if (!valid) return std::pair<casa::SquareMatrix<casa::Complex, 2>, bool>(casa::SquareMatrix<casa::Complex, 2>(),valid);

  casa::SquareMatrix<casa::Complex, 2> result(anyLeakageValid ? casa::SquareMatrix<casa::Complex, 2>::General :
      casa::SquareMatrix<casa::Complex, 2>::Diagonal);

  result(0,0) = gTerm.g1IsValid() ? gTerm.g1() : casa::Complex(1.,0.);
  result(1,1) = gTerm.g2IsValid() ? gTerm.g2() : casa::Complex(1.,0.);

  // Note: here we assume that only one of leakage and bpleakage is valid for the cases of interest.
  // The maths will get more complicated otherwise
  if (leakageValid) {
     result(0,1) = (dTerm.d12IsValid() ? dTerm.d12() : 0.) * result(0,0);
     result(1,0) = (dTerm.d21IsValid() ? -dTerm.d21() : 0.) * result(1,1);
  } else if (bpleakageValid) {
     result(0,1) = (bpdTerm.d12IsValid() ? bpdTerm.d12() : 0.) * result(0,0);
     result(1,0) = (bpdTerm.d21IsValid() ? -bpdTerm.d21() : 0.) * result(1,1);
  }

  if (bpTerm.g1IsValid()) {
      result(0,0) *= bpTerm.g1();
      if (anyLeakageValid) result(1,0) *= bpTerm.g1();
  }
  if (bpTerm.g2IsValid()) {
      if (anyLeakageValid) result(0,1) *= bpTerm.g2();
      result(1,1) *= bpTerm.g2();
  }
  return std::pair<casa::SquareMatrix<casa::Complex, 2>, bool>(result,valid);
}

std::pair<casa::SquareMatrix<casa::Complex, 2>, bool> ICalSolutionConstAccessor::jonesAndValidity(const casa::uInt ant, const::casa::uInt beam, const casa::uInt chan) const
{
    return jonesAndValidity(JonesIndex(ant,beam),chan);
}

} // namespace accessors
} // namespace askap
