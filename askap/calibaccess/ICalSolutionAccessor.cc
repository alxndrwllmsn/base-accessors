/// @file
/// @brief An interface for accessing calibration solutions for reading and writing.
/// @details This interface is used to access calibration parameters for both
/// reading and writing. It is derived from read-only version of the interface.
/// Various implementations are possible, i.e. parset-based,
/// table-based and working via database ice service.
///
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

#include <askap/calibaccess/ICalSolutionAccessor.h>
#include <askap/askap/AskapError.h>

namespace askap {

namespace accessors {

/// @brief set a single element of the Jones matrix (i.e gains or leakages)
/// @details This method simplifies writing both gains and leakages solution. It reads the current
/// gains and leakages and then replaces one element with the given value setting the validity flag.
/// The stokes parameter is controlling which element of the Jones matrix is replaced. We assume that
/// only linear polarisation products are used (an exception is thrown if it is not the case).
/// XX and YY represent parallel-hand gains (two elements of JonesJTerm) and XY and YX represent
/// cross-pol leakages (two elements of JonesDTerm).
/// @param[in] index ant/beam index
/// @param[in] stokes what element to update (choose from XX,XY,YX and YY)
/// @param[in] elem value to set
void ICalSolutionAccessor::setJonesElement(const JonesIndex &index, const casacore::Stokes::StokesTypes stokes, const casacore::Complex &elem)
{
  if ( (stokes == casacore::Stokes::XX) || (stokes == casacore::Stokes::YY) ) {
      // parallel-hand case
      const JonesJTerm oldJTerm = gain(index);
      setGain(index, stokes == casacore::Stokes::XX ? JonesJTerm(elem, true, oldJTerm.g2(), oldJTerm.g2IsValid()) :
                     JonesJTerm(oldJTerm.g1(), oldJTerm.g1IsValid(), elem, true));
  } else if ( (stokes == casacore::Stokes::XY) || (stokes == casacore::Stokes::YX) ) {
      // cross-pol case (need to implement validity flags at some stage)
      const JonesDTerm oldDTerm = leakage(index);
      setLeakage(index, stokes == casacore::Stokes::XY ? JonesDTerm(elem, true, oldDTerm.d21(), oldDTerm.d21IsValid()) :
                     JonesDTerm(oldDTerm.d12(), oldDTerm.d12IsValid(), elem, true));
  } else {
      ASKAPTHROW(AskapError, "Only XX, YY, XY and YX stokes are supported by setJonesElement, you passed stokes="<<stokes);
  }
}

/// @brief set a single element of the Jones matrix (i.e. gains or leakages)
/// @details This version of the method gets explicitly defined antenna and beam indices.
/// @param[in] ant ant index
/// @param[in] beam beam index
/// @param[in] stokes what element to update (choose from XX,XY,YX and YY)
/// @param[in] elem value to set
void ICalSolutionAccessor::setJonesElement(casacore::uInt ant, casacore::uInt beam, const casacore::Stokes::StokesTypes stokes, const casacore::Complex &elem)
{
  setJonesElement(JonesIndex(ant, beam), stokes, elem);
}

/// @brief set a single element of bandpass
/// @details This method simplifies writing bandpass solution. It reads the current frequency-dependent
/// gains for the given channel and then replaces one of the elements with the given value setting
/// the validity flag. We assume that only linear polarisation frame is to be used with this method
/// (an exception is thrown if it is not the case).
/// @param[in] index ant/beam index
/// @param[in] stokes what element to update (choose either XX or YY)
/// @param[in] chan spectral channel of interest
/// @param[in] elem value to set
void ICalSolutionAccessor::setBandpassElement(const JonesIndex &index, const casacore::Stokes::StokesTypes stokes,
                                              casacore::uInt chan, const casacore::Complex &elem)
{
  if (stokes == casacore::Stokes::XX || stokes == casacore::Stokes::YY) {
      const JonesJTerm oldBP = bandpass(index,chan);
      if (stokes == casacore::Stokes::XX) {
          const JonesJTerm newBP(elem, true, oldBP.g2(), oldBP.g2IsValid());
          setBandpass(index, newBP, chan);
      } else if (stokes == casacore::Stokes::YY) {
          const JonesJTerm newBP(oldBP.g1(), oldBP.g1IsValid(),elem,true);
          setBandpass(index, newBP, chan);
      }
  } else if (stokes == casacore::Stokes::XY || stokes == casacore::Stokes::YX) {
      const JonesDTerm oldBPL = bpleakage(index,chan);
      if (stokes == casacore::Stokes::XY) {
          const JonesDTerm newBPL(elem, true, oldBPL.d21(), oldBPL.d21IsValid());
          setBPLeakage(index, newBPL, chan);
      } else if (stokes == casacore::Stokes::YX) {
          const JonesDTerm newBPL(oldBPL.d12(), oldBPL.d12IsValid(),elem,true);
          setBPLeakage(index, newBPL, chan);
      }
  } else {
      ASKAPTHROW(AskapError, "Only XX, XY, YX and YY stokes are supported by setBandpassElement, you passed stokes="<<stokes);
  }
}


/// @brief set a single element of bandpass
/// @details This version of the method uses explicitly defined antenna and beam indices.
/// @param[in] ant ant index
/// @param[in] beam beam index
/// @param[in] stokes what element to update (choose either XX or YY)
/// @param[in] chan spectral channel of interest
/// @param[in] elem value to set
void ICalSolutionAccessor::setBandpassElement(casacore::uInt ant, casacore::uInt beam,
                                              const casacore::Stokes::StokesTypes stokes, casacore::uInt chan,
                                              const casacore::Complex &elem)
{
  ASKAPCHECK(chan < 20736, "Channel number is supposed to be less than 20736");
  setBandpassElement(JonesIndex(ant, beam), stokes, chan, elem);
}

/// @brief set a single ionospheric parameter
/// @details This method simplifies the writing of ionospheric solutions. It reads the current
/// solution and then replaces one element with the given value setting the validity flag.
/// @param[in] index ant/beam index
/// @param[in] elem value to set
void ICalSolutionAccessor::setIonosphericElement(const JonesIndex &index, const casacore::Complex &elem)
{
  //const IonoTerm oldITerm = ionoparam(index);
  setIonosphere(index, IonoTerm(elem, true));

}
/// @brief set a single ionospheric parameter
/// @details This version of the method uses explicitly defined antenna and beam indices.
/// @param[in] param parameter index
/// @param[in] dir direction index
/// @param[in] elem value to set
void ICalSolutionAccessor::setIonosphericElement(const casacore::uInt param, const casacore::uInt dir,
                                                 const casacore::Complex &elem)
{
  setIonosphericElement(JonesIndex(param, dir), elem);
}

} // namespace accessors

} // namespace askap
