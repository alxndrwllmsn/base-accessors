/// @file DopplerConverter.h
/// @brief A class for interconversion between frequencies
/// and velocities
/// @details This is an implementation of a relatively low-level
/// interface, which is used within the implementation of the data
/// accessor. The end user interacts with the IDataConverter class only.
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

#ifndef ASKAP_ACCESSORS_DOPPLER_CONVERTER_H
#define ASKAP_ACCESSORS_DOPPLER_CONVERTER_H

// CASA includes
#include <casacore/measures/Measures/MDoppler.h>
#include <casacore/measures/Measures/MCDoppler.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/casa/Quanta/MVFrequency.h>

// own includes
#include <askap/dataaccess/IDopplerConverter.h>
#include <askap/dataaccess/DataAccessError.h>


namespace askap {

namespace accessors {

/// @brief A class for interconversion between frequencies
/// and velocities
/// @details This is an implementation of a relatively low-level
/// interface, which is used within the implementation of the data
/// accessor. The end user interacts with the IDataConverter class only.
///
/// The idea behind this class is very similar to CASA's VelocityMachine,
/// but we require a bit different interface to use the class efficiently
/// (and the interface conversion would be equivalent in complexity to
/// the transformation itself). Hence, we will use this class 
/// instead of the VelocityMachine
/// @ingroup dataaccess_conv
struct DopplerConverter : virtual public IDopplerConverter {

    /// constructor
    /// @param[in] restFreq The rest frequency used for interconversion between
    ///                 frequencies and velocities
    /// @param[in] velType velocity (doppler) type (i.e. radio, optical)
    /// Default is radio definition.
    explicit DopplerConverter(const casacore::MVFrequency &restFreq,
                              casacore::MDoppler::Types velType =
                              casacore::MDoppler::RADIO);
    
    /// convert specified frequency to velocity in the same reference
    /// frame. Velocity definition (i.e. optical or radio, etc) is
    /// determined by the implementation class.
    ///
    /// @param[in] freq an MFrequency measure to convert.
    /// @return a reference on MRadialVelocity object with the result
    virtual const casacore::MRadialVelocity& operator()(const casacore::MFrequency &freq) const;

    /// convert specified velocity to frequency in the same reference
    /// frame. Velocity definition (i.e. optical or radio, etc) is
    /// determined by the implementation class.
    ///
    /// @param[in] vel an MRadialVelocity measure to convert.
    /// @return a reference on MFrequency object with the result
    virtual const casacore::MFrequency& operator()(const casacore::MRadialVelocity &vel) const;
protected:
    /// setting the measure frame doesn't make sense for this class
    /// because we're not doing conversions here. This method is empty.
    /// Defined here to make the compiler happy.
    ///
    /// @param[in] frame  MeasFrame object (can be constructed from
    ///               MPosition or MEpoch on-the-fly). Not used.
    virtual void setMeasFrame(const casacore::MeasFrame &frame);

    /// convert frequency frame type to velocity frame type
    /// @param[in] type frequency frame type to convert
    /// @return resulting velocity frame type
    ///
    /// Note, an exception is thrown if the the frame type is
    /// MFrequency::REST (it doesn't make sense to always return zero
    /// velocity).
    static casacore::MRadialVelocity::Types
        freqToVelType(casacore::MFrequency::Types type);

    /// convert velocity frame type to frequency frame type
    /// @param[in] type velocity frame type to convert
    /// @return resulting frequency frame type
    static casacore::MFrequency::Types
      velToFreqType(casacore::MRadialVelocity::Types type);
		 
private:
    /// doppler converters:
    /// from own velocity type specified in the constructor
    ///   to BETA (true velocity)
    mutable casacore::MDoppler::Convert itsToBettaConv;
    /// from true velocity to velocity type specified in the constructor
    mutable casacore::MDoppler::Convert itsFromBettaConv;

    /// rest frequency required for conversion in Hz
    casacore::Double itsRestFrequency;

    /// result buffers
    mutable casacore::MRadialVelocity itsRadialVelocity;
    mutable casacore::MFrequency itsFrequency;
};

} // namespace accessors

} // namespace askap



#endif // #ifndef DOPPLER_CONVERTER_H
