/// @file
/// @brief an implementation of IConstDataAccessor
///
/// @details TableConstDataAccessor is an implementation of the
/// DataAccessor working with TableConstDataIterator. It is currently
/// derived from DataAccessorStub as most of the
/// methods are stubbed. However, in the future
/// it should become a separate class derived
/// directly from its interface
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
#ifndef ASKAP_ACCESSORS_TABLE_CONST_DATA_ACCESSOR_H
#define ASKAP_ACCESSORS_TABLE_CONST_DATA_ACCESSOR_H

// own includes
#include <askap/dataaccess/IConstDataAccessor.h>
#include <askap/dataaccess/DataAccessError.h>
#include <askap/dataaccess/CachedAccessorField.h>
#include <askap/dataaccess/UVWRotationHandler.h>

namespace askap {
	
namespace accessors {

/// to be able to link this class to appropriate iterator
/// @ingroup dataaccess_tab
class TableConstDataIterator;


/// @brief an implementation of IConstDataAccessor in the table-based case
///
/// @details TableConstDataAccessor is an implementation of the
/// DataAccessor working with TableConstDataIterator.
/// It is currently
/// derived from DataAccessorStub as most of the
/// methods are stubbed. However, in the future
/// it should become a separate class derived
/// directly from its interface
/// @ingroup dataaccess_tab
class TableConstDataAccessor : virtual public IConstDataAccessor
{
public:
  /// construct an object linked with the given iterator
  /// @param iter a reference to associated iterator
  explicit TableConstDataAccessor(const TableConstDataIterator &iter);

  /// The number of rows in this chunk
  /// @return the number of rows in this chunk
  virtual casacore::uInt nRow() const throw();

  /// The number of spectral channels (equal for all rows)
  /// @return the number of spectral channels
  virtual casacore::uInt nChannel() const throw();

  /// The number of polarization products (equal for all rows)
  /// @return the number of polarization products (can be 1,2 or 4)
  virtual casacore::uInt nPol() const throw();
  
  /// Return pointing centre directions of the first antenna/feed
  /// @return a vector with direction measures (coordinate system
  /// is set via IDataConverter), one direction for each
  /// visibility/row
  virtual const casacore::Vector<casacore::MVDirection>& pointingDir1() const;

  /// Return pointing centre directions of the second antenna/feed
  /// @return a vector with direction measures (coordinate system
  /// is set via IDataConverter), one direction for each
  /// visibility/row
  virtual const casacore::Vector<casacore::MVDirection>& pointingDir2() const;

  /// pointing direction for the centre of the first antenna 
  /// @details The same as pointingDir1, if the feed offsets are zero
  /// @return a vector with direction measures (coordinate system
  /// is is set via IDataConverter), one direction for each
  /// visibility/row
  virtual const casacore::Vector<casacore::MVDirection>& dishPointing1() const;

  /// pointing direction for the centre of the first antenna 
  /// @details The same as pointingDir2, if the feed offsets are zero
  /// @return a vector with direction measures (coordinate system
  /// is is set via IDataConverter), one direction for each
  /// visibility/row
  virtual const casacore::Vector<casacore::MVDirection>& dishPointing2() const;


  /// Visibilities (a cube is nRow x nChannel x nPol; each element is
  /// a complex visibility)
  /// @return a reference to nRow x nChannel x nPol cube, containing
  /// all visibility data
  virtual const casacore::Cube<casacore::Complex>& visibility() const;

  /// Cube of flags corresponding to the output of visibility() 
  /// @return a reference to nRow x nChannel x nPol cube with flag 
  ///         information. If True, the corresponding element is flagged.
  virtual const casacore::Cube<casacore::Bool>& flag() const;

  /// UVW
  /// @return a reference to vector containing uvw-coordinates
  /// packed into a 3-D rigid vector
  virtual const casacore::Vector<casacore::RigidVector<casacore::Double, 3> >&uvw() const;

  /// @brief uvw after rotation
  /// @details This method calls UVWMachine to rotate baseline coordinates 
  /// for a new tangent point. Delays corresponding to this correction are
  /// returned by a separate method.
  /// @param[in] tangentPoint tangent point to rotate the coordinates to
  /// @return uvw after rotation to the new coordinate system for each row
  virtual const casacore::Vector<casacore::RigidVector<casacore::Double, 3> >&
	           rotatedUVW(const casacore::MDirection &tangentPoint) const;
	         
  /// @brief delay associated with uvw rotation
  /// @details This is a companion method to rotatedUVW. It returns delays corresponding
  /// to the baseline coordinate rotation. An additional delay corresponding to the 
  /// translation in the tangent plane can also be applied using the image 
  /// centre parameter. Set it to tangent point to apply no extra translation.
  /// @param[in] tangentPoint tangent point to rotate the coordinates to
  /// @param[in] imageCentre image centre (additional translation is done if imageCentre!=tangentPoint)
  /// @return delays corresponding to the uvw rotation for each row
  virtual const casacore::Vector<casacore::Double>& uvwRotationDelay(
	       const casacore::MDirection &tangentPoint, const casacore::MDirection &imageCentre) const;
  
  /// Frequency for each channel
  /// @return a reference to vector containing frequencies for each
  ///         spectral channel (vector size is nChannel). Frequencies
  ///         are given as Doubles, the frame/units are specified by
  ///         the DataSource object
  virtual const casacore::Vector<casacore::Double>& frequency() const;

  /// Timestamp for each row
  /// @return a timestamp for this buffer (it is always the same
  ///         for all rows. The timestamp is returned as 
  ///         Double w.r.t. the origin specified by the 
  ///         DataSource object and in that reference frame
  virtual casacore::Double time() const;
  
  /// First antenna IDs for all rows
  /// @return a vector with IDs of the first antenna corresponding
  /// to each visibility (one for each row)
  virtual const casacore::Vector<casacore::uInt>& antenna1() const;

  /// Second antenna IDs for all rows
  /// @return a vector with IDs of the second antenna corresponding
  /// to each visibility (one for each row)
  virtual const casacore::Vector<casacore::uInt>& antenna2() const;

  /// First feed IDs for all rows
  /// @return a vector with IDs of the first feed corresponding
  /// to each visibility (one for each row)
  virtual const casacore::Vector<casacore::uInt>& feed1() const;

  /// Second feed IDs for all rows
  /// @return a vector with IDs of the second feed corresponding
  /// to each visibility (one for each row)
  virtual const casacore::Vector<casacore::uInt>& feed2() const;
  
  /// Position angles of the first feed for all rows
  /// @return a vector with position angles (in radians) of the
  /// first feed corresponding to each visibility
  virtual const casacore::Vector<casacore::Float>& feed1PA() const;
     
  /// Position angles of the second feed for all rows
  /// @return a vector with position angles (in radians) of the
  /// second feed corresponding to each visibility
  virtual const casacore::Vector<casacore::Float>& feed2PA() const;
  
  /// Noise level required for a proper weighting
  /// @return a reference to nRow x nChannel x nPol cube with
  ///         complex noise estimates
  virtual const casacore::Cube<casacore::Complex>& noise() const;
  
  /// Velocity for each channel
  /// @return a reference to vector containing velocities for each
  ///         spectral channel (vector size is nChannel). Velocities
  ///         are given as Doubles, the frame/units are specified by
  ///         the DataSource object (via IDataConverter).
  virtual const casacore::Vector<casacore::Double>& velocity() const;

  /// @brief polarisation type for each product
  /// @return a reference to vector containing polarisation types for
  /// each product in the visibility cube (nPol() elements).
  /// @note All rows of the accessor have the same structure of the visibility
  /// cube, i.e. polarisation types returned by this method are valid for all rows.
  virtual const casacore::Vector<casacore::Stokes::StokesTypes>& stokes() const;
  

  /// @brief invalidate fields updated on each iteration
  /// @details Such caches like visibility, uvw, noise and flags are updated
  /// on each new iteration. These are invalidated by call to this methid.
  /// Caches of frequency/velocity axis are updated less regularly (may be
  /// only once if the is just one spectral window in the measurement set).
  /// These are invalidated by a call to notifyNewSpectralWindow(), if
  /// the new window is not the same as the cached one
  virtual void invalidateIterationCaches() const throw();
  
  /// @brief invalidate fields corresponding to the spectral axis
  /// @details See invalidateIterationCaches for more details
  void invalidateSpectralCaches() const throw();
  
  /// @brief invalidate cache of rotated uvw and delays
  /// @details Cache of rotated uvw and delays is kept per accessor, need this
  /// method to access private field
  void invalidateRotatedUVW() const throw();

  /// @brief Obtain a const reference to associated iterator.
  /// @details This method is mainly intended to be used in the derived
  /// non-const implementation, which works with a different type of the
  /// iterator.
  /// @return a const reference to the associated iterator
  const TableConstDataIterator& iterator() const;
private:  
  /// a helper adapter method to set the time via non-const reference
  /// @param[in] time a reference to buffer to fill with the current time 
  void readTime(casacore::Double &time) const;
  
  /// a reference to iterator managing this accessor
  const TableConstDataIterator& itsIterator;
  
  /// internal buffer for visibility
  CachedAccessorField<casacore::Cube<casacore::Complex> > itsVisibility;
  
  /// internal buffer for flag
  CachedAccessorField<casacore::Cube<casacore::Bool> > itsFlag;
 
  /// internal buffer for uvw
  CachedAccessorField<casacore::Vector<casacore::RigidVector<casacore::Double, 3> > > itsUVW;
  
  /// internal buffer for rotated uvw and associated delay
  UVWRotationHandler itsRotatedUVW; 
 
  /// internal buffer for frequency
  CachedAccessorField<casacore::Vector<casacore::Double> > itsFrequency;

  /// internal buffer for time
  CachedAccessorField<casacore::Double> itsTime;
  
  /// internal buffer for the first antenna ids
  CachedAccessorField<casacore::Vector<casacore::uInt> > itsAntenna1;

  /// internal buffer for the second antenna ids
  CachedAccessorField<casacore::Vector<casacore::uInt> > itsAntenna2;
  
  /// internal buffer for the first feed ids
  CachedAccessorField<casacore::Vector<casacore::uInt> > itsFeed1;
  
  /// internal buffer for the second feed ids
  CachedAccessorField<casacore::Vector<casacore::uInt> > itsFeed2;
  
  /// internal buffer for the position angles of the first feed
  CachedAccessorField<casacore::Vector<casacore::Float> > itsFeed1PA;
  
  /// internal buffer for the position angles of the second feed
  CachedAccessorField<casacore::Vector<casacore::Float> > itsFeed2PA;
  
  /// internal buffer for the pointing directions of the first antenna/feed
  CachedAccessorField<casacore::Vector<casacore::MVDirection> > itsPointingDir1;

  /// internal buffer for the pointing directions of the second antenna/feed
  CachedAccessorField<casacore::Vector<casacore::MVDirection> > itsPointingDir2;

  /// internal buffer for the pointing directions of the centre of the first antenna
  CachedAccessorField<casacore::Vector<casacore::MVDirection> > itsDishPointing1;

  /// internal buffer for the pointing directions of the centre of the second antenna
  CachedAccessorField<casacore::Vector<casacore::MVDirection> > itsDishPointing2;
  
  /// internal buffer for the noise figures
  CachedAccessorField<casacore::Cube<casacore::Complex> > itsNoise;
  
  /// internal buffer for the polarisation types
  CachedAccessorField<casacore::Vector<casacore::Stokes::StokesTypes> > itsStokes;
};


} // namespace accessors

} // namespace askap

#endif // #ifndef TABLE_CONST_DATA_ACCESSOR_H
