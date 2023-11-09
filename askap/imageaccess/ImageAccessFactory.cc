/// @file
/// @brief Build an appropriate image access class
/// @details This file contains a factory method generating a shared pointer to the image
/// accessor from the parset file
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
///

#include <askap/imageaccess/ImageAccessFactory.h>
#include <askap/imageaccess/CasaImageAccess.h>
#include <askap/imageaccess/CasaADIOSImageAccess.h>
#include <askap/imageaccess/CasaADIOSImageAccessParallel.h>
#include <askap/imageaccess/FitsImageAccess.h>
#include <askap/imageaccess/FitsImageAccessParallel.h>

#include <askap/askap/AskapError.h>

#include <string>

using namespace askap;
using namespace askap::accessors;

/// @brief Build an appropriate image access class
/// @details This is a factory method generating a shared pointer to the image
/// accessor from the parset file
/// @param[in] parset parameters containing description of image accessor to be constructed
/// @return shared pointer to the image access object
/// @note CASA images are used by default
boost::shared_ptr<IImageAccess<casacore::Float> > askap::accessors::imageAccessFactory(const LOFAR::ParameterSet &parset)
{
   const std::string imageType = parset.getString("imagetype","casa");

   boost::shared_ptr<IImageAccess<> > result;
   if (imageType == "casa") {
       boost::shared_ptr<CasaImageAccess<casacore::Float> > iaCASA(new CasaImageAccess<casacore::Float>());
       // optional parameter setting may come here
       result = iaCASA;
   } else if (imageType == "adios"){
       std::string config = parset.getString("adiosconfig", "");
       boost::shared_ptr<CasaADIOSImageAccess<casacore::Float> > iaADIOS(new CasaADIOSImageAccess<casacore::Float>(config));
       result = iaADIOS; 
   } else if (imageType == "fits"){
       boost::shared_ptr<FitsImageAccess> iaFITS(new FitsImageAccess());
       const bool fast = (parset.getString("imagealloc","fast") == "fast");
       iaFITS->useFastAlloc(fast);
       result = iaFITS;
  } else {
      throw AskapError(std::string("Unsupported image type ")+imageType+" has been requested");
   }
   return result;

}/// @brief Build an appropriate image access class
/// @details This is a factory method generating a shared pointer to the image
/// accessor from the parset file
/// @param[in] parset parameters containing description of image accessor to be constructed
/// @param[in] comms, MPI communicator
/// @return shared pointer to the image access object
/// @note CASA images are used by default
boost::shared_ptr<IImageAccess<casacore::Float> > askap::accessors::imageAccessFactory(
    const LOFAR::ParameterSet &parset, askapparallel::AskapParallel &comms, size_t comm_index)
{
   const std::string imageType = parset.getString("imagetype","casa");
   const std::string imageAccessType = parset.getString("imageaccess","individual");
   ASKAPCHECK(imageType!="casa" || imageAccessType=="individual","Collective I/O not supported for imagetype casa")

   boost::shared_ptr<IImageAccess<> > result;
   if (imageType == "casa") {
       boost::shared_ptr<CasaImageAccess<casacore::Float> > iaCASA(new CasaImageAccess<casacore::Float>());
       // optional parameter setting may come here
       result = iaCASA;
   } else if (imageType == "adios"){
       std::string config = parset.getString("adiosconfig", "");
       boost::shared_ptr<CasaADIOSImageAccessParallel<casacore::Float> > iaADIOS(new CasaADIOSImageAccessParallel<casacore::Float>(comms, config, comm_index));
       result = iaADIOS; 
   } else if (imageType == "fits"){
       const bool fast = (parset.getString("imagealloc","fast") == "fast");
       if (imageAccessType == "collective") {
           uint axis = parset.getUint("imageaccess.axis",0);
           boost::shared_ptr<FitsImageAccessParallel> iaFITS(new FitsImageAccessParallel(comms,axis));
           iaFITS->useFastAlloc(fast);
           result = iaFITS;
       } else {
           boost::shared_ptr<FitsImageAccess> iaFITS(new FitsImageAccess());
           iaFITS->useFastAlloc(fast);
           result = iaFITS;
       }
   }
   else {
      throw AskapError(std::string("Unsupported image type ")+imageType+" has been requested");
   }
   return result;
}
