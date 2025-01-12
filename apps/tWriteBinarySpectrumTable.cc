//
// @file tImageWrite.cc : functional test to exercise image accessor and
//                        write image cube with fake data
//
/// @copyright (c) 2020 CSIRO
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

// std includes
#include <string>
#include <vector>


// ASKAPSoft includes
#include "askap_accessors.h"
#include "askap/AskapLogging.h"
#include "askap/AskapError.h"
#include "askap/Application.h"
#include "askap/StatReporter.h"
#include "askapparallel/AskapParallel.h"
#include "askap/imageaccess/ImageAccessFactory.h"
#include "askap/imageaccess/FitsAuxImageSpectra.h"
#include <askap/scimath/utils/ComplexGaussianNoise.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/IPosition.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <askap/scimath/utils/MultiDimPosIter.h>
#include <casacore/coordinates/Coordinates/DirectionCoordinate.h>
#include <casacore/coordinates/Coordinates/SpectralCoordinate.h>
#include <casacore/coordinates/Coordinates/StokesCoordinate.h>
#include <casacore/coordinates/Coordinates/Projection.h>
#include <askap/scimath/utils/PolConverter.h>
#include <casacore/coordinates/Coordinates/CoordinateSystem.h>
#include <casacore/casa/Containers/Record.h>
#include <ctime>
#include <algorithm>


ASKAP_LOGGER(logger, ".tImageWriteBinaryTable");



// casa
#include <casacore/casa/OS/Timer.h>


// boost
#include <boost/shared_ptr.hpp>

/// 3rd party
#include <Common/ParameterSet.h>

#include <memory>

namespace askap {

namespace accessors {

class TestWriteSpectrumTableApp : public askap::Application {
public:

    void setup()
    {
        std::cout << "Testing" << std::endl;
        LOFAR::ParameterSet parset;
        //parset.add("imagetype","fits");
        itsCurrentRow = 0;
        itsCol = 288;
        casacore::Record record;
        record.define("Stoke","I");
        remove("spectrum_table.fits");
        casacore::CoordinateSystem coord;
        itsFitsAuxImageSpectraTable.reset(new FitsAuxImageSpectra("spectrum_table",itsCol,0,coord,record));
        srand((unsigned int)time(NULL));
    }

   // random number beween id and id+1
   float generate(unsigned int id)
   {
        // generate a random number between 0 and 1 (inclusively)
        float rnum = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        return rnum + id;
   }

   void addNRow(unsigned int nrows)
   {
      std::vector<float> randomSpectrum;

      for (unsigned int row = 0; row < nrows; row++) {
        for(int n = 0; n < itsCol; n++) {
            float r = generate(row);
            randomSpectrum.push_back(r);
        }
        itsCurrentRow += 1;
        //std::string id = std::string("Source_")  + std::to_string(itsCurrentRow);
        std::string id = "SB69705_Sun_AnnularObs_11.0deg_0_238.0_Beam35_component_117a";
        itsFitsAuxImageSpectraTable->add(id,casacore::Vector<float>(randomSpectrum));
        randomSpectrum.resize(0);
     }
   }

   void addNRow2(unsigned int nrows)
   {
      //std::vector<std::vector<float>> arrayOfRandomSpectrum;
      casacore::Matrix<float> arrayOfRandomSpectrum(nrows,itsCol);
      std::vector<std::string> ids;

      //for (unsigned int row = 0; row < nrows; row++) {
      unsigned int num = itsCurrentRow;
      unsigned int j = 0;
      for(int n = 0; n < itsCol; n++) {
        for (unsigned int row = 0; row < nrows; row++) {
            float r = generate(num);
            arrayOfRandomSpectrum(casacore::IPosition(2,row,n)) = r;
            if ( j == itsCol - 1) {
                num += 1;
                j = 0;
            } else {
                j += 1;
            }
        }
     }
     for (unsigned int row = 0; row < nrows; row++) {
        unsigned int r = row + itsCurrentRow;
        //std::string id = std::string("Source_")  + std::to_string(r);
        std::string id = std::string("Source_")  + "SB69705_Sun_AnnularObs_11.0deg_0_238.0_Beam35_component_117a";
        id.append(std::to_string(r));
        ids.push_back(id);
     }
     itsFitsAuxImageSpectraTable->add(ids,arrayOfRandomSpectrum);
     itsCurrentRow += nrows;
   }

   void readSpectrum(long row,casacore::Vector<float>& spectrum)
   {
        itsFitsAuxImageSpectraTable->get(row,spectrum);
        //std::string id = std::string("Source_") + std::to_string(row);
        //itsFitsAuxImageSpectraTable->get(id,spectrum);
   }

   int run(int argc, char* argv[]) final
   {
     try {
        askap::StatReporter stats;
        setup();
        addNRow2(10);
        //addNRow2(10);
        casacore::Vector<float> spectrum;
        readSpectrum(3,spectrum);

        std::vector<float> stdVect  =spectrum.tovector();
        // Because of the way we insert the spectrum to the table, we know
        // the spectrum for row 3 is between 2.0 and 3.0
        bool status = std::all_of(stdVect.begin(),stdVect.end(),
                                  [](float v) { return (v >= 2.0 && v <= 3);});
        //ASKAPCHECK(status, "Error: spectrum in row 3 is not between 2 and 3");
        std::cout << std::endl << "[ ";
        for_each(stdVect.begin(),stdVect.end(),
                    [](float v) {std::cout << v << " ";});
        std::cout << "]" << std::endl;

        stats.logSummary();
        return 0;
     }
     catch (const askap::AskapError& e) {
        ASKAPLOG_FATAL_STR(logger, "Askap error in " << argv[0] << ": " << e.what());
        std::cerr << "Askap error in " << argv[0] << ": " << e.what() << std::endl;
        return 1;
     } catch (const std::exception& e) {
        ASKAPLOG_FATAL_STR(logger, "Unexpected exception in " << argv[0] << ": " << e.what());
        std::cerr << "Unexpected exception in " << argv[0] << ": " << e.what()
                  << std::endl;
        return 1;
     }
   }
private:

   std::string getVersion() const override {
      const std::string pkgVersion = std::string("base-accessor:") + ASKAP_PACKAGE_VERSION;
      return pkgVersion;
   }
   std::unique_ptr<FitsAuxImageSpectra> itsFitsAuxImageSpectraTable;
   unsigned long itsCurrentRow;
   unsigned long itsCol;
   /// @brief image accessor
   //boost::shared_ptr<accessors::IImageAccess<casacore::Float> > itsImageAccessor;

   /// @brief name of the image cube to write
   //std::string itsName;
};


} // namespace accessors

} // namespace askap

int main(int argc, char *argv[])
{
    askap::accessors::TestWriteSpectrumTableApp app;
    return app.main(argc, argv);
}
