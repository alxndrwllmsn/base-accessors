//
// @file tDataAccess.cc : evolving test/demonstration program of the
//                        data access layer
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
///
/// @author Max Voronkov <maxim.voronkov@csiro.au>


#include <askap/dataaccess/TableDataSource.h>
#include <askap_accessors.h>
#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".tDataAccess");

#include <askap/askap/AskapError.h>
#include <askap/dataaccess/SharedIter.h>
#include <askap/dataaccess/ParsetInterface.h>

#include <askap/dataaccess/TableManager.h>
#include <askap/dataaccess/IDataConverterImpl.h>

// casa
#include <casacore/measures/Measures/MFrequency.h>
#include <casacore/tables/Tables/Table.h>
#include <casacore/casa/OS/Timer.h>
#include <casacore/casa/Arrays/ArrayMath.h>


// std
#include <stdexcept>
#include <iostream>

using std::cout;
using std::cerr;
using std::endl;

using namespace askap;
using namespace accessors;

void timeDependentSubtableTest(const string &ms, const IConstDataSource &ds) 
{
  IDataConverterPtr conv=ds.createConverter();  
  //conv->setEpochFrame(casacore::MEpoch(casacore::Quantity(53635.5,"d"),
  //                    casacore::MEpoch::Ref(casacore::MEpoch::UTC)),"s");
  IDataSelectorPtr sel=ds.createSelector();
  //sel->chooseFeed(1);
  //sel<<LOFAR::ParameterSet("test.in").makeSubset("TestSelection.");
  const IDataConverterImpl &dci=dynamic_cast<const IDataConverterImpl&>(*conv);
  const TableManager tm(casacore::Table(ms),true);
  const IFeedSubtableHandler &fsh = tm.getFeed();  
  for (IConstDataSharedIter it=ds.createConstIterator(sel,conv);it!=it.end();++it) {  
       cout<<"direction: "<<it->pointingDir2()<<endl;
       cout<<"time: "<<it->time()<<" "<<dci.epochMeasure(it->time())<<" "<<
             fsh.getAllBeamOffsets(dci.epochMeasure(it->time()),0)<<endl;
  }
}

void doReadOnlyTest(const IConstDataSource &ds) {
  IDataSelectorPtr sel=ds.createSelector();
  //sel->chooseFeed(1);
  //sel<<LOFAR::ParameterSet("test.in").makeSubset("TestSelection.");
  IDataConverterPtr conv=ds.createConverter();  
  conv->setFrequencyFrame(casacore::MFrequency::Ref(casacore::MFrequency::BARY),"MHz");
  conv->setEpochFrame(casacore::MEpoch(casacore::Quantity(53635.5,"d"),
                      casacore::MEpoch::Ref(casacore::MEpoch::UTC)),"s");
  conv->setDirectionFrame(casacore::MDirection::Ref(casacore::MDirection::AZEL));                    
    
  for (IConstDataSharedIter it=ds.createConstIterator(sel,conv);it!=it.end();++it) {  
       //cout<<"this is a test "<<it->visibility().nrow()<<" "<<it->frequency()<<endl;
       //cout<<"flags: "<<it->flag()<<endl;
       //cout<<"feed1 pa: "<<it->feed1PA()<<endl;
       cout<<"w: [";
       for (casacore::uInt row = 0; row<it->nRow(); ++row) {
            cout<<it->uvw()[row](2);
	    if (row + 1 != it->nRow()) {
	        cout<<", ";
	    }
       }
       cout<<"]"<<endl;
       //cout<<"noise: "<<it->noise().shape()<<endl;
       //cout<<"direction: "<<it->pointingDir2()<<endl;
       //cout<<"direction: "<<it->dishPointing2()<<endl;
       //cout<<"ant1: "<<it->antenna1()<<endl;
       //cout<<"ant2: "<<it->antenna2()<<endl;
       cout<<"time: "<<it->time()<<endl;
  }
}

void doReadWriteTest(const IDataSource &ds) {
  IDataSelectorPtr sel=ds.createSelector();
  //sel->chooseFeed(1);  
  IDataConverterPtr conv=ds.createConverter();
  conv->setFrequencyFrame(casacore::MFrequency::Ref(casacore::MFrequency::TOPO),"MHz");
  conv->setEpochFrame(casacore::MEpoch(casacore::Quantity(53635.5,"d"),
                      casacore::MEpoch::Ref(casacore::MEpoch::UTC)),"s");
  IDataSharedIter it=ds.createIterator(sel,conv);
  //for (size_t run=0;run<10;++run)
  for (it.init();it!=it.end();it.next()) {
       //cout<<it.buffer("TEST").rwVisibility()<<endl;
       it->frequency();
       it->pointingDir1();
       it->time();
       it->antenna1();
       it->feed1();
       it->uvw();
       //it.buffer("TEST").rwVisibility()=it->visibility();
       //it.chooseBuffer("MODEL_DATA");
       //it->rwVisibility()=it.buffer("TEST").visibility();
       it.chooseOriginal();
       it->rwVisibility().set(casacore::Complex(1.,0.0));
       const double l=0., m=0.003975472185;
       for (casacore::uInt row = 0; row<it->nRow(); ++row) {
            for (casacore::uInt chan=0; chan<it->nChannel(); ++chan) {
                 const double phase = 2.*casacore::C::pi*(it->uvw()(row)(0)*l+it->uvw()(row)(1)*m)/casacore::C::c*it->frequency()[chan]*1e6;
                 const casacore::Complex phasor(cos(phase),sin(phase));
                 casacore::Array<casacore::Complex> tmp = it->rwVisibility().yzPlane(row).row(chan);
                 tmp *= phasor;
            }
       }
  }
}

int main(int argc, char **argv) {
  try {
     if (argc!=2) {
         cerr<<"Usage "<<argv[0]<<" measurement_set"<<endl;
	 return -2;
     }

     casacore::Timer timer;

     timer.mark();
     //TableDataSource ds(argv[1],TableDataSource::REMOVE_BUFFERS |
     //                           TableDataSource::MEMORY_BUFFERS);     
     //TableDataSource ds(argv[1],TableDataSource::MEMORY_BUFFERS | TableDataSource::WRITE_PERMITTED);     
     TableDataSource ds(argv[1],TableDataSource::MEMORY_BUFFERS);     
     std::cerr<<"Initialization: "<<timer.real()<<std::endl;
     //timeDependentSubtableTest(argv[1],ds);
     timer.mark();
     doReadOnlyTest(ds);
     //doReadWriteTest(ds);    
     std::cerr<<"Job: "<<timer.real()<<std::endl;
     
  }
  catch(const AskapError &ce) {
     cerr<<"AskapError has been caught. "<<ce.what()<<endl;
     return -1;
  }
  catch(const std::exception &ex) {
     cerr<<"std::exception has been caught. "<<ex.what()<<endl;
     return -1;
  }
  catch(...) {
     cerr<<"An unexpected exception has been caught"<<endl;
     return -1;
  }
  return 0;
}
