/// @file TableScalarFieldSelector.cc
/// @brief An implementation of ITableDataSelectorImpl for simple (scalar) fields, like feed ID.
/// @details This class represents a selection of visibility
///         data according to some criterion. This is an
///         implementation of the part of the IDataSelector
///         interface, which can be done with the table selection
///         mechanism in the table based case. Only simple
///         (scalar) fields are included in this selection.
///         Epoch-based selection is done via a separate class
///         because a fully defined converter is required to
///         perform such selection.
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

/// own include
#include <askap/dataaccess/TableScalarFieldSelector.h>
#include <askap/dataaccess/DataAccessError.h>
#include <casacore/tables/TaQL/ExprNodeSet.h>


using namespace askap;
using namespace askap::accessors;
using namespace casa;


/// Choose a single feed, the same for both antennae
/// @param feedID the sequence number of feed to choose
void TableScalarFieldSelector::chooseFeed(casacore::uInt feedID)
{
   if (itsTableSelector.isNull()) {
       itsTableSelector= (table().col("FEED1") ==
                  static_cast<casacore::Int>(feedID)) && (table().col("FEED2") ==
                  static_cast<casacore::Int>(feedID));
   } else {
       itsTableSelector=itsTableSelector && (table().col("FEED1") ==
                  static_cast<casacore::Int>(feedID)) && (table().col("FEED2") ==
                  static_cast<casacore::Int>(feedID));
   }
}

/// @brief choose user-defined index
/// @param[in] column column name in the measurement set for a user-defined index
/// @param[in] value index value
void TableScalarFieldSelector::chooseUserDefinedIndex(const std::string &column, const casacore::uInt value)
{
   if (itsTableSelector.isNull()) {
       itsTableSelector = (table().col(column) == value);
   } else {
       itsTableSelector = itsTableSelector && (table().col(column) == value);
   }
}

/// Choose a single baseline
/// @param ant1 the sequence number of the first antenna
/// @param ant2 the sequence number of the second antenna
/// Which one is the first and which is the second is not important
void TableScalarFieldSelector::chooseBaseline(casacore::uInt ant1,
                                              casacore::uInt ant2)
{
   if (itsTableSelector.isNull()) {
       itsTableSelector= (table().col("ANTENNA1") ==
           static_cast<casacore::Int>(ant1)) && (table().col("ANTENNA2") ==
	   static_cast<casacore::Int>(ant2));
   } else {
       itsTableSelector=itsTableSelector && (table().col("ANTENNA1") ==
           static_cast<casacore::Int>(ant1)) && (table().col("ANTENNA2") ==
	   static_cast<casacore::Int>(ant2));
   }
}

/// Choose all baselines to given antenna
/// @param[in] ant the sequence number of antenna
void TableScalarFieldSelector::chooseAntenna(casacore::uInt ant)
{
   if (itsTableSelector.isNull()) {
       itsTableSelector=(table().col("ANTENNA1") ==
           static_cast<casacore::Int>(ant)) || (table().col("ANTENNA2") ==
	   static_cast<casacore::Int>(ant));
   } else {
       itsTableSelector=itsTableSelector && ((table().col("ANTENNA1") ==
           static_cast<casacore::Int>(ant)) || (table().col("ANTENNA2") ==
	   static_cast<casacore::Int>(ant)));
   }
}

/// @brief Choose samples corresponding to a uv-distance larger than threshold
/// @details This effectively rejects the baselines giving a smaller
/// uv-distance than the specified threshold
/// @param[in] uvDist threshold
void TableScalarFieldSelector::chooseMinUVDistance(casacore::Double uvDist)
{
  TableExprNode uvwExprNode = table().col("UVW");
  const TableExprNode uExprNode = uvwExprNode(IPosition(1,0));
  const TableExprNode vExprNode = uvwExprNode(IPosition(1,1));
  
  if (itsTableSelector.isNull()) {
      itsTableSelector = (ndim(uvwExprNode) == 1) && (nelements(uvwExprNode) >= 2) 
                    && (sqrt(square(uExprNode)+square(vExprNode)) >= uvDist);
  } else {
      itsTableSelector = itsTableSelector && (ndim(uvwExprNode) == 1) && 
                 (nelements(uvwExprNode) >= 2) && 
                 (sqrt(square(uExprNode) + square(vExprNode)) >= uvDist);
  }
}

/// @brief Choose samples corresponding to either zero uv-distance or larger than threshold
/// @details This effectively rejects the baselines giving a smaller
/// uv-distance than the specified threshold (in metres), but unlike chooseMinUVDistance
/// preserve samples with uvw equal to exatly zero. One example of such zero uvw samples is
/// auto-correlation (which can be filtered out separately by another selector call), but the 
/// main motivation behind such method is to preserve completely flagged samples which may not
/// have uvw defined (and therefore it could be set to zero)
/// @param[in] uvDist threshold
void TableScalarFieldSelector::chooseMinNonZeroUVDistance(casacore::Double uvDist)
{
  TableExprNode uvwExprNode = table().col("UVW");
  const TableExprNode uExprNode = uvwExprNode(IPosition(1,0));
  const TableExprNode vExprNode = uvwExprNode(IPosition(1,1));
  const TableExprNode wExprNode = uvwExprNode(IPosition(1,2));
  
  if (itsTableSelector.isNull()) {
      itsTableSelector = (ndim(uvwExprNode) == 1) && (nelements(uvwExprNode) >= 3) 
                    && ((sqrt(square(uExprNode)+square(vExprNode)) >= uvDist) || 
                    ((uExprNode == 0.) && (vExprNode == 0.) && (wExprNode == 0.)));
  } else {
      itsTableSelector = itsTableSelector && (ndim(uvwExprNode) == 1) && 
                 (nelements(uvwExprNode) >= 3) && 
                 ((sqrt(square(uExprNode) + square(vExprNode)) >= uvDist) || 
                 ((uExprNode == 0.) && (vExprNode == 0.) && (wExprNode == 0.)));
  }
}


/// @brief Choose samples corresponding to a uv-distance smaller than threshold
/// @details This effectively rejects the baselines giving a larger
/// uv-distance than the specified threshold
/// @param[in] uvDist threshold
void TableScalarFieldSelector::chooseMaxUVDistance(casacore::Double uvDist)
{
  TableExprNode uvwExprNode = table().col("UVW");
  const TableExprNode uExprNode = uvwExprNode(IPosition(1,0));
  const TableExprNode vExprNode = uvwExprNode(IPosition(1,1));
  
  if (itsTableSelector.isNull()) {
      itsTableSelector = (ndim(uvwExprNode) == 1) && (nelements(uvwExprNode) >= 2) 
                  && (sqrt(square(uExprNode) + square(vExprNode)) <= uvDist);
  } else {
      itsTableSelector = itsTableSelector && (ndim(uvwExprNode) == 1) &&
                 (nelements(uvwExprNode) > 2) && 
                 (sqrt(square(uExprNode) + square(vExprNode)) <= uvDist);
  }
}


/// Choose a single scan number
/// @param[in] scanNumber the scan number to choose
void TableScalarFieldSelector::chooseScanNumber(casacore::uInt scanNumber)
{
    if (itsTableSelector.isNull()) {
        itsTableSelector = (table().col("SCAN_NUMBER") ==
                static_cast<casacore::Int>(scanNumber));
    } else {
        itsTableSelector = itsTableSelector && (table().col("SCAN_NUMBER") ==
                static_cast<casacore::Int>(scanNumber));
    }
}

/// @brief Choose autocorrelations only
void TableScalarFieldSelector::chooseAutoCorrelations()
{
   if (itsTableSelector.isNull()) {
       itsTableSelector = (table().col("ANTENNA1") ==
                           table().col("ANTENNA2")) &&
                           (table().col("FEED1") ==
                           table().col("FEED2"));
   } else {
       itsTableSelector = itsTableSelector && (table().col("ANTENNA1") ==
                                               table().col("ANTENNA2")) &&
                           (table().col("FEED1") == table().col("FEED2"));
   }
}
  
/// @brief Choose crosscorrelations only
void TableScalarFieldSelector::chooseCrossCorrelations()
{
   if (itsTableSelector.isNull()) {
       itsTableSelector = (table().col("ANTENNA1") !=
                           table().col("ANTENNA2")) || 
                           (table().col("FEED1") != table().col("FEED2"));
   } else {
       itsTableSelector = itsTableSelector && ((table().col("ANTENNA1") !=
                                               table().col("ANTENNA2")) ||
                           (table().col("FEED1") != table().col("FEED2"))) ;
   }
}

/// Choose a single spectral window (also known as IF).
/// @param spWinID the ID of the spectral window to choose
void TableScalarFieldSelector::chooseSpectralWindow(casacore::uInt spWinID)
{
   // one spectral window can correspond to multiple data description IDs
   // We need to obtain this information from the DATA_DESCRIPTION table
   std::vector<size_t> dataDescIDs = subtableInfo().
        getDataDescription().getDescIDsForSpWinID(static_cast<int>(spWinID));
   if (dataDescIDs.size()) {
       std::vector<size_t>::const_iterator ci=dataDescIDs.begin();
       TableExprNode tempNode=(table().col("DATA_DESC_ID") ==
                  static_cast<casacore::Int>(*ci));
       for(++ci;ci!=dataDescIDs.end();++ci) {           
	   tempNode = tempNode || (table().col("DATA_DESC_ID") ==
                  static_cast<casacore::Int>(*ci));
       }
       if (itsTableSelector.isNull()) {
           itsTableSelector=tempNode;
       } else {
           itsTableSelector = itsTableSelector && tempNode;
       }       
   } else {
     // required spectral window is not present in the measurement set
     // we have to insert a dummy expression, otherwise an exception
     // is thrown within the table selection.
     itsTableSelector=(table().col("DATA_DESC_ID") == -1) && False;
   }   
}
 
/// @brief Obtain a table expression node for selection. 
/// @details This method is
/// used in the implementation of the iterator to form a subtable
/// obeying the selection criteria specified by the user via
/// IDataSelector interface
/// @return a const reference to table expression node object
const casacore::TableExprNode& TableScalarFieldSelector::getTableSelector(const
            boost::shared_ptr<IDataConverterImpl const> &) const
{ 
  return itsTableSelector;
}

/// @brief get read-write access to expression node
/// @return a reference to the cached table expression node
///
casacore::TableExprNode& TableScalarFieldSelector::rwTableSelector() const
{ 
  return itsTableSelector;
}

