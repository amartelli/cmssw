#include "RecoLocalCalo/HGCalRecAlgos/interface/HGCalDirect3DClustering.h"

// Geometry
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"

#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"
//
#include "DataFormats/CaloRecHit/interface/CaloID.h"

void HGCalDirect3DClustering::populate(const HGCRecHitCollection& hits){

  std::cout << "verbosity = " << verbosity<< std::endl;

  if(verbosity < hgcal::pINFO){
    std::cout << " **** HGCalDirect3DClustering::populate **** " << std::endl;
  }


  //loop over all hits and create the Hexel structure, skip energies below ecut

  if (dependSensor) {
  // for each layer and wafer calculate the thresholds (sigmaNoise and energy) once
    computeThreshold();
  }

  std::vector<bool> firstHit(2*(maxlayer+1), true);
  for (unsigned int i=0;i<hits.size();++i) {

    const HGCRecHit& hgrh = hits[i];
    DetId detid = hgrh.detid();
    unsigned int layer = rhtools_.getLayerWithOffset(detid);
    float thickness = 0.;
    // set sigmaNoise default value 1 to use kappa value directly in case of sensor-independent thresholds
    float sigmaNoise = 1.;
    if(dependSensor){
      if (layer<= lastLayerFH) // only EE and FH have silicon sensors
        thickness = rhtools_.getSiThickness(detid);
      double storedThreshold=thresholds[layer-1][layer<=lastLayerFH ? rhtools_.getWafer(detid) : 0];
      sigmaNoise = v_sigmaNoise[layer-1][layer<=lastLayerFH ? rhtools_.getWafer(detid) : 0];

      if(hgrh.energy() <  storedThreshold)
        continue; //this sets the ZS threshold at ecut times the sigma noise for the sensor
    }
    if(!dependSensor && hgrh.energy() < ecut) continue;

    // map layers from positive endcap (z) to layer + maxlayer+1 to prevent mixing up hits from different sides
    layer += int(rhtools_.zside(detid)>0)*(maxlayer+1);

    // determine whether this is a half-hexagon
    bool isHalf = rhtools_.isHalfCell(detid);
    const GlobalPoint position( std::move( rhtools_.getPosition( detid ) ) );

    //here's were the KDNode is passed its dims arguments - note that these are *copied* from the Hexel
    //define x,y,z coordinates as in the Hexel

    if(thickness<0.)thickness = 0.;
    unsigned int zside = unsigned(position.z()>0);
    float reducedz = float(rhtools_.getLayerWithOffset(detid))/2.;
    //RA
    //need to change reducedz as wrt dE/dX weighting*
    points[zside].pts.emplace_back(hgrh,detid,isHalf,sigmaNoise,thickness,&rhtools_);

    if(points[zside].pts.size()==0){
      minpos[zside][0] = position.x(); minpos[zside][1] = position.y(); minpos[zside][2] = reducedz;
      maxpos[zside][0] = position.x(); maxpos[zside][1] = position.y(); maxpos[zside][2] = reducedz;
    }else{
      minpos[zside][0] = std::min((float)position.x(),minpos[zside][0]);
      minpos[zside][1] = std::min((float)position.y(),minpos[zside][1]);
      minpos[zside][2] = std::min((float)reducedz,minpos[zside][2]);
      maxpos[zside][0] = std::max((float)position.x(),maxpos[zside][0]);
      maxpos[zside][1] = std::max((float)position.y(),maxpos[zside][1]);
      maxpos[zside][1] = std::max((float)reducedz,maxpos[zside][2]);
    }
  } // end loop hits
}



// Create a vector of Hexels associated to one cluster from a collection of HGCalRecHits - this can be used
// directly to make the final cluster list - this method can be invoked multiple times for the same event
// with different input (reset should be called between events)
void HGCalDirect3DClustering::makeClusters()
{

  if(verbosity < hgcal::pINFO){
    std::cout << " **** building HGCalDirect3DClustering **** " << std::endl;
    std::cout << " delta_c " << vecDeltas[0] << ", " << vecDeltas[1] << ", " << vecDeltas[2] << " kappa " << kappa << std::endl;
  }


  //assign all hits in each layer to a cluster core or halo
  // just 2 blocks one per endcap
  for (unsigned int i = 0; i < 2; ++i) {

    my_kd_tree_t hit_kdtree(3, points[i], nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
    hit_kdtree.buildIndex();

    double maxdensity = calculateLocalDensity(points[i], hit_kdtree); // also stores rho (energy density) for each point (node)
    // calculate distance to nearest point with higher density storing distance (delta) and point's index
    calculateDistanceToHigher(points[i]);
    findAndAssignClusters(points[i], hit_kdtree, maxdensity);
  }
  //make the cluster vector
}


//RA tocheck
std::vector<reco::BasicCluster> HGCalDirect3DClustering::getClusters(bool doSharing){

  if(verbosity < hgcal::pINFO){
    std::cout << " **** HGCalDirect3DClustering::getClusters **** " << std::endl;
  }

  reco::CaloID caloID = reco::CaloID::DET_HGCAL_ENDCAP;
  std::vector< std::pair<DetId, float> > thisCluster;
  for (unsigned int i = 0; i < current_v.size(); ++i){
    double energy = 0;
    Point position;

    if( doSharing ) {
      std::vector<unsigned> seeds = findLocalMaximaInCluster(current_v[i] );
      // sharing found seeds.size() sub-cluster seeds in cluster i

      std::vector<std::vector<double> > fractions;
      // first pass can have noise it in
      shareEnergy(current_v[i],seeds,fractions);

      // reset and run second pass after vetoing seeds
      // that result in trivial clusters (less than 2 effective cells)

      for( unsigned isub = 0; isub < fractions.size(); ++isub ) {
	double effective_hits = 0.0;
	double energy  = calculateEnergyWithFraction(current_v[i],fractions[isub]);
	Point position = calculatePositionWithFraction(current_v[i],fractions[isub]);

	for( unsigned ihit = 0; ihit < fractions[isub].size(); ++ihit ) {
	  const double fraction = fractions[isub][ihit];
	  if( fraction > 1e-7 ) {
	    effective_hits += fraction;
	    thisCluster.emplace_back(current_v[i][ihit].detid,fraction);
	  }
	}

	if (verbosity < hgcal::pINFO)
	  {
	    std::cout << "\t******** NEW CLUSTER (SHARING) ********" << std::endl;
	    std::cout << "\tEff. No. of cells = " << effective_hits << std::endl;
	    std::cout << "\t     Energy       = " << energy << std::endl;
	    std::cout << "\t     Phi          = " << position.phi() << std::endl;
	    std::cout << "\t     Eta          = " << position.eta() << std::endl;
	    std::cout << "\t*****************************" << std::endl;
	  }
	clusters_v.push_back(reco::BasicCluster(energy, position, caloID, thisCluster,
						algoId));
	thisCluster.clear();
      }
    }
    else{
      position = calculatePosition(current_v[i]); // energy-weighted position
      
      for (auto& it: current_v[i])
	{
	  energy += it.isHalo ? 0. : it.weight;
	  // use fraction to store whether this is a Halo hit or not
	  thisCluster.emplace_back(std::pair<DetId, float>(it.detid,(it.isHalo?0.:1.)));
	};
      if (verbosity < hgcal::pINFO)
	{
	  std::cout << "******** NEW CLUSTER (HGCIA) ********" << std::endl;
	  std::cout << "Index          " << i                   << std::endl;
	  std::cout << "No. of cells = " << current_v[i].size() << std::endl;
	  std::cout << "     Energy     = " << energy << std::endl;
	  std::cout << "     Phi        = " << position.phi() << std::endl;
	  std::cout << "     Eta        = " << position.eta() << std::endl;
	  std::cout << "*****************************" << std::endl;
	}
      clusters_v.push_back(reco::BasicCluster(energy, position, caloID, thisCluster,
					      algoId));
      thisCluster.clear();
    }
  }
  return clusters_v;
}



math::XYZPoint HGCalDirect3DClustering::calculatePosition(std::vector<Hexel> &v){


  if(verbosity < hgcal::pINFO){
    std::cout << " **** HGCalDirect3DClustering::calculatePosition **** " << std::endl;
  }

  float total_weight = 0.;
  float x = 0.;
  float y = 0.;
  float z = 0.;
  unsigned int v_size = v.size();
  unsigned int maxEnergyIndex = 0;
  float maxEnergyValue = 0;
  bool haloOnlyCluster = true;
  
  // loop over hits in cluster candidate building up weight for
  // energy-weighted position calculation and determining the maximum
  // energy hit in case this is a halo-only cluster
  for (unsigned int i = 0; i < v_size; i++){
    if(!v[i].isHalo){
      haloOnlyCluster = false;
      total_weight += v[i].weight;
      x += v[i].x*v[i].weight;
      y += v[i].y*v[i].weight;
      z += v[i].z*v[i].weight;
    }
    else {
      if (v[i].weight > maxEnergyValue) {
        maxEnergyValue = v[i].weight;
        maxEnergyIndex = i;
      }
    }
  }

  if (!haloOnlyCluster) {
    if (total_weight != 0) {
      return math::XYZPoint( x/total_weight,
			     y/total_weight,
			     z/total_weight );
    }
  }
  else if (v_size > 0) {
    // return position of hit with maximum energy
    return math::XYZPoint(v[maxEnergyIndex].x, v[maxEnergyIndex].y, v[maxEnergyIndex].z);
  }
  return math::XYZPoint(0, 0, 0);
}


//compute local energy density for each point node (nd_i) in the tree structure (lp)
double HGCalDirect3DClustering::calculateLocalDensity(HexelCloud &nd, my_kd_tree_t &lp){

  if(verbosity < hgcal::pINFO){
    std::cout << " **** HGCalDirect3DClustering::calculateLocalDensity **** " << std::endl;
  }

  double maxdensity = 0.;

  // delta_c is maximum search distance (critical distance) for local density calculation
  // make it subdetector dependent
  float delta_c; 

  for(unsigned int i = 0; i < nd.pts.size(); ++i){
    //RA would give reduced or effective z instead of lz
    const double query_pt[3] = {nd.pts[i].x,nd.pts[i].y,nd.pts[i].lz};

    if( nd.pts[i].lz*2.<=lastLayerEE ) delta_c = vecDeltas[0];
    else if( nd.pts[i].lz*2.<=lastLayerFH ) delta_c = vecDeltas[1];
    else delta_c = vecDeltas[2];

    double search_radius = static_cast<double>(delta_c*delta_c);
    nanoflann::SearchParams params;

    // prepare container for index and distance of neighbors within radius                                                                        
    std::vector<std::pair<size_t,double> >  ret_matches; 
    const size_t nMatches = lp.radiusSearch(&query_pt[0],search_radius, ret_matches, params);

    if (verbosity < hgcal::pINFO){
      std::cout << "hit " << i << " found " << nMatches << " nearest neighbors " << std::endl;
    }

    nd.pts[i].nNeighbors = nMatches;
    //neighbours = number of points in the search area
    for(unsigned int k = 0; k < ret_matches.size(); k++){
      unsigned int j = ret_matches[k].first;
      if (verbosity < hgcal::pINFO){ 
	std::cout << "hit i(node),j(found) " << i << "," << j
		  << " distance by Hexel " << distance2(nd.pts[i],nd.pts[j]) 
		  << " distance by kdtree " <<  ret_matches[k].second << std::endl;
      std::cout << "delta x " << nd.pts[i].x-nd.pts[j].x << std::endl;
      std::cout << "delta y " << nd.pts[i].y-nd.pts[j].y << std::endl;
      std::cout << "delta lz " << nd.pts[i].lz-nd.pts[j].lz << std::endl;
    }
      //here compute the local energy density => possible to specify the Kernel for rho
      //by construction distance is < delta_c
      //if(distance(nd[i].data,found[j].data) < delta_c)
      // flat kernel
      nd.pts[i].rho += nd.pts[j].weight;
      //RA
      //Gaus kernel 
    }//end loop over found neighbours

    if(nd.pts[i].rho > maxdensity) maxdensity = nd.pts[i].rho;

    if (verbosity < hgcal::pINFO){
      std::cout << "hit " << i << " local density is " << nd.pts[i].rho << std::endl;
    }
  } //end loop over nodes
  return maxdensity;
}


//order nodes by density
//for each node (start from the 2nd) save
//in delta the smallest distance wrt the hits with higher density
//in nearestHigher the hit with higher density having the smallest distance
double HGCalDirect3DClustering::calculateDistanceToHigher(HexelCloud &nd){

  if(verbosity < hgcal::pINFO){
    std::cout << " **** HGCalDirect3DClustering::calculateDistanceToHigher **** " << std::endl;
  }


  //sort vector of Hexels by decreasing local density
  std::vector<size_t> rs = sorted_indices(nd.pts);

  double maxdensity = 0.0;
  int nearestHigher = -1;

  //either first element or there are no hits
  if(!rs.empty())
    maxdensity = nd.pts[rs[0]].rho;
  else
    return maxdensity;

  double dist2 = 0.;

  //start by setting delta for the highest density hit to
  //the most distant hit - this is a convention
  for (auto& j: nd.pts) {
    double tmp = distance2(nd.pts[rs[0]], j);
    if(tmp > dist2) dist2 = tmp;
  }
  nd.pts[rs[0]].delta = std::sqrt(dist2);
  nd.pts[rs[0]].nearestHigher = nearestHigher;

  //now we save the largest distance as a starting point
  const double max_dist2 = dist2;
  const unsigned int nd_size = rs.size();

  if(verbosity < hgcal::pINFO){
    std::cout << " number 0 delta = " << nd.pts[rs[0]].delta << " nearestHigher = " << nearestHigher << std::endl;
  }

  //start from second-highest density
  for(unsigned int oi = 1; oi < nd_size; ++oi){
    dist2 = max_dist2;
    unsigned int i = rs[oi];
    // we only need to check up to oi since hits
    // are ordered by decreasing density
    // and all points coming BEFORE oi are guaranteed to have higher rho
    // and the ones AFTER to have lower rho
    for(unsigned int oj = 0; oj < oi; ++oj){
      unsigned int j = rs[oj];
      double tmp = distance2(nd.pts[i], nd.pts[j]);
      if(tmp <= dist2){ //this "<=" instead of "<" addresses the (rare) case when there are only two hits
	dist2 = tmp;
	nearestHigher = j;
      }
    }
    nd.pts[i].delta = std::sqrt(dist2);
    nd.pts[i].nearestHigher = nearestHigher; //this uses the original unsorted hitlist

    if (verbosity < hgcal::pINFO){
      std::cout << "hit " << i << " distance to higher is " << nd.pts[i].delta << std::endl;
      std::cout << "hit " << i << " nearest higher is " << nd.pts[i].nearestHigher << std::endl;
    }
  }
  return maxdensity;
}



int HGCalDirect3DClustering::findAndAssignClusters(HexelCloud &nd,my_kd_tree_t &lp, double maxdensity){

  if(verbosity < hgcal::pINFO){
    std::cout << " **** HGCalDirect3DClustering::findAndAssignClusters **** " << std::endl;
  }

  nanoflann::SearchParams params;
  params.sorted = true;

  //this is called once per side... 
  //so when filling the cluster temporary vector of Hexels we resize each time by the number
  //of clusters found. This is always equal to the number of cluster centers...

  unsigned int clusterIndex = 0;
  float delta_c; // critical distance

  std::vector<size_t> rs = sorted_indices(nd.pts); // indices sorted by decreasing rho
  std::vector<size_t> ds = sort_by_delta(nd.pts); // sort in decreasing distance to higher

  const unsigned int nd_size = nd.pts.size();
  for(unsigned int i=0; i < nd_size; ++i){

    if( nd.pts[ds[i]].lz*2.<=lastLayerEE ) delta_c = vecDeltas[0];
    else if( nd.pts[ds[i]].lz*2.<=lastLayerFH ) delta_c = vecDeltas[1];
    else delta_c = vecDeltas[2];


    //look for isolated centers =>  need minimum distance delta_C
    if(nd.pts[ds[i]].delta < delta_c) break; // no more cluster centers to be looked at

    //look for high energy centers => minimum required rho_c
    if(dependSensor){

      float rho_c = kappa*nd.pts[ds[i]].sigmaNoise;
      if(nd.pts[ds[i]].rho < rho_c ) continue; // set equal to kappa times noise threshold

    }
    else if(nd.pts[ds[i]].rho*kappa < maxdensity)
      continue;

    //for each node clusterIndex is the number of the cluster to which it belongs
    //-1 if does not 
    nd.pts[ds[i]].clusterIndex = clusterIndex;
    if (verbosity < hgcal::pINFO)
      {
	std::cout << "Adding new cluster with index " << clusterIndex+cluster_offset << std::endl;
	std::cout << "Cluster center is hit " << ds[i] << std::endl;
      }
    clusterIndex++;
  }

  //at this point clusterIndex is equal to the number of cluster centers - if it is zero we are
  //done
  if(clusterIndex==0) return clusterIndex;

  //assign remaining points (that have clusterIndex -1 = initialization value) to clusters, 
  //using the nearestHigher set from previous step (always set except
  // for top density hit that is skipped...)
  for(unsigned int oi =1; oi < nd_size; ++oi){
    unsigned int i = rs[oi];
    int ci = nd.pts[i].clusterIndex;
    if(ci == -1){ // clusterIndex is initialised with -1 if not yet used in cluster
      nd.pts[i].clusterIndex =  nd.pts[nd.pts[i].nearestHigher].clusterIndex;
    }
  }

  //make room in the temporary cluster vector for the additional clusterIndex clusters
  // from this layer
  if (verbosity < hgcal::pINFO)
    {
      std::cout << "resizing cluster vector by "<< clusterIndex << std::endl;
    }
  current_v.resize(cluster_offset+clusterIndex);

  //border hits used to define border density and then halo  
  //for each node-center assign the points closer than dc to other clusters to border region
  //also assign isolated hits as border ??? why? 
  //and find critical border density
  std::vector<double> rho_b(clusterIndex,0.);
  // no need to clear previous results lp.clear();
  // within radiusSearch: Previous contents of \a IndicesDists are cleared

  //now loop on all hits again :( and check: if there are hits from another cluster within d_c -> flag as border hit
  for(unsigned int i = 0; i < nd_size; ++i){
    int ci = nd.pts[i].clusterIndex;
    bool flag_isolated = true;
    if(ci != -1){

      double search_radius;
      if( nd.pts[ds[i]].lz*2.<=lastLayerEE ) delta_c = vecDeltas[0];
      else if( nd.pts[ds[i]].lz*2.<=lastLayerFH ) delta_c = vecDeltas[1];
      else delta_c = vecDeltas[2];
      search_radius = static_cast<double>(delta_c*delta_c);

      //RA update lz to reducedz 
      const double query_pt[3] = {nd.pts[i].x,nd.pts[i].y,nd.pts[i].lz};
      std::vector<std::pair<size_t,double> >  ret_matches; //index and distance of neighbors within radius
      const size_t nMatches = lp.radiusSearch(&query_pt[0],search_radius, ret_matches, params);

      if(nMatches!=nd.pts[i].nNeighbors && verbosity < hgcal::pERROR){
	std::cout << "for hit " << i << " of cluster index "<< ci 
		  << " neighbors found differ " << nMatches << " - " << nd.pts[i].nNeighbors
                  << std::endl;
      }

      //the matches are sorted (closest point first) 
      //nMatches == 1 means 1 cluster with 1! hit
      if(nMatches != 1){
	//security check: each hit is the closest to itself
        assert(i==ret_matches[0].first);

	//start from closest hit (exclude the hit itself)
        for(unsigned int k = 1; k < ret_matches.size(); k++){
          unsigned int j = ret_matches[k].first;
	  //check existance of another hit within delta_c associated to the same cluster => non_isolated
          if(nd.pts[j].clusterIndex!=-1 && nd.pts[j].clusterIndex==ci){
            flag_isolated = false;
            break;
          }
	  //in case of hits within delta_c, belonging to other clusters, and with higher density => mark as border
          if(nd.pts[j].clusterIndex!=-1 && nd.pts[j].clusterIndex!=ci && nd.pts[j].rho>nd.pts[i].rho){
            if (verbosity < hgcal::pINFO)
              {
		std::cout << "for hit " << i << " of cluster index "<< ci << " found nearest hit at "
                          << nd.pts[j].x << ","
                          << nd.pts[j].y << ","
                          << nd.pts[j].lz << " with cluster index " << nd.pts[j].clusterIndex
                          << " rho[i] " << nd.pts[i].rho
                          << " rho[j] " << nd.pts[j].rho << std::endl;
              }
            //in which case we assign it to the border                                                               
            nd.pts[i].isBorder = true;
          }
        }

      }
      //if single hit also mark as border
      if(flag_isolated) nd.pts[i].isBorder = true; //the hit is more than delta_c from any of its brethren
    
      //check if this border hit has density larger than the current rho_b and update
      if(nd.pts[i].isBorder && rho_b[ci] < nd.pts[i].rho){
        if (verbosity < hgcal::pINFO)
          {
	    std::cout << "for hit " << i << " of cluster index "<< ci << " updating rho_b: before " <<  rho_b[ci]
                      << " after " << nd.pts[i].rho << std::endl;
          }
        rho_b[ci] = nd.pts[i].rho;
      }      
    }
  }

  //flag points in cluster with density < rho_b as halo points, then fill the cluster vector
  for(unsigned int i = 0; i < nd.pts.size(); ++i){
    int ci = nd.pts[i].clusterIndex;
    if(ci != -1){
      if(nd.pts[i].rho <= rho_b[ci]) nd.pts[i].isHalo = true;
      current_v[ci+cluster_offset].push_back(nd.pts[i]);
      if (verbosity < hgcal::pINFO)
	{
	  std::cout << "Pushing hit " << i << " into cluster with index " << ci+cluster_offset << std::endl;
	  std::cout << "Size now " << current_v[ci+cluster_offset].size() << std::endl;
	}
    }
  }

  //prepare the offset for the next layer if there is one
  if (verbosity < hgcal::pINFO)
    {
      std::cout << "moving cluster offset by " << clusterIndex << std::endl;
    }
  cluster_offset += clusterIndex;

  return clusterIndex;
}


// find local maxima within delta_c, marking the indices in the cluster
std::vector<unsigned> HGCalDirect3DClustering::findLocalMaximaInCluster(const std::vector<Hexel> &v){

  if(verbosity < hgcal::pINFO){
    std::cout << " **** HGCalDirect3DClustering::findLocalMaximaInCluster **** " << std::endl;
  }


  std::vector<unsigned> result;
  std::vector<bool> seed(v.size(),true);
  float delta_c = 2.;

  for(unsigned int i = 0; i < v.size(); ++i){
    for( unsigned j = 0; j < v.size(); ++j ) {
      if( distance(v[i],v[j]) < delta_c && i != j) {
	if(v[i].weight < v[j].weight ) {
	  seed[i] = false;
	  break;
	}
      }
    }
  }

  for(unsigned int i = 0; i < v.size(); ++i){
    if( seed[i] && v[i].weight > 5e-4) {
      // seed at i with energy cluster[i].weight
      result.push_back(i);
    }
  }
  // Found result.size() sub-clusters in input cluster of length cluster.size()
  return result;
}


math::XYZPoint HGCalDirect3DClustering::calculatePositionWithFraction(const std::vector<Hexel> &hits,
								      const std::vector<double>& fractions) {
  double norm(0.0), x(0.0), y(0.0), z(0.0);
  for( unsigned i = 0; i < hits.size(); ++i ) {
    const double weight = fractions[i]*hits[i].weight;
    norm += weight;
    x += weight*hits[i].x;
    y += weight*hits[i].y;
    z += weight*hits[i].z;
  }
  math::XYZPoint result(x,y,z);
  double norm_inv = 1.0/norm;
  result *= norm_inv;
  return result;
}


double HGCalDirect3DClustering::calculateEnergyWithFraction(const std::vector<Hexel> &hits,
							    const std::vector<double>& fractions) {


  if(verbosity < hgcal::pINFO){
    std::cout << " **** HGCalDirect3DClustering::calculateEnergyWithFraction **** " << std::endl;
  }

  double result = 0.0;
  for( unsigned i = 0 ; i < hits.size(); ++i ) {
    result += fractions[i]*hits[i].weight;
  }
  return result;
}


void HGCalDirect3DClustering::shareEnergy(const std::vector<Hexel> &incluster, 
					  const std::vector<unsigned>& seeds,
					  std::vector<std::vector<double> >& outclusters) {

  if(verbosity < hgcal::pINFO){
    std::cout << " **** HGCalDirect3DClustering::shareEnergy **** " << std::endl;
  }


  std::vector<bool> isaseed(incluster.size(),false);
  outclusters.clear();
  outclusters.resize(seeds.size());
  std::vector<Point> centroids(seeds.size());
  std::vector<double> energies(seeds.size());
  
  if( seeds.size() == 1 ) { // short circuit the case of a lone cluster
    outclusters[0].clear();
    outclusters[0].resize(incluster.size(),1.0);
    return;
  }

  // saving seeds

  // create quick seed lookup
  for( unsigned i = 0; i < seeds.size(); ++i ) {
    isaseed[seeds[i]] = true;
  }

  // initialize clusters to be shared
  // centroids start off at seed positions
  // seeds always have fraction 1.0, to stabilize fit
  // initializing fit
  for( unsigned i = 0; i < seeds.size(); ++i ) {
    outclusters[i].resize(incluster.size(),0.0);
    for( unsigned j = 0; j < incluster.size(); ++j ) {
      if( j == seeds[i] ) {
	outclusters[i][j] = 1.0;
	centroids[i] = math::XYZPoint(incluster[j].x,incluster[j].y,incluster[j].z);
	energies[i]  = incluster[j].weight;
      }
    }
  }

  // run the fit while we are less than max iterations, and clusters are still moving
  const double minFracTot = 1e-20;
  unsigned iter = 0;
  const unsigned iterMax = 50;
  double diff = std::numeric_limits<double>::max();
  const double stoppingTolerance = 1e-8;
  const double toleranceScaling = std::pow(std::max(1.0,seeds.size()-1.0),2.0);
  std::vector<Point> prevCentroids;
  std::vector<double> frac(seeds.size()), dist2(seeds.size());
  while( iter++ < iterMax && diff > stoppingTolerance*toleranceScaling ) {
    for( unsigned i = 0; i < incluster.size(); ++i ) {
      const Hexel& ihit = incluster[i];
      double fracTot(0.0);
      for( unsigned j = 0; j < seeds.size(); ++j ) {
	double fraction = 0.0;
	double d2 = ( std::pow(ihit.x - centroids[j].x(),2.0) +
		      std::pow(ihit.y - centroids[j].y(),2.0) +
		      std::pow(ihit.z - centroids[j].z(),2.0)   )/sigma2;
	dist2[j] = d2;
	// now we set the fractions up based on hit type
	if( i == seeds[j] ) { // this cluster's seed
	  fraction = 1.0;
	} else if( isaseed[i] ) {
	  fraction = 0.0;
	} else {
	  fraction = energies[j]*std::exp( -0.5*d2 );
	}
	fracTot += fraction;
	frac[j] = fraction;
      }
      // now that we have calculated all fractions for all hits
      // assign the new fractions
      for( unsigned j = 0; j < seeds.size(); ++j ) {
	if( fracTot > minFracTot ||
	    ( i == seeds[j] && fracTot > 0.0 ) ) {
	  outclusters[j][i] = frac[j]/fracTot;
	} else {
	  outclusters[j][i] = 0.0;
	}
      }
    }
    
    // save previous centroids
    prevCentroids = std::move(centroids);
    // finally update the position of the centroids from the last iteration
    centroids.resize(seeds.size());
    double diff2 = 0.0;
    for( unsigned i = 0; i < seeds.size(); ++i ) {
      centroids[i] = calculatePositionWithFraction(incluster,outclusters[i]);
      energies[i]  = calculateEnergyWithFraction(incluster,outclusters[i]);
      // calculate convergence parameters
      const double delta2 = (prevCentroids[i]-centroids[i]).perp2();
      if( delta2 > diff2 ) diff2 = delta2;
    }
    //update convergance parameter outside loop
    diff = std::sqrt(diff2);
  }
}


void HGCalDirect3DClustering::computeThreshold() {

  if(verbosity < hgcal::pINFO){
    std::cout << " **** HGCalDirect3DClustering::computeThreshold **** " << std::endl;
  }


  if(initialized) return; // only need to calculate thresholds once
  const std::vector<DetId>& listee(rhtools_.getGeometry()->getValidDetIds(DetId::Forward,ForwardSubdetector::HGCEE));
  const std::vector<DetId>& listfh(rhtools_.getGeometry()->getValidDetIds(DetId::Forward,ForwardSubdetector::HGCHEF));

  std::vector<double> dummy;
  dummy.resize(maxNumberOfWafersPerLayer, 0);
  thresholds.resize(maxlayer, dummy);
  v_sigmaNoise.resize(maxlayer, dummy);
  int previouswafer=-999;

  for(unsigned icalo=0;icalo<2;++icalo)
    {
      const std::vector<DetId>& listDetId( icalo==0 ? listee : listfh);

      for(auto& detid: listDetId)
	{
	  int wafer = rhtools_.getWafer(detid);
	  if(wafer==previouswafer) continue;
	  previouswafer = wafer;
	  // no need to do it twice
	  if(rhtools_.zside(detid)<0) continue;
	  int layer = rhtools_.getLayerWithOffset(detid);
	  float thickness = rhtools_.getSiThickness(detid);
	  int thickIndex = -1;
	  if( thickness>99. && thickness<101.) thickIndex=0;
	  else if( thickness>199. && thickness<201. ) thickIndex=1;
	  else if( thickness>299. && thickness<301. ) thickIndex=2;
	  else assert( thickIndex>0 && "ERROR - silicon thickness has a nonsensical value" );
	  float sigmaNoise = 0.001 * fcPerEle * nonAgedNoises[thickIndex] * dEdXweights[layer] / (fcPerMip[thickIndex] * thicknessCorrection[thickIndex]);
	  thresholds[layer-1][wafer]=sigmaNoise*ecut;
	  v_sigmaNoise[layer-1][wafer] = sigmaNoise;
	}
    }

  // now BH, much faster
  for ( unsigned ilayer=lastLayerFH+1;ilayer<=maxlayer;++ilayer)
    {
      float sigmaNoise = 0.001 * noiseMip * dEdXweights[ilayer];
      std::vector<double> bhDummy_thresholds;
      std::vector<double> bhDummy_sigmaNoise;
      bhDummy_thresholds.push_back(sigmaNoise*ecut);
      bhDummy_sigmaNoise.push_back(sigmaNoise);
      thresholds[ilayer-1]=bhDummy_thresholds;
      v_sigmaNoise[ilayer-1]=bhDummy_sigmaNoise;
    }
  initialized=true;
}
