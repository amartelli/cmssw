#ifndef __RecoHGCal_TICL_SeedingRegion_H__
#define __RecoHGCal_TICL_SeedingRegion_H__

#include <array>
#include <vector>
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

struct SeedingRegion {

  int etaStart;
  int etaLength;
  int phiStart;
  int phiLength;
  int zSide;

  GlobalPoint point;
  GlobalVector direction;
  int index;
  int etaBin;
  int phiBin;

  SeedingRegion(){}

  SeedingRegion( int etaS, int etaL, int phiS, int phiL, int zS,
		 GlobalPoint p, GlobalVector d, int idx, int etaB, int phiB){
    etaStart = etaS;
    etaLength = etaL;
    phiStart = phiS;
    phiLength = phiL;
    zSide = zS;

    point = p;
    direction = d;
    index = idx;
    etaBin = etaB;
    phiBin = phiB;
  }

};

#endif
