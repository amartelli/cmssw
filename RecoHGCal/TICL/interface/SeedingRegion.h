#ifndef __RecoHGCal_TICL_SeedingRegion_H__
#define __RecoHGCal_TICL_SeedingRegion_H__

#include <array>
#include <vector>

struct SeedingRegion {

  int etaStart;
  int etaLength;
  int phiStart;
  int phiLength;
  int zSide;

  SeedingRegion( int etaS, int etaL, int phiS, int phiL, int zS){
    etaStart = etaS;
    etaLength = etaL;
    phiStart = phiS;
    phiLength = phiL;
    zSide = zS;
  }

};

#endif
