#ifndef __RecoHGCal_TICL_PropagationSeedingPoint_H__
#define __RecoHGCal_TICL_PropagationSeedingPoint_H__

#include <array>
#include <vector>

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

struct PropagationSeedingPoint {

  GlobalPoint p;
  GlobalVector v;
  int index;
  int etaBin;
  int phiBin;
};

#endif
