#ifndef _FILL_INFO_H_
#define _FILL_INFO_H_

enum class FILLER { Single, Dipart, None };

struct FillVals {
  CUTS ePos;
  FILLER type;
  Particle* part = NULL;
  Particle* part2 = NULL;

  FillVals(): type(FILLER::None) {}
  FillVals(CUTS _ePos): ePos(_ePos), type(FILLER::None) {}
  FillVals(CUTS _ePos, FILLER _type, Particle* _part, Particle* _part2=NULL) : 
     ePos(_ePos), type(_type), part(_part), part2(_part2) { }

};

#endif
