#ifndef WCSimRandomParameters_h
#define WCSimRandomParameters_h 1
#include "WCSimRandomMessenger.hh"
#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/RanluxEngine.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RanecuEngine.h"

class WCSimRandomParameters
{

public:
  WCSimRandomParameters() 
    {
      seed=31415; 
      RandomMessenger = new WCSimRandomMessenger(this);
    }
  ~WCSimRandomParameters() {delete RandomMessenger;}

  int GetSeed() {return CLHEP::HepRandom::getTheSeed();}
  void SetSeed(int iseed) 
    { 
      CLHEP::HepRandom::setTheSeed(iseed);
      printf("Setting the Random Seed to: %d\n",iseed); 
      seed = iseed;
    }


private:
  int seed;
  WCSimRandomMessenger *RandomMessenger;
};

#endif
