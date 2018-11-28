#include <random>
extern
    "C" {
#include <mc.h>
}

static std::mt19937 rng;

double get_rand(system_t *system)
{
  if (!system->initialized){
    system->initialized = 1;
    if (system->preset_seeds_on){
      rng.seed(system->preset_seeds);
    }
    else {
      static std::random_device rd;
      rng.seed(rd());
    }
  }
    static std::uniform_real_distribution<double> uid(0.0,1.0); // random dice
    return uid(rng); // use rng as a generator
}
