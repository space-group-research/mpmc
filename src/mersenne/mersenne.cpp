#include <random>
extern
    "C" {
#include <mc.h>
}
 
static std::random_device rd; // random device engine, usually based on /dev/random on UNIX-like systems
// initialize Mersennes' twister using rd to generate the seed
static std::mt19937 rng(rd());

double get_rand()
{
    static std::uniform_real_distribution<double> uid(0.0,1.0); // random dice
    return uid(rng); // use rng as a generator
}
