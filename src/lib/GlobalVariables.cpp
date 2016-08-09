#include <lib/GlobalVariables.hpp>

using namespace gcm;


real Clock::time = 0;
real Clock::timeStep = 0;
int Mpi::rank = 0;
int Mpi::size = 1;
bool Mpi::forceSequence = false;


