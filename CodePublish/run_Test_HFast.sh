#!/bin/bash
chmod u=rwx,g=r,o=r run_Test_HFast.sh
module load intel
# module load impi # Use this module only if you want the same rank for every system
mpic++ -o Test_Fast Main.cpp -O3
mpirun -np 1 ./Test_Fast ammo_0.in rajout
