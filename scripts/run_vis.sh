#!/bin/bash      
#SBATCH -J vis -p development -N 1 --ntasks-per-node 56 -t 2:00:00 -A AST21002
source $HOME/.bashrc
python ~/scripts/CrunchSnaps/scripts/make_movie_from_camerafile.py camerafile.txt ../../M2e4_R10_S0_T1_B0.01_Res271_n2_sol0.5_42/output --np=14 --np_render=4 --res=1920 --no_timestamp --fresco_stars
