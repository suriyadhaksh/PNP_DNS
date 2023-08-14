rm *.plt
#touch Energy_2D.plt
#echo " VARIABLES = \"t" "mass" "bulk" "interface" "total\" "  > Energy_2D.plt;
mpirun -n 8 /home/suriya/PROJECTS/PNP_FullResidual_MMSGF/build/pnp ksp_rtol 1E-8 -pc_type asm -ksp_type bcgs -snes_monitor -snes_converged_reason
#-ksp_monitor