#!/bin/bash

scf=slurm_script
DATA_DIR_OUT=NE/Vg_25/ #NE/Katanin/Vary_U/ #NE/Katanin/Vg_25/ #NE/Katanin/Vary_U/ #NE/Katanin/Vg_25/
hl=`seq -0.00001 0.00001 0.00001`
height=`seq .22 .005 0.29`
sites=`seq 0 1 60`
interaction=`seq -5 1 -1`
#interaction=`seq -1 1 -1`
chpot=`seq -1.475 -0.005 -1.575`
delta=`seq .0 0.001 .0`
len=`seq 61 20 161`

#mkdir -p $HOME/DATA/$DATA_DIR_OUT
mkdir -p /naslx/projects/uh3o1/ri26yad/DATA/$DATA_DIR_OUT
for h in .000 #hl #`./my_exp $hl`             #magnetic field
 do
  for a in .25 #.94   #height of potential
   do
    for U0 in .55 #`./my_exp $interaction`
     do
      for chavg in -1.475 #$chpot
       do
        for deltaCh in $delta #$chpot #.1 #$delta #$chpot #-1.475 #$chpot
         do
          #cat options > $scf
          cat SLURM_OPTIONS > $scf

          muL="$(echo "$chavg + $deltaCh" | bc)"
          muR="$(echo "$chavg - $deltaCh" | bc)"

          echo PROGRAM=\"$HOME/bin/dfRG2_general_B0_clean\" >> $scf
          #echo PROGRAM=\"$HOME/bin/dfRG_lin_mod\" >> $scf
          echo OPTIONS=\"-muL $muL -muR $muR -N 31 -Nff 1400 -Nfb 1400 -U $U0 -h $h -Vg $a -TL 0.001 -TR 0.001 -pot_type 0 -save_all 1 -save_dens 1\" >> $scf
          #echo OPTIONS=\"-mu $mu -N 101 -Nff 150 -Nfb 150 -U $U0 -h $h -Vg $a -site 30 -save_all 1 -save_dens 1\" >> $scf
          echo DIR_DATA_OUT=\"$WORK/DATA/$DATA_DIR_OUT\" >> $scf

          cat debuginfo >> $scf
          
          echo "cd \$DIR_DATA_OUT/" >> $scf
          echo "date" >> $scf
          echo "time \$PROGRAM \$OPTIONS" >> $scf
          echo "date" >> $scf

          sbatch slurm_script
          cd ..
        done
      done
    done
  done
done
