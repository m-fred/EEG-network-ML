# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 23:02:50 2019

@author: Matthew
"""

###############################################################################

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=2:mem=2gb

module load anaconda3/personal
cp -a $PBS_O_WORKDIR/. $TMPDIR/

declare -i i

for i in {0..100000..1000}; do
    python bashrun6140.py ${i};
    cp $TMPDIR/T_* $PBS_O_WORKDIR/;
    rm $TMPDIR/T_*;
done;



##############################################################################


#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=2:mem=2gb

module load matlab
cp -a $PBS_O_WORKDIR/. $TMPDIR/

declare -i i

for i in {0..100000..1000}; do
    matlab -nodisplay -nodesktop -r "bashrun6140($i)";
    cp $TMPDIR/P_* $PBS_O_WORKDIR/;
    rm $TMPDIR/P_*;
done;

