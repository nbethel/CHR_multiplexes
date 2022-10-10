#!/bin/bash
#SBATCH -p cpu
#SBATCH -t 00:12:00
#SBATCH --mem=4g
#SBATCH -o task_log

LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p two_helix_tasks)
${LINE}
/software/rosetta/latest/bin/rosetta_scripts.hdf5.linuxgccrelease  -nstruct 10 -out:prefix rep"$SLURM_ARRAY_TASK_ID"_ -parser:protocol ../gen_backbones.xml  @../flags -remodel:staged_sampling true

#sbatch -a 1-$(cat two_helix_tasks | wc -l) genBB.sh

