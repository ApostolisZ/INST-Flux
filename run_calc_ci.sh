#!/bin/bash

NTASK=21

#for i in `seq 1 $NTASK`
#do
#    srun matlab -r "kfp_emu('syn7002_wt_average.csv', 'syn7002_kfp.csv', 10); exit" &
#    sleep 1
#done

for i in `seq 1 $NTASK`
do

    sbatch << 'EOF'
#!/bin/bash
#SBATCH -J calc_ci
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=4096
#SBATCH --cpus-per-task=2
#SBATCH --partition=onegbnet
#SBATCH --qos=micro_dslun

module purge
module load matlab
module load cplex

#matlab -r "calc_ci_samp('cw15+N', 'kfp_emu_20161220T032353', 2); exit"
#matlab -r "calc_ci_samp('cw15-N', 'kfp_emu_20170104T051846', 1); exit"
#matlab -r "calc_ci_samp('sta6+N', 'kfp_emu_20170107T201512', 2); exit"
matlab -r "calc_ci_samp('sta6-N', 'kfp_emu_20170105T045638', 4); exit"
EOF

    sleep 2s;
    echo $i;
done;
