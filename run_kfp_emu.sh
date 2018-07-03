#!/bin/bash

NTASK=1

for i in `seq 1 $NTASK`
do

    sbatch << 'EOF'
#!/bin/bash
#SBATCH -J kfp_emu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=4096
#SBATCH --cpus-per-task=2
#SBATCH --partition=onegbnet
#SBATCH --qos=micro_dslun

module purge
module load matlab
module load cplex

matlab -r "kfp_emu('sta6-N', 1); exit"

#matlab -r "myCluster = parcluster('local'); delete(myCluster.Jobs); parpool(8); kfp_emu_ps('syn7002_wt_average-3d', 1); exit"
EOF

    sleep 2s;
    echo $i;
done;
