#!/bin/bash
#SBATCH --job-name=ntm_sweep
#SBATCH --partition=nodes
#SBATCH --array=0-319
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/sweep_%A_%a.out

cd /home/bimstepf/endor/vault/NTM_RU_AB/Transport_Network_Neumann
/opt/matlab/r2022b/bin/matlab -nodisplay -nosplash -batch "call_NTM_brandontest_sweep_cluster"