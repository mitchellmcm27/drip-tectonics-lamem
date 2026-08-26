#!/bin/bash

#SBATCH -J shortening_lamem
#SBATCH --account=gts-ssim33-atlas
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --time=12:00:00
#SBATCH --output=%A\_%a.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mitchellmcm27@gmail.com
#SBATCH --array=0-80

module purge
module load gcc/12 openmpi/4.1.5
module load julia/1.11.3

cd $SLURM_SUBMIT_DIR

# create list of model versions
source scripts/versions_list_20260305.sh
local versions=()
create_versions versions

# optionally override defaults by passing arguments to the slurm script
# append arguments to the end of each version string, so they override the defaults
for arg; do
    for i in "${!versions[@]}"; do
        versions[$i]="${versions[$i]} $arg"
    done
done

# NOTE: Make sure job array matches total number

printf ' - %q\n' "${versions[@]}"

# get the specific version for this job
ver=${versions[$SLURM_ARRAY_TASK_ID]}
date_slug=$(date +"%Y%m%d")
# parse the version string and return a consistent slug (Julia script)
ver_slug=$(./scripts/version_slug "$ver")

out_dir=output/pace/${date_slug}_${ver_slug}

# generate input files
julia 2d_plateau.jl -n $SLURM_NTASKS --version "$ver" --out-dir "$out_dir" --dry-run --yse --aps

# run the model
cd $out_dir
srun $HOME/shared/LaMEM/bin/opt/LaMEM -ParamFile fixed.dat
