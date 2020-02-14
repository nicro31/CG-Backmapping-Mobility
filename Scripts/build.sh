#!/usr/bin/env bash
#SBATCH -U "snic2017-12-59" -N 1 -n 32
#SBATCH -t  05:00:00


LOE-CTP-FRAG build AA.pdb config.txt AA.grp

