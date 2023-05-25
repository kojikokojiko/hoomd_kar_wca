#!/bin/bash


# cd /home/isobelab2022/projects/hoomd_kar_wca1

source ~/anaconda3/etc/profile.d/conda.sh
conda activate iwase_env1


rho=0.8
ave_flow=8.0
static_dia=30.0 
reduced_speed=0.0 
rotational_diffusion=0.0

python create_state_karman.py $rho $ave_flow $static_dia $reduced_speed  $rotational_diffusion