#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate iwase_env1

rho_list=(0.8)
ave_flow_list=( 5.0 8.0 10.0 15.0 20.0 )
static_dia=40.0
reduced_speed_list=(0.0 )
rotational_diffusion_list=(0.0)
for rho in "${rho_list[@]}"
do 
    for ave_flow in  "${ave_flow_list[@]}"
    do
        for reduced_speed in  "${reduced_speed_list[@]}"
        do

            for rotational_diffusion in  "${rotational_diffusion_list[@]}"
            do
                echo $rho $ave_flow $static_dia $reduced_speed  $rotational_diffusion
                name="${rho}_${ave_flow}_${static_dia}_${reduced_speed}_${rotational_diffusion}"
                nohup time python create_state_karman.py $rho $ave_flow $static_dia $reduced_speed  $rotational_diffusion > ${name}.log 2>&1 &
                # nohup time python anim.py $rho $ave_flow $static_dia $reduced_speed  $rotational_diffusion > ${name}.log 2>&1 &
                echo $!>${name}.txt
            done
            
        done
    done
done

# g++ -std=c++11 -o iabp2 iABP2.cpp
# nohup time ./iabp2 1.0 1.0 10000 # Pe M N
conda activate iwase_env1
python create_state_karman.py 0.8 





# rho=float(sys.argv[1])
# ave_flow=float(sys.argv[2])
# static_dia=float(sys.argv[3])
# reduced_speed=float(sys.argv[5])
# rotational_diffusion=float(sys.argv[6])

