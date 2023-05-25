import sys
import os
sys.path.append('/home/isobelab2022/build3/hoomd')

import matplotlib.pyplot as plt
import matplotlib.patches as pat
import gsd.hoomd
import hoomd
import sys
from PIL import Image
import os
import math
import  numpy as np
import seaborn as sns
import sys
import cv2



rho=float(sys.argv[1])
ave_flow=float(sys.argv[2])
static_dia=float(sys.argv[3])
reduced_speed=float(sys.argv[4])
rotational_diffusion=float(sys.argv[5])
ver=str(rho)+"_"+str(ave_flow)+"_"+str(static_dia)+"_"+str(reduced_speed)+"_"+str(rotational_diffusion)


main_dir="./"+ver


temp_dir="./"+ver+"/log_pos_"+ver+".gsd"
rho_dir=main_dir+"/rho_map"
# i_phi_dir=main_dir+"/i_phi"
if not os.path.exists(rho_dir): os.makedirs(rho_dir)
# if not os.path.exists(i_phi_dir): os.makedirs(i_phi_dir)


print(temp_dir)
traj = gsd.hoomd.open(temp_dir, 'rb')
# traj = gsd.hoomd.open(dir, 'rb')


# plt.figure(figsize=(12.5,5.0))

# print(traj[0].particles.position)
# print(len(traj))
lx=traj[0].configuration.box[0]
ly=traj[0].configuration.box[1]

figsize=(10,10*ly/lx)
NP=len(traj[0].particles.position)
print(NP)

sigma=1.0
r=0.5
# あんま頭よくない処理
size=sigma


bo=(lx/2.0)
l_gx_d=lx/bo
l_gx=l_gx_d
l_gy=l_gx
n_gx=math.ceil(lx/l_gx)
n_gy=math.ceil(ly/l_gy)
n_g_len=n_gx*n_gy




for t in range(len(traj)-1,0,-10):
    print(t)
    pos=traj[t].particles.position.T
    rx=pos[0]
    ry=pos[1]
    rho_map=np.zeros((n_gy,n_gx))
    rho_map2=np.zeros((n_gy,n_gx))
    for i in  range(NP):
        gx_map=int((rx[i]+lx/2)/l_gx)
        gy_map=int((ry[i]+ly/2)/l_gy)
        rho_map[gy_map][gx_map]+=1

    for i in range(n_gy):
        for j in range(n_gx):
            rho_map2[i][j]=(rho_map[i][j])/(l_gx_d*l_gx_d)

    for i in range(n_gy):
        for j in range(n_gx):
            rho_map[i][j]=rho_map2[i][j]


    D_MAP = cv2.GaussianBlur(rho_map, (5, 5), sigmaX=1)

    plt.figure(figsize=figsize)
    plt.xticks([])
    plt.yticks([])
    ax = plt.gca()

    im = ax.imshow(D_MAP, interpolation = 'nearest', cmap = "Reds")

    plt.colorbar(im,  cmap = "Reds")
    plt.savefig(main_dir+"/rho_map/local_dense{0}.png".format(t))

    plt.close()