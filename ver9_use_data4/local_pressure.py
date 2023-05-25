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


def make_neighbor_list(around_grid_num):
    if around_grid_num==5:
        neighbor_col=[
            -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
            -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
            -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
            -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
            -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
            -5,-4,-3,-2,-1,    1, 2, 3, 4, 5,
            -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
            -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
            -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
            -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
            -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5
        ]
        neighbor_row=[
            -5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,
            -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
            -3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,
            -2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
            -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
             0, 0, 0, 0, 0,    0, 0, 0, 0, 0,
             1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
             2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
             3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
             4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
             5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5

        ]
    return neighbor_row,neighbor_col





rho=float(sys.argv[1])
ave_flow=float(sys.argv[2])
static_dia=float(sys.argv[3])
reduced_speed=float(sys.argv[4])
rotational_diffusion=float(sys.argv[5])
ver=str(rho)+"_"+str(ave_flow)+"_"+str(static_dia)+"_"+str(reduced_speed)+"_"+str(rotational_diffusion)


main_dir="./"+ver


# temp_dir="./"+ver+"/log_pos_"+ver+".gsd"
temp_dir="./"+ver+"/log_pos_"+ver+".gsd"
local_p_dir=main_dir+"/local_p"
# i_phi_dir=main_dir+"/i_phi"
if not os.path.exists(local_p_dir): os.makedirs(local_p_dir)
# if not os.path.exists(i_phi_dir): os.makedirs(i_phi_dir)


print(temp_dir)
traj = gsd.hoomd.open(temp_dir, 'rb')
# traj = gsd.hoomd.open(dir, 'rb')




lx=traj[0].configuration.box[0]
ly=traj[0].configuration.box[1]
figsize=(10,10*ly/lx)
plt.figure(figsize=figsize)

ppi=72 # points per inche 

print(lx)
NP=len(traj[0].particles.position)
print(NP)


# 検査体積 binn
l_grid_c=1.5
n_grid_c=math.ceil(lx/l_grid_c)
Vc=l_grid_c*ly
bin_x=[i*l_grid_c for i in range(n_grid_c)]
print(n_grid_c)
##################



sigma=1.0
# あんま頭よくない処理
size=sigma
sigma2=sigma*sigma
sigma6=sigma2*sigma2*sigma2
sigma12=sigma6*sigma6
rc=2**(1.0/6.0)*sigma
rc2=rc*rc
m=1.0

epsilon=1.0
around_grid_num=5
neighbor_row,neighbor_col= make_neighbor_list(around_grid_num)
neighbor_len=len(neighbor_row)


bo=(lx/0.4)
l_gx_d=lx/bo
l_gx=l_gx_d
l_gy=l_gx
n_gx=math.ceil(lx/l_gx)
n_gy=math.ceil(ly/l_gy)
n_g_len=n_gx*n_gy

r_cut_sann=l_gx_d*around_grid_num
pair_length_g=50
NN_MAX=20
Nmax=80


# 描画設定##############################################
# https://qiita.com/stanaka2/items/c40841f858d7083aad4e
fig = plt.figure(figsize=figsize, dpi=100.0)
# ax  = fig.add_axes((0.1,  0.1, 0.8, 0.8))

# # 枠の範囲指定
# xmin,xmax=-lx/2.0, lx/2.0
# ymin,ymax=-ly/2.0, ly/2.0

# # 描写範囲の長さを取得(dpi単位)
# # x軸をベースに計算しているがy軸でも同じ。アスペクト比が違えばおかしくなる
# ax_length=ax.bbox.get_points()[1][1]-ax.bbox.get_points()[0][1]

# # dpi単位の長さをポイント単位に変換(dpiからインチ単位にし、インチからポイント単位に変換)
# ax_point = ax_length*ppi/fig.dpi

# # x軸の実スケールをポイント単位に変換
# xsize=xmax-xmin
# fact=ax_point/xsize

# # scatterのマーカーサイズは直径のポイントの二乗を描くため、実スケールの半径をポイントに変換し直径にしておく
# size*=2*fact

####################################################
for t in range(len(traj)-1,0,-10):
    
    print(t)
    pos=traj[t].particles.position.T
    rx=pos[0]
    ry=pos[1]

    vel=traj[t].particles.velocity.T
    vx=vel[0]
    vy=vel[1]


    ax=np.zeros(NP)
    ay=np.zeros(NP)

    w_list_c=np.zeros(n_grid_c)

    # Gmap_create########################
    g_map=np.full((n_gy,n_gx),-1,dtype="int")
    pair_list_g=np.full((NP,pair_length_g),-1,dtype="int")

    for i in range(NP):
        gx_map=int((rx[i]+lx/2)/l_gx)
        gy_map=int((ry[i]+ly/2)/l_gy)
        if(g_map[gy_map][gx_map]!=-1):
            print("TOO LARGE GRID")
            sys.exit()
        g_map[gy_map][gx_map]=i
    # np.savetxt(main_dir+"/g_map{0}.txt".format(t),g_map)

    ##################################################################
    # make_pairlist_g####################################

    for i in range(n_gy):
        for j in range(n_gx):
            select_index=g_map[i][j]
            if select_index==-1:
                continue
            particle_counter=0
            for k in range(neighbor_len):
                search_gx=j+neighbor_col[k]
                search_gy=i+neighbor_row[k]
                if(search_gx>=n_gx):
                    search_gx-=n_gx
                elif(search_gx<0):
                    search_gx+=n_gx
                if (search_gy>=n_gy):
                    search_gy-=n_gy
                elif(search_gy<0):
                    search_gy+=n_gy
                
                search_index=g_map[search_gy][search_gx]
                if search_index==select_index:
                    print("what wrong")
                    print(select_index)
                if search_index!=-1:
                    pair_list_g[select_index][particle_counter]=search_index
                    particle_counter+=1

            if particle_counter==pair_length_g-1:
                print("ERROR HAPPEN :PAIR_LENGH_G is few")
                sys.exit()
                
            pair_list_g[select_index][-1]=particle_counter
    ############################################

    # 力の計算############################
    for i in range(NP):
        roop_num=pair_list_g[i][-1]
        for j in range(roop_num):

            gx_c=int((rx[i]+lx/2.0)/l_grid_c)
     
            pair_index=pair_list_g[i][j]
            rxij=rx[i]-rx[pair_index]
            # minimum image convention
            if (rxij>=lx/2):
                rxij=rxij-lx
            elif (rxij<-lx/2):
                rxij=rxij+lx
            else:
                rxij=rxij
            
            ryij=ry[i]-ry[pair_index]
            # minimum image convention
            if (ryij>=ly/2):
                ryij=ryij-ly
            elif (ryij<-ly/2):
                ryij=ryij+ly
            else:
                ryij=ryij

            r2=rxij*rxij+ryij*ryij
            ir2=1.0/r2
            ir6=ir2*ir2*ir2
            ir12=ir6*ir6
            if r2>rc2:
                fx=0.0
                fy=0.0
            else:
                fx=24.0*epsilon*(2.0*sigma12*ir12-sigma6*ir6)*ir2*rxij
                fy=24.0*epsilon*(2.0*sigma12*ir12-sigma6*ir6)*ir2*ryij

            # ax[i]+=fx
            # ay[i]+=fy
            # ax[pair_index]-=fx
            # ay[pair_index]-=fy
            w_list_c[gx_c]+=rxij*fx+ryij*fy
            # w_list_c[pair_index]+=rxij*fx+ryij*fy

            
    #############################運動エネルギーを計算
    K_list_c=np.zeros(n_grid_c)
    # print(len(K_list_c))
    for i in range(NP):
        gx_c=int((rx[i]+lx/2.0)/l_grid_c)
        K_list_c[gx_c]+=m*(vx[i]*vx[i]+vy[i]*vy[i])

    p_list_c=np.zeros(n_grid_c)

    for i in range(n_grid_c):
        p_list_c[i]=(K_list_c[i]+0.5*w_list_c[i])/(3*Vc)
    



# 描画############################
    plt.scatter(bin_x,p_list_c)
    plt.savefig(local_p_dir+"/figure{0}.png".format(t))
    plt.cla()
    plt.clf()
    # plt.close()







# Real_phi6#################################



# dpiは任意の値を設定できる。デフォルトは100
    # fig = plt.figure(figsize=(12.5,5.0), dpi=100.0)
    # ax  = fig.add_axes((0.1,  0.1, 0.8, 0.8))

    # ax.set_aspect("equal") # アスペクト比を等しくする
    #                     # 今は枠のサイズ設定で最初から等しいのであってもなくてもよい

    # # 二乗にして与える 
    # im=ax.scatter(rx,ry,s=size**2,c=phi6_2,cmap="Reds", linewidths=0)
    # # for i in range(0,NP,50):
    # #     c=pat.Circle(xy=(rx[i],ry[i]),radius=r_sann_1[i],fc="g",alpha=0.2)
    # #     ax.add_patch(c)

    # # for i in range(100,200):
    # #     ax.text(rx[i], ry[i], i,fontsize="small")

    # c=pat.Circle(xy=(lx/4-lx/2,ly/2-ly/2),radius=static_dia/2,fc="b")
    # ax.add_patch(c) 
    # plt.colorbar(im)

    # ax.set_xlim(xmin,xmax)
    # ax.set_ylim(ymin,ymax)

    # ax.set_xticks(np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 5))
    # ax.set_yticks(np.linspace(ax.get_ylim()[0], ax.get_ylim()[1], 5))

    # # ax.grid(which='both', axis='both')

    # # plt.colorbar()
    # # plt.colorbar(sc, ax=ax)
    # plt.title("phi6_2")
    # plt.savefig(phi6_2_dir+"/phi6_2_{0}.png".format(t))
    # # plt.savefig(main_dir+"/figure.png")
    # # plt.cla()
    # plt.close()
