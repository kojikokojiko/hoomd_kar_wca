
import sys
import os
sys.path.append('/home/isobelab2022/build3/hoomd')
import matplotlib.pyplot as plt
import gsd.hoomd
import hoomd
import sys
from PIL import Image
import os
import math
import  numpy as np
rho=float(sys.argv[1])
ave_flow=float(sys.argv[2])
static_dia=float(sys.argv[3])
reduced_speed=float(sys.argv[4])
rotational_diffusion=float(sys.argv[5])
ver=str(rho)+"_"+str(ave_flow)+"_"+str(static_dia)+"_"+str(rotational_diffusion)+"_"+str(reduced_speed)


main_dir="./"+ver


traj = gsd.hoomd.open('./'+ver+'/log_pos_'+ver+'.gsd', 'rb')


print(traj[0].particles.position)
# print(len(traj))
lx=traj[0].configuration.box[0]
ly=traj[0].configuration.box[1]
print(lx)
print(ly)
output_dir=main_dir+"/vor_cg"
cg_l_g=5.0
cg_ngx=math.ceil(lx/cg_l_g)
cg_ngy=math.ceil(ly/cg_l_g)

vor_map=np.zeros((cg_ngy,cg_ngx))

plt.figure(figsize=(10,10*ly/lx))


def scale(arr):
    arr=np.array(arr)
    return  (arr-np.min(arr))/(np.max(arr)-np.min(arr))*2-1
def cg(pos,vel,cg_l_g,cg_ngx,cg_ngy,ave_flow):
    p_count_map=np.zeros((cg_ngy,cg_ngx))
    vx_map=np.zeros((cg_ngy,cg_ngx))
    vy_map=np.zeros((cg_ngy,cg_ngx))
    rx_map=np.zeros((cg_ngy,cg_ngx))
    ry_map=np.zeros((cg_ngy,cg_ngx))
    Rx=pos[0]
    Ry=pos[1]
    Vx=vel[0]
    Vy=vel[1]

    for i in range(len(Rx)):
        cg_gx_index=int((Rx[i]+lx/2)/cg_l_g)
        cg_gy_index=int((Ry[i]+ly/2)/cg_l_g)
        p_count_map[cg_gy_index][cg_gx_index]+=1
        vx_map[cg_gy_index][cg_gx_index]+=(Vx[i]+ave_flow)
        vy_map[cg_gy_index][cg_gx_index]+=Vy[i]
    
    for i in range(cg_ngy):
        for j in range(cg_ngx):
            rx_map[i][j]=cg_l_g/2+cg_l_g*j-lx/2
            ry_map[i][j]=cg_l_g/2+cg_l_g*i-ly/2
            
            p_count=p_count_map[i][j]
            if (p_count==0):
                p_count=1
            vx_map[i][j]/=p_count
            vy_map[i][j]/=p_count

    # rx_list=np.ravel(rx_map)
    # ry_list=np.ravel(ry_map)
    # vx_list=np.ravel(vx_map)
    # vy_list=np.ravel(vy_map)

    return rx_map,ry_map,vx_map,vy_map




if not os.path.exists(output_dir): os.makedirs(output_dir)
print(len(traj))
for t in range(len(traj)-1,0 ,-4):
    print(t)
    bx=plt.axes()
    # plt.axis([-lx/2,lx/2,-ly/2,ly/2])
    pos=traj[t].particles.position.T
    vel=traj[t].particles.velocity.T
    
    rx_map,ry_map,vx_map,vy_map=cg(pos,vel,cg_l_g,cg_ngx,cg_ngy,ave_flow)
    vor_map=np.zeros((cg_ngy,cg_ngx))


    # 速度場のx,y方向の偏微分を計算する
    dvy_dx = np.gradient(vy_map, axis=1)
    dvx_dy = np.gradient(vx_map, axis=0)

    # 渦度を計算する
    vor_map = dvy_dx - dvx_dy


    # for i in range(1,cg_ngy-1):
    #     for j in range(1,cg_ngx-1):
    #         vor_map[i][j]=(vy_map[i+1][j]-vy_map[i-1][j]-vx_map[i][j+1]+vx_map[i][j-1])/2


    vor_map=scale(vor_map)
    

    # for i in range(len(vor_map)):
    #     for j in range(len(vor_map[i])):
    #         if vor_map[i][j] < -1 or vor_map[i][j] > 1:
    #             raise 





    # plt.quiver(rx_list,ry_list,vx_list,vy_list)
    # c=pat.Circle(xy=(lx/4-lx/2,ly/2-ly/2),radius=static_dia/2,fc="b")
    # bx.add_patch(c) 
    # for i in range(N):
    #     c=pat.Circle(xy=(position[i][0],position[i][1]),radius=0.5,fc="r")
    #     bx.add_patch(c)
    plt.imshow(vor_map,cmap="seismic",vmin=-1,vmax=1)
    plt.colorbar()

    plt.title("step"+str(t))
    plt.savefig(output_dir+"/v_vor_cg{0}.png".format(t))
    plt.cla()
    plt.clf()



    ############アニメーション################    
# images=[]
# # image_num=sum(os.path.isfile(os.path.join(pic_output name)) for name in os.listdir(pic_output))
# image_num=sum(os.path.isfile(os.path.join(output_dir,name))for name in os.listdir(output_dir))
# print(image_num)
# for i in range(0,image_num):
#     file_name=output_dir+"/v_vec_cg"+str(i)+".png"
#     im=Image.open(file_name)
#     images.append(im)

# gif_output_dir=output_dir+"/gif"

# if not os.path.exists(gif_output_dir): os.makedirs(gif_output_dir)
# images[0].save(gif_output_dir+"/v_vec_cg.gif",save_all=True,append_images=images[1:],loop=0,duration=10)
    
    
# for data in traj:
#     v=data[0].particles.velocity
#     r=data[0].particles.position



# print(v0)