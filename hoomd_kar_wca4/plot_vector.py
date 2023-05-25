import matplotlib.pyplot as plt
import gsd.hoomd
import hoomd
import sys
from PIL import Image
import os
rho=float(sys.argv[1])
ave_flow=float(sys.argv[2])
static_dia=float(sys.argv[3])
reduced_speed=float(sys.argv[4])
rotational_diffusion=float(sys.argv[5])
ver=str(rho)+"_"+str(ave_flow)+"_"+str(static_dia)+"_"+str(rotational_diffusion)+"_"+str(reduced_speed)

main_dir="./"+ver


traj = gsd.hoomd.open('./'+ver+'/log_force2d_'+ver+'.gsd', 'rb')


plt.figure(figsize=(12.5,5.0))

print(traj[0].particles.position)
# print(len(traj))
lx=traj[0].configuration.box[0]
ly=traj[0].configuration.box[1]
# print(lx)
# print(ly)
output_dir=main_dir+"/v_vec"

if not os.path.exists(output_dir): os.makedirs(output_dir)
for t in range(len(traj)-1):
    print(t)
    bx=plt.axes()
    plt.axis([-lx/2,lx/2,-ly/2,ly/2])
    pos=traj[t].particles.position.T
    vel=traj[t].particles.velocity.T



    plt.quiver(pos[0],pos[1],vel[0],vel[1])
    # c=pat.Circle(xy=(lx/4-lx/2,ly/2-ly/2),radius=static_dia/2,fc="b")
    # bx.add_patch(c) 
    # for i in range(N):
    #     c=pat.Circle(xy=(position[i][0],position[i][1]),radius=0.5,fc="r")
    #     bx.add_patch(c)

    plt.title("step"+str(t))
    plt.savefig(output_dir+"/v_vec{0}.png".format(t))
    plt.cla()



    ############アニメーション################    
images=[]
# image_num=sum(os.path.isfile(os.path.join(pic_output name)) for name in os.listdir(pic_output))
image_num=sum(os.path.isfile(os.path.join(output_dir,name))for name in os.listdir(output_dir))
print(image_num)
for i in range(0,image_num):
    file_name=output_dir+"/v_vec"+str(i)+".png"
    im=Image.open(file_name)
    images.append(im)

gif_output_dir=output_dir+"/gif"

if not os.path.exists(gif_output_dir): os.makedirs(gif_output_dir)
images[0].save(gif_output_dir+"/v_vec.gif",save_all=True,append_images=images[1:],loop=0,duration=10)
    
    
# for data in traj:
#     v=data[0].particles.velocity
#     r=data[0].particles.position



# print(v0)