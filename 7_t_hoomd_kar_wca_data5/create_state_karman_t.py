import sys
import os
sys.path.append('/home/isobelab2022/build3/hoomd')
import itertools
import math

import gsd.hoomd
import hoomd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
from numba import jit, f8, i8, b1, void,njit
import matplotlib.patches as pat
import sys
import os
import fresnel
from PIL import Image
import IPython
import packaging.version
import random
import time
 

first_time=time.time()
def initialize_ra(Nx,Ny,dx,dy,lx,ly,remove_dx_num):
    rx=[]
    ry=[]
#     rx.append(lx/4.0-lx/2)
#     ry.append(ly/2.0-ly/2)
    
    for i in range(Nx):
        for j in range(Ny):
            if (j>=int(Ny/2)-remove_dx_num and j<=int(Ny/2)+remove_dx_num):
                if(i>=int(Nx/4)-remove_dx_num and i<=int(Nx/4) +remove_dx_num):
                    continue


            # rx.append(i*dx+dx/2)
            # ry.append(j*dy+dy/2)
            
            rx.append(i*dx-lx/2+dx/2)
            ry.append(j*dy-ly/2+dy/2)



    return rx,ry

class RelativeFlow(hoomd.custom.Action):
    def __init__(self, ave_flow,h):
        self.ave_flow=ave_flow
        self.h=h
        
    def act (self,timestep):
        snap=self._state.get_snapshot()
        if snap.communicator.rank==0:
            snap.particles.position[:,0]+=self.ave_flow*self.h
        self._state.set_snapshot(snap)
# class ChangeVelocity(hoomd.custom.Action):
#     def __init__(self, velocity,lx):
#         self.velocity=velocity
#         self.lx=lx
        
#     def act (self,timestep):
#         snap=self._state.get_snapshot()
#         if snap.communicator.rank==0:
#             position=snap.particles.position
#             position_x=position.T[0]
#             # change_index=np.where(position_x>lx-1)
#             # snap.particles.velocity[change_index,0]=self.velocity
#             snap.particles.velocity[:,0]+=self.velocity
#         self._state.set_snapshot(snap)
rho=float(sys.argv[1])
ave_flow=float(sys.argv[2])
static_dia=float(sys.argv[3])
reduced_speed=float(sys.argv[4])
rotational_diffusion=float(sys.argv[5])

ver=str(rho)+"_"+str(ave_flow)+"_"+str(static_dia)+"_"+str(reduced_speed)+"_"+str(rotational_diffusion)

init_ver=str(rho)+"_"+str(0.0)+"_"+str(static_dia)+"_"+str(0.0)+"_"+str(5.0)

main_dir="./"+ver


input_traj = gsd.hoomd.open('../7_hoomd_kar_wca_data5/{0}/log_pos_{0}.gsd'.format(ver), 'rb')

input_pos=input_traj[-1].particles.position.T
rx=input_pos[0]
ry=input_pos[1]


dx=np.sqrt(1.0/rho)
dy=dx
lx=input_traj[0].configuration.box[0]
ly=input_traj[0].configuration.box[1]
Nx=int(lx/dx)
Ny=int(ly/dy)


figsize=(10.0,10.0*ly/lx)
remove_dx_num=math.ceil(static_dia*0.5/dx)


# temp=1.0
m=1.0
epsilon=1.0


# 刻み幅小さすぎの可能性もあるから大きめにしてみてもいいかも

# h=5e-4
# dt = 1e-6
real_time=80
dt = 1e-4
nsteps=int(real_time/dt)
pos_hout=int(nsteps/500)
thermo_hout=int(nsteps/2.0)

# t_end=nsteps*h
# h_rev=1/h 
# h2=0.5*h*h



# rx,ry=initialize_ra(Nx,Ny,dx,dy,lx,ly,remove_dx_num)
# rx=np.loadtxt("../data_file/square/rx_{0}_{1}.txt".format(rho,static_dia))
# ry=np.loadtxt("../data_file/square/ry_{0}_{1}.txt".format(rho,static_dia))

N=len(rx)
print(N)
position=[]
for i in range(len(rx)):
    position.append((rx[i],ry[i],0.0))



snapshot = gsd.hoomd.Snapshot()
snapshot.particles.N = N
snapshot.particles.position = position[0:N]


# (1,0,0,0,)
# snapshot.particles.orientation = orientation
snapshot.particles.typeid = (N)*[0]
snapshot.particles.types = ['Move']
snapshot.particles.diameter=(N)*[1.0]
snapshot.particles.mass = np.ones((N))
snapshot.configuration.box = [lx, ly, 0, 0, 0, 0]


# with gsd.hoomd.open(name='_init'+str(static_dia)+'.gsd', mode='xb') as f:
#     f.append(snapshot)





print(main_dir)
output_dir=main_dir+"/figure_2d"
if not os.path.exists(output_dir): os.makedirs(output_dir)
os.chdir(main_dir)

sim = hoomd.Simulation(device=hoomd.device.GPU(), seed=12)


# 結局ここで最後のsnap挿入
last_snap=input_traj[-1]
sim.create_state_from_snapshot(last_snap)
# Integration information


sigma_aa=1.0
sigma_ab=(static_dia/2+0.5)
epsilon_aa=1.0
epsilon_ab=1.0

cell = hoomd.md.nlist.Cell(buffer=0.4)
lj = hoomd.md.pair.LJ(nlist=cell, mode="shift")
lj.params[("Move", "Move")] = dict(epsilon=epsilon_aa, sigma=sigma_aa)
lj.r_cut[("Move", "Move")] = 2**(1/6)*sigma_aa


walls=[hoomd.wall.Cylinder(origin=(lx/4-lx/2,ly/2-ly/2,0),radius=static_dia/2,inside=False,axis=(0,0,1))]
ljw = hoomd.md.external.wall.LJ(walls=walls)
# ljw.params['Move'] = {"epsilon": epsilon_ab, "sigma": sigma_ab, "r_cut": sigma_ab*2**(1/6)}

ljw.params['Move'] = {"epsilon": 0.1, "sigma": 1.0, "r_cut": sigma_ab*2**(1/6)-static_dia/2}

ktemp = 0.5

# # rotational_diffusion = 20.0
# # Apply Stokes-Einstein
# traslational_diffusion = 3.0 * rotational_diffusion
# #???????
# brownian = hoomd.md.methods.Brownian(filter=hoomd.filter.All(), kT=ktemp)
# brownian.gamma.default = ktemp / traslational_diffusion
# brownian.gamma_r.default =[ktemp / rotational_diffusion,ktemp / rotational_diffusion,0]
# # brownian.gamma_r.default = np.full((3,), ktemp / rotational_diffusion)
# active = hoomd.md.force.Active(filter=hoomd.filter.All())

# act_force = reduced_speed * brownian.gamma.default
# active.active_force["Move"] = (act_force, act_force,0)
# active.active_torque["Move"] = (0, 0, 0)


# flow_force=ave_flow* brownian.gamma.default
# constant = hoomd.md.force.Constant(
#     filter=hoomd.filter.All()
#     )
# constant.constant_force['Move'] = (flow_force,0,0)

# rotational_diffusion_updater = active.create_diffusion_updater(
#  trigger=hoomd.trigger.Periodic(1), rotational_diffusion=rotational_diffusion
# )

nve = hoomd.md.methods.NVE(filter=hoomd.filter.All())
integrator = hoomd.md.Integrator(
 dt=dt,
 methods=[nve],
 forces=[lj, ljw],
)

velocity_operation=hoomd.update.CustomUpdater(action=RelativeFlow(ave_flow,dt),trigger=1)
sim.operations+=velocity_operation
# sim.operations += rotational_diffusion_updater
sim.operations.integrator = integrator

# カスタムクラス
class PrintTimestep(hoomd.custom.Action):
    def act (self,timestep):
        print(timestep)
custom_action = PrintTimestep()
custom_op = hoomd.write.CustomWriter(action=custom_action,
                                 trigger=hoomd.trigger.Periodic(10000))
sim.operations.writers.append(custom_op)






# logger定義
pos_logger = hoomd.logging.Logger()
pos_logger.add(sim, quantities=['timestep', 'walltime'])
pos_logger.add(lj ,quantities=["forces"])
# pos_logger.add(thermodynamic_properties,quantities=["kinetic_temperature","kinetic_energy","potential_energy","volume"])
gsd_writer_pos = hoomd.write.GSD(filename="log_pos_"+ver+".gsd",
                             trigger=hoomd.trigger.Periodic(pos_hout),
                             mode='xb',
                             filter=hoomd.filter.All())
gsd_writer_pos.log = pos_logger

sim.operations.writers.append(gsd_writer_pos)

##########################thermo_log#######################
thermodynamic_properties = hoomd.md.compute.ThermodynamicQuantities(
    filter=hoomd.filter.All())
sim.operations.computes.append(thermodynamic_properties)
thermo_logger= hoomd.logging.Logger()
thermo_logger.add(thermodynamic_properties)
gsd_writer_thermo = hoomd.write.GSD(filename="log_thermo_"+ver+".gsd",
                             trigger=hoomd.trigger.Periodic(thermo_hout),
                             mode='xb',
                             filter=hoomd.filter.Null())
gsd_writer_thermo.log = thermo_logger

sim.operations.writers.append(gsd_writer_thermo)

########################################################

# sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=ktemp)
# 初期速度変化
# sim.run(0)

print("-----run--------")
second_time=time.time()

sim.run(nsteps)

print(time.time()-first_time)

print(time.time()-second_time)
os.chdir("../")
traj = gsd.hoomd.open('./'+ver+'/log_pos_'+ver+'.gsd', 'rb')


# traj = gsd.hoomd.open('log_force2d_'+ver+'.gsd', 'rb')
plt.figure(figsize=figsize)

print(len(traj))

for t in range(len(traj)-1,0,-2):
    print(t)
    bx=plt.axes()
    plt.axis([-lx/2,lx/2,-ly/2,ly/2])
    position=traj[t].particles.position   
    c=pat.Circle(xy=(lx/4-lx/2,ly/2-ly/2),radius=static_dia/2,fc="b")
    bx.add_patch(c) 
    for i in range(N):
        c=pat.Circle(xy=(position[i][0],position[i][1]),radius=0.5,fc="r")
        bx.add_patch(c)

    plt.title("step"+str(t))
    plt.savefig(output_dir+"/figure{0}.png".format(t))
    plt.cla()

    ############アニメーション################    
images=[]
# image_num=sum(os.path.isfile(os.path.join(pic_output name)) for name in os.listdir(pic_output))
image_num=sum(os.path.isfile(os.path.join(output_dir,name))for name in os.listdir(output_dir))
print(image_num)
for i in range(0,image_num):
    file_name=output_dir+"/figure"+str(i)+".png"
    im=Image.open(file_name)
    images.append(im)

gif_output_dir=main_dir+"/abpgif2"

if not os.path.exists(gif_output_dir): os.makedirs(gif_output_dir)
images[0].save(gif_output_dir+"/out_ela2.gif",save_all=True,append_images=images[1:],loop=0,duration=10)
    
    