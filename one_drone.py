from utils.pkg_func import *
from utils.Trajectory import *
from utils.Drone import *
from utils.Planner import Planner
import utils.toolbox_visualize as vis

import casadi as ca


"""
task: 
    drone 1 [0,0,0] --> [5,5,5]
    corridor [2,2,2] --> [2,2,3.5], r=0.15
    
goal:
    1. reach the target
    2. through the corridor with radius 0.5m
    3. minimum snap
"""


"""
bezier curve B:
    degree = 10
    T = 4
    dt = .08
    m = int(T/dt)+1=51
"""

# ---------------- configuration -------------
start = np.array([0,0,0])
target = np.array([5,5,5])


corridor = [np.linspace([2.5,2.5,3.5],[2.5,3.5,3.5], 2)]

corridor = [np.linspace([1,1,1],[1,1,2.5],3),
            np.linspace([2.5,2.5,3.5],[2.5,3.5,3.5],3)]

corridor_r = 0.15  

corridor_info = {}
corridor_info['end_points'] = corridor
corridor_info['radius'] = corridor_r

# ------------- parameters ---------------
optimal_ord_idx = [3, 4]
kT=70
optiT = True

order = 5
dim = 3
T = 6


# -------------------------------- old methods with single segements -----------------

# drone = Drone(start, target, corridor_info, T, 0)
# drone.ctrlPt.set_degree_dim(order, dim)

# ts = time.time()
# drone.get_primary_traj()
# print(f"QP takes {time.time()-ts:.4f} sec")

# planner = Planner(corridor_r)
# starts = [start]
# targets = [target]
# corridors = [corridor]
# drones = [drone]

# view_angle=[-31, 34]
# planner.plot_final_frame(drones, corridors, view_angle)

# plt.show()

# -------------------------------- new methods with multiple segements -----------------

planner = Planner(corridor_r)

drone = Drone_traj(start, target, corridor_info, order, T, priority=0, idx_w=optimal_ord_idx)

# get the primary trajectory results
T_orig = [sp.T for sp in drone.traj_bezier.splines]

if optiT:
    planner.get_primary_traj_optiT(drone, kT=kT)
else:
    planner.get_primary_traj(drone)

# compare time for each segments
T_now = [sp.T for sp in drone.traj_bezier.splines]
formatted_T_orig = ', '.join([f"{value:.2f}" for value in T_orig])
formatted_now_values = ', '.join([f"{value:.2f}" for value in T_now])
print(f"Original T: {formatted_T_orig}, total: {sum(T_orig):.2f} ")
print(f"Now T: {formatted_now_values}, total: {sum(T_now):.2f} ")

# starts = [start]
# targets = [target]
# corridors = [corridor]
drones = [drone]


view_angle=[-31, 34]

traj_pt, t_span = drone.get_traj_pt(derivative=3)
vis.visualize_splines_1D(traj_pt, len(start), t_span)
vis.plot_final_frame_3D(drones, corridor_info, view_angle, verbose={'show':True, 'save':False})
# vis.plot_animation_3D(drones, corridor_info, view_angle, verbose={'show':True, 'save':False})



