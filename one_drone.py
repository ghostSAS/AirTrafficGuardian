from utils.pkg_func import *
from utils.Bezier import Bezier
from utils.Drone import Drone
from utils.Planner import Planner

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
corridor = [np.array([[1,1,1],[1,1,2.5]]),
            np.array([[2.5,2.5,3.5],[2.5,3.5,3.5]])]
corridor_r = 0.15

degree = 10
dim = 3
T = 4
dt = .08
m = int(T/dt)+1


# t_sample = np.linspace(0,T,50)[[0,20,24,28,32,-1]]
# corridor_sample = np.linspace(corridor[0], corridor[1], len(t_sample)-2)

corridor_info = {}
corridor_info['end_points'] = corridor
corridor_info['radius'] = corridor_r

drone = Drone(start, target, corridor_info, T, 0)
drone.ctrlPt.set_degree_dim(degree, dim)

ts = time.time()
drone.get_primary_traj()
print(f"QP takes {time.time()-ts:.4f} sec")

planner = Planner()
starts = [start]
targets = [target]
corridors = [corridor]
trajs = [drone.traj_pt]

view_angle=[-31, 34]
planner.plot_final_frame(trajs, starts, targets, corridors, corridor_r, view_angle)

plt.show()



