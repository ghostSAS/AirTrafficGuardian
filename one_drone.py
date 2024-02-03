from utils.pkg_func import *
from utils.Trajectory import *
from utils.Drone import *
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
start = np.array([0.1,0.1,0.1])
target = np.array([5,5,5])
# corridor = [np.array([[1,1,1],[1,1,2.5]]),
#             np.array([[2.5,2.5,3.5],[2.5,3.5,3.5]])]

corridor = [np.linspace([1,1,1],[1,1,2.5],3),
            np.linspace([2.5,2.5,3.5],[2.5,3.5,3.5], 3)]


# start      = [0,0,0]
# target    = [0,0,5]
# corridor   = [np.array([[2,2,2],[2,2,3.5]])]

corridor_r = 0.15  

order = 6
dim = 3
T = 6
dt = .08
m = int(T/dt)+1


corridor_info = {}
corridor_info['end_points'] = corridor
corridor_info['radius'] = corridor_r



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

drone = Drone_traj(start, target, corridor_info, order, T, priority=0)
drone.traj_bezier.set_degree_dim(order, dim, idx_w=[4])
planner = Planner(corridor_r)

ts = time.time()
planner.get_primary_traj(drone)
print(f"QP takes {time.time()-ts:.4f} sec")

starts = [start]
targets = [target]
corridors = [corridor]
drones = [drone]

view_angle=[-31, 34]
planner.plot_final_frame(drones, corridors, view_angle)



