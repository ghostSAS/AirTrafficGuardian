from utils.pkg_func import *
from utils.Bezier import Bezier
from utils.Drone import Drone
from utils.Planner import Planner



"""
task: 
    drone 1 [0,0,0] --> [0,5,5]
    drone 2 [0,5,0] --> [0,0,5]
    corridor [2,2,2] --> [2,2,3.5], r=0.15
    
goal:
    1. reach the target
    2. they shall not collide
    3. minimum snap
"""

num_drones = 2

starts      = [[0,0,0], [0,5,0]]
targets     = [[0,5,5], [0,0,5]]
corridors   = [[np.array([[2,2,2],[2,2,3.5]])],  
               [np.array([[2,2,3.5],[2,2,2]])]]
corridor_r = 0.15
T = 4

planner = Planner()
prim_trajs = []
scheduled_trajs = []

"""
step1
    based on specified constraints for each drone, solve multiple QPs to get individual trajectory
"""
for i in range(num_drones):
    start = starts[i]
    target = targets[i]
    corridor = corridors[i]
    
    corridor_info = {}
    corridor_info['end_points'] = corridor
    corridor_info['radius'] = corridor_r
    
    drone = Drone(start, target, corridor_info, T, priority=i)
    drone.get_primary_traj()
    prim_trajs.append(drone.traj_pt)

"""
step2 
    find the colliding time points
"""
timer_map = planner.get_timer_map(prim_trajs, T, corridor_r*3, method=0, spy=0)


"""
step3
    modify timeline based on priority
"""


"""
plot trajectories
    either final frame, or animation
"""
view_angle=[-31, 34]
# fig, ax = planner.plot_final_frame(fig, ax, trajs, starts, targets, corridors, corridor_r, view_angle)
ani = planner.plot_animation(prim_trajs, starts, targets, corridors, corridor_r, T, view_angle)

plt.show()

