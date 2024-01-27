from utils.pkg_func import *
from utils.Bezier import Bezier
from utils.Drone import Drone
from utils.Planner import Planner


num_drones = 2
corridor_r = 0.15
collision_r = corridor_r*3
T = 4

"""
task: fly towards each other
    drone 1 [0,0,0] --> [0,5,5]
    drone 2 [0,5,0] --> [0,0,5]
    corridor [2,2,2] --> [2,2,3.5], r=0.15
    
goal:
    1. reach the target
    2. they shall not collide
    3. minimum snap
"""
starts      = [[0,0,0], [0,5,0]]
targets     = [[0,5,5], [0,0,5]]
corridors   = [[np.array([[2,2,2],[2,2,3.5]])],  
               [np.array([[2,2,3.5],[2,2,2]])]]
critical_time_pt = np.array([[0, 0],
                    [42, 11],
                   [51, 51]])

# """
# task: fly along the same direction
#     drone 1 [0,0,0] --> [4,4,0]
#     drone 2 [0,0,5] --> [4,4,5]
#     corridor [2,2,2] --> [2,2,3.5], r=0.15
    
# goal:
#     1. reach the target
#     2. they shall not collide
#     3. minimum snap
# """
# starts      = [[0,0,0], [4,4,0]]
# targets     = [[0,0,5], [4,4,5]]
# corridors   = [[np.array([[2,2,2],[2,2,3.5]])],  
#                [np.array([[2,2,2],[2,2,3.5]])]]
# critical_time_pt = np.array([[0, 0],
#                     [21,17],
#                     [34,30]])


planner = Planner(corridor_r)
prim_trajs = []
scheduled_trajs = []
drones = []

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
    drones.append(drone)
    
"""
step2 
    find the colliding time points
"""
timer_map = planner.get_timer_map(prim_trajs, T, collision_r, method=0, spy=0)


"""
step3
    modify timeline based on priority
"""
trajs = planner.update_timer_map(drones, timer_map, critical_time_pt, update_traj = 1)

timer_map = planner.get_timer_map(trajs, T, collision_r, method=0, spy=0)


"""
step4
    validate no collision
"""
for (i, j) in list(combinations_of_2(range(len(drones)))):
    dis = np.linalg.norm(trajs[i]-trajs[j], axis=1)
    collide_idx = np.where(dis<=collision_r)[0]
    if len(collide_idx)>0:
        print(f"Collision found, algorithm fails! \n \t\t at t={collide_idx}, with distance of {dis[collide_idx]}")
    else:
        print("No collision found")
    
    
"""
plot trajectories
    either final frame, or animation
"""
# drones[0].traj_pt = prim_trajs[0][42:,:]
# drones[1].traj_pt = prim_trajs[1][12:,:]
view_angle=[-31, 34]

# fig, ax = planner.plot_final_frame(drones, corridors, view_angle)
ani = planner.plot_animation(drones, corridors, view_angle, show_now=1)


