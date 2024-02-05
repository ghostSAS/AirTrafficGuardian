from utils.pkg_func import *
from utils.Trajectory import Bezier
from utils.Drone import *
from utils.Planner import Planner
import utils.toolbox_visualize as vis



num_drones = 2
corridor_r = 0.15
collision_r = corridor_r*3

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
corridors   = [[np.linspace([2,2,2],[2,2,3.5],3)],  
               [np.linspace([2,2,3.5],[2,2,2],3)]]
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


# ------------- parameters ---------------
optimal_ord_idx = [3, 4]
kT=70
optiT = True

order = 5
dim = 3
T = 6

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
    
    drone = Drone_traj(start, target, corridor_info, order, T, priority=i, idx_w=optimal_ord_idx)

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

    traj_pt, t_span = drone.get_traj_pt(derivative=2, resolution=150)
    prim_trajs.append(drone.pos_pt)
    drones.append(drone)
    vis.visualize_splines_1D(traj_pt, len(start), t_span)
    
view_angle=[-31, 34]
vis.plot_final_frame_3D(drones, corridor_info, view_angle, verbose={'show':False, 'save':False})
# vis.plot_animation_3D(drones, corridor_info, view_angle, verbose={'show':True, 'save':False})


"""
step2 
    find the colliding time points
"""
timer_map = planner.get_timer_map(prim_trajs, T, collision_r, method=0, spy=1)


# """
# step3
#     modify timeline based on priority
# """
# trajs = planner.update_timer_map(drones, timer_map, critical_time_pt, update_traj = 1)

# timer_map = planner.get_timer_map(trajs, T, collision_r, method=0, spy=0)


# """
# step4
#     validate no collision
# """
# for (i, j) in list(combinations_of_2(range(len(drones)))):
#     dis = np.linalg.norm(trajs[i]-trajs[j], axis=1)
#     collide_idx = np.where(dis<=collision_r)[0]
#     if len(collide_idx)>0:
#         print(f"Collision found, algorithm fails! \n \t\t at t={collide_idx}, with distance of {dis[collide_idx]}")
#     else:
#         print("No collision found")
    
    
# """
# plot trajectories
#     either final frame, or animation
# """
# # drones[0].traj_pt = prim_trajs[0][42:,:]
# # drones[1].traj_pt = prim_trajs[1][12:,:]
# view_angle=[-31, 34]

# # fig, ax = planner.plot_final_frame(drones, corridors, view_angle)
# ani = planner.plot_animation(drones, corridors, view_angle, verbose=1)


