from utils.pkg_func import *
from utils.Trajectory import Bezier
import utils.toolbox_visualize as vis
import utils.toolbox_casdi as caF
import casadi as ca


class Planner():
    def __init__(self, corridor_r) -> None:
        self.optimized = False
        self.corridor_r = corridor_r
        
    def get_primary_traj(self, drone):
        uneq_con = False

        traj = drone.traj_bezier
        corridor_r = self.corridor_r

        opti = ca.Opti()
        X = opti.variable(traj.d*(traj.n+1)*traj.m)
        Tspan = [sp.T for sp in traj.splines]
        
        H_full = traj.get_snap_cost()
        Aeq  = np.zeros((0,(traj.n+1)*traj.d*traj.m))
        beq  = np.zeros(0)
        Auneq = np.zeros_like(Aeq) 
        buneq = np.zeros_like(beq)
        if traj.constraints['useFEP']:
            if not uneq_con:
                Aeq_ep, beq_ep = traj.get_EPC(Tspan)
                Aeq  = np.r_[Aeq, Aeq_ep]
                beq  = np.r_[beq, beq_ep]
            else:
                Aeq_ep, beq_ep, Aun_eq, bun_eq = traj.get_EPC_uneq(Tspan)
                Aeq  = np.r_[Aeq, Aeq_ep]
                beq  = np.r_[beq, beq_ep]
                Auneq = np.r_[Auneq, Aun_eq]
                buneq = np.r_[buneq, bun_eq]
        if traj.constraints['useCON']:
            A_con,b_con = traj.get_CC(Tspan)
            Aeq = np.r_[Aeq, A_con]
            beq = np.r_[beq, b_con]
            
        # ----------------- construct optimization solver ----------------

        opti.minimize(X.T@H_full@X)
        opti.subject_to(Aeq@X == beq)

        # opti.minimize(X.T@H_full@X + np.ones((1,traj.m))@Tspan)
        # opti.subject_to(Aeq@X == beq)
        # opti.subject_to(np.eye(traj.m)@Tspan >= np.ones((traj.m, 1))*.1)
        # opti.subject_to(np.ones((1,traj.m))@Tspan >= traj.T/2)

        if uneq_con:
            dis = (Auneq@X-buneq).reshape((-1, traj.d))
            for i in range(dis.shape[0]):
                opti.subject_to(ca.norm_2(dis[i, :]) <= corridor_r)

        critical_pts = np.c_[traj.constraints['P0'][:traj.d,:], traj.constraints['PT'][:traj.d,-1]]
        P_init = np.zeros(0)
        for i in range(critical_pts.shape[1]-1):
            P_init = np.r_[P_init, np.linspace(critical_pts[:,i], critical_pts[:,i+1], traj.n+1).reshape((-1,))]
        # P_init = np.linspace(drone.start, drone.target, (traj.n+1)*traj.m).reshape((-1,1))

        opti.set_initial(X, P_init)

        opts = {}
        opts = {"ipopt.print_level":0, "print_time": False, 'ipopt.max_iter':100}
        opti.solver("ipopt", opts)
        
        ts = time.time()
        sol = opti.solve()
        print(f"\n\nQP takes {time.time()-ts:.4f} sec")

        traj.set_P(sol.value(X))
        # drone.traj_pt, t_span = traj.evaluate_in_time([0, traj.get_T_cum()[-1]], derivative=0)
        # drone.traj_pt = drone.traj_pt.T  


    def get_primary_traj_optiT(self, drone, kT):
        traj = drone.traj_bezier
        corridor_r = self.corridor_r

        opti = ca.Opti()
        X = opti.variable(traj.d*(traj.n+1)*traj.m)
        Tspan = opti.variable(traj.m)
        
        H_full = caF.get_snap_cost(traj)
        if traj.constraints['useFEP']:
            Aeq_ep, beq_ep = caF.get_EPC(traj, Tspan)
        if traj.constraints['useCON']:
            A_con, b_con = caF.get_CC(traj, Tspan)
            Aeq = ca.vcat([Aeq_ep, A_con])
            beq = ca.vcat([beq_ep, b_con])
            
        # ----------------- construct optimization solver ----------------

        # opti.minimize(X.T@H_full@X)
        # opti.subject_to(Aeq@X == beq)

        T_init = np.array([sp.T for sp in traj.splines])
        opti.minimize(X.T@H_full@X + kT*np.ones((1,traj.m))@Tspan)
        opti.subject_to(Aeq@X == beq)
        opti.subject_to(np.eye(traj.m)@Tspan >= np.ones((traj.m,1))*.2)
        # opti.subject_to(np.ones((1,traj.m))@Tspan >= traj.T/3)

        critical_pts = np.c_[traj.constraints['P0'][:traj.d,:], traj.constraints['PT'][:traj.d,-1]]
        P_init = np.zeros(0)
        for i in range(critical_pts.shape[1]-1):
            P_init = np.r_[P_init, np.linspace(critical_pts[:,i], critical_pts[:,i+1], traj.n+1).reshape((-1,))]
        # P_init = np.linspace(drone.start, drone.target, (traj.n+1)*traj.m).reshape((-1,1))

        opti.set_initial(X, P_init)
        opti.set_initial(Tspan, T_init)


        opts = {}
        opts = {"ipopt.print_level":0, "print_time": False, 'ipopt.max_iter':100}
        opti.solver("ipopt", opts)
        
        ts = time.time()
        sol = opti.solve()
        print(f"\n\nQP takes {time.time()-ts:.4f} sec")
        

        traj.set_P(sol.value(X))
        drone.set_Ts(sol.value(Tspan))
        

    @timer
    def get_timer_map(self, trajs, T, collide_r, idx=[0,1],  method=0, spy=1):
        N = len(trajs[0])
        t_span = np.linspace(0,T,N)
        timer_map = np.zeros((N,N))
        
        assert(len(trajs)==2)
        if method == 0:     # brute force 1: O(N^2)
            for i in range(N):
                for j in range(N):
                    if np.linalg.norm(trajs[0][i]-trajs[1][j])<=collide_r:
                        timer_map[i,j] = 1
        elif method == 1:   # brute force 2: O(N), worst O(N^2)
            for offset in range(10):
                
                dis = np.linalg.norm(trajs[0]-trajs[1], axis=1)
                idx = np.where(dis<collide_r)[0]
                print(idx)
        if spy:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            
            ax.spy(timer_map, origin='lower')
            ax.set_xlabel(f'drone {idx[0]} time line')
            ax.set_ylabel(f'drone {idx[1]} time line')

        return timer_map

    

    def update_timer_map(self, drones, timer_map, pt, update_traj = 1):
        # TODO: 1. find critical time points using optimization
        N = drones[0].num_pt
        T = drones[0].T
        pt = np.r_[pt, [[N,N]]]
        time_dif = np.diff(pt,axis=0)
        num_per_seg = np.max(time_dif, axis=1)
        # num_per_seg = [50, 50]
        trajs = []
        # dilution
        for j in range(2):
            time_samples = np.empty((0,))
            for i in range(len(pt)-1):
                time_samples = np.r_[time_samples, np.linspace(pt[i,j], pt[i+1,j], num_per_seg[i])[:-1]*T/N]
            
            time_samples = np.r_[time_samples, T]
            
            trajs.append(drones[j].ctrlPt.evaluate_in_time(time_samples))

        if update_traj:
            for i in range(len(trajs)):
                drones[i].traj_pt = trajs[i] 
                
        return trajs
        

    