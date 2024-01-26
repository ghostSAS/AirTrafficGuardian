from utils.pkg_func import *
from utils.Bezier import Bezier
import casadi as ca



class Drone():
    
    def __init__(self, start, target, corridor_info, T, priority) -> None:
        self.corridor = corridor_info['end_points']
        self.corridor_r = corridor_info['radius']
        # self.process_corridor(corridor_info)
        self.start = np.array(start)
        self.target = np.array(target)
        self.T = T
        self.optimized = False
        self.ctrlPt = Bezier(7,len(start),T)
        self.num_pt = 51
        self.traj_pt = np.empty((self.num_pt, len(start)))
        self.priority = priority
        
    def process_corridor(self, corridor_info):
        self.corridor_axis = [np.nonzero(c[0,:]-c[1,:])[0] for c in self.corridor]
        
        
    def get_samples(self):
        critical_pts = np.zeros((2+len(self.corridor)*2, self.start.size))
        critical_pts[0,:] = self.start
        for i in range(len(self.corridor)):
            critical_pts[i*2+1:(i+1)*2+1,:] = self.corridor[i]
        critical_pts[-1,:] = self.target

        total_dis = 0
        total_dis_list = []
        for i in range(len(critical_pts)-1):
            dis = np.linalg.norm(critical_pts[i+1] - critical_pts[i])
            total_dis += dis
            total_dis_list.append(total_dis)

        num_pt_seg = 5
        time_samples = np.zeros((0))
        corridor_samples = np.zeros((0,self.corridor[0].shape[1]))
        for c in self.corridor:
            corridor_samples = np.r_[corridor_samples, np.linspace(c[0], c[1], num_pt_seg)]
        total_dis = total_dis_list[-1]
        for i, dis in enumerate(total_dis_list):
            if i%2:
                Ts = total_dis_list[i-1]/total_dis*self.T
                Te = (total_dis_list[i]-total_dis_list[i-1])/total_dis*self.T
                time_samples = np.r_[time_samples, np.linspace(Ts, Te+Ts, num_pt_seg)]
        
        return corridor_samples, time_samples

        
    def get_primary_traj(self):
        traj = self.ctrlPt
        start = self.start
        target = self.target
        dim = self.ctrlPt.d
        corridor_r = self.corridor_r
        
        corridor_sample, t_sample = self.get_samples()
        
        # ---------------------------
        # t_sample = [0]
        # t_sample.extend(t_sample_tmp)
        # t_sample.append(self.T)
        # # t_sample = np.linspace(0,self.T,50)[[0,20,24,28,32,-1]]
        # A = traj.get_A(np.array(t_sample))


        # eq_idx, ineq_idx = traj.get_idx_constr(t_sample)

        # # A_eq@P = [start[x,y,z].T, target[x,y,z].T] = b_eq
        # A_eq = A[eq_idx, :]
        # b_eq = np.r_[start, target]

        # # A_ineq@P <= b_ineq
        # A_ineq = A[ineq_idx, :]

        # ---------------------------
        # corridor_sample = np.linspace(corridor[0], corridor[1], len(t_sample)-2)

        A_eq = traj.get_A(np.array([0, self.T]))
        b_eq = np.r_[start, target]
        A_ineq = traj.get_A(t_sample)
        

        # ----------------- construct optimization solver ----------------
        opti = ca.Opti()

        X = opti.variable(traj.d*(traj.n+1))

        # constraints 1 within corridor
        opti.subject_to(A_eq@X==b_eq)
        for i in range(len(corridor_sample)):
            opti.subject_to(ca.norm_2(A_ineq[i*dim:(i+1)*dim]@X - corridor_sample[i]) <= corridor_r)

            
        # constraints 2 == corridor sample
        # opti.subject_to(A_eq@X==b_eq)
        # opti.subject_to(A_ineq@X==corridor_sample.reshape((-1)))

        opti.minimize(X.T@traj.get_Q(traj.weights)@X)

        P_init = np.linspace(start, target, traj.n+1).reshape((-1,1))
        opti.set_initial(X,P_init)

        opts = {"ipopt.print_level":0, "print_time": False, 'ipopt.max_iter':100}
        opti.solver("ipopt", opts)
        sol = opti.solve()

        self.ctrlPt.set_P(sol.value(X))
        self.traj_pt = traj.evaluate_in_time(np.linspace(0,self.T,self.num_pt))  



        
        

        
    