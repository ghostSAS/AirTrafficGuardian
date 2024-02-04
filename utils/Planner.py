from utils.pkg_func import *
from utils.Trajectory import Bezier
import casadi as ca


class Planner():
    def __init__(self, corridor_r) -> None:
        self.cm = plt.cm.colors.CSS4_COLORS
        self.colors = ['cyan', 'green', 'navy', 'black', 'tomato', 'red', 'darkkhaki']
        self.optimized = False
        self.corridor_r = corridor_r
        
    def get_primary_traj(self, drone):
        uneq_con = False

        traj = drone.traj_bezier
        corridor_r = self.corridor_r
        
        H = traj.get_snap_cost()
        Aeq  = np.zeros((0,(traj.n+1)*traj.d*traj.m))
        beq  = np.zeros(0)
        Auneq = np.zeros_like(Aeq) 
        buneq = np.zeros_like(beq)
        if traj.constraints['useFEP']:
            if not uneq_con:
                Aeq_ep, beq_ep = traj.get_EPC()
                Aeq  = np.r_[Aeq, Aeq_ep]
                beq  = np.r_[beq, beq_ep]
            else:
                Aeq_ep, beq_ep, Aun_eq, bun_eq = traj.get_EPC_uneq()
                Aeq  = np.r_[Aeq, Aeq_ep]
                beq  = np.r_[beq, beq_ep]
                Auneq = np.r_[Auneq, Aun_eq]
                buneq = np.r_[buneq, bun_eq]
        if traj.constraints['useCON']:
            A_con,b_con = traj.get_CC()
            Aeq = np.r_[Aeq, A_con]
            beq = np.r_[beq, b_con]
            
        # ----------------- construct optimization solver ----------------
        opti = ca.Opti()

        X = opti.variable(traj.d*(traj.n+1)*traj.m)
        Tspan = opti.variable(traj.m)

        # opti.minimize(X.T@H@X)
        # opti.subject_to(Aeq@X == beq)

        opti.minimize(X.T@H@X + np.ones((1,traj.m))@Tspan)
        opti.subject_to(Aeq@X == beq)
        opti.subject_to(np.eye(traj.m)@Tspan >= np.ones((traj.m, 1))*.1)
        opti.subject_to(np.ones((1,traj.m))@Tspan >= traj.T/2)

        if uneq_con:
            dis = (Auneq@X-buneq).reshape((-1, traj.d))
            for i in range(dis.shape[0]):
                opti.subject_to(ca.norm_2(dis[i, :]) <= corridor_r)

        P_init = np.linspace(drone.start, drone.target, (traj.n+1)*traj.m).reshape((-1,1))
        T_init = [sp.T for sp in traj.splines]

        opti.set_initial(X, P_init)
        opti.set_initial(Tspan, T_init)


        opts = {"ipopt.print_level":0, "print_time": False, 'ipopt.max_iter':100}
        opts = {}
        opti.solver("ipopt", opts)
        sol = opti.solve()

        traj.set_P(sol.value(X))
        drone.traj_pt, t_span = traj.evaluate_in_time([0, traj.get_T_cum()[-1]], derivative=0)
        drone.traj_pt = drone.traj_pt.T  
        
        
    def corridor_geo(self, corridor, radius):
        height = np.linalg.norm((corridor[0,:]-corridor[-1,:]))
        center = np.mean(corridor, axis=0)
        phi = np.linspace(0, 2 * np.pi, 25)
        z = np.linspace(-height/2, height/2, 2)

        Z, PHI = np.meshgrid(z, phi)
        X = radius * np.cos(PHI)
        Y = radius * np.sin(PHI)
        
        a = np.array([0,0,1])
        b = (corridor[0,:]-corridor[-1,:])/height
        c = a@b     # cosine of angle
        # s = np.linalg.norm(v)    # sine of angle
        if np.abs(c-1)<1e-5 or np.abs(c+1)<1e-5:
            return X+center[0], Y+center[1], Z+center[2]
        else:
            v = np.cross(a,b)
            vx = np.array([[0, -v[2], v[1]],
                           [v[2], 0, -v[0]],
                           [-v[1], v[0], 0]])
            R = np.eye(3) + vx + vx@vx/(1+c)
            
            XYZ = np.c_[X.reshape(-1), Y.reshape(-1), Z.reshape(-1)].T
            XYZ_rot = R@XYZ
            return XYZ_rot[0].reshape(X.shape)+center[0], XYZ_rot[1].reshape(Y.shape)+center[1], XYZ_rot[2].reshape(Z.shape)+center[2]

        
    def drone_geo(self):
        """
        geometry of drone
        # TODO: finish function 
        """
        
        return 0
    
    @timer
    def get_timer_map(self, trajs, T, collide_r, method=0, spy=1):
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
            ax.set_xlabel('drone 1 time line')
            ax.set_ylabel('drone 2 time line')

            plt.show()
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
        

    
    def plot_final_frame(self, drones, corridors, view_angle, show_now=True):
        
        trajs = [d.traj_pt for d in drones]
        targets = [d.target for d in drones]
        starts = [d.start for d in drones]
        corridor_r = self.corridor_r
                
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        cm, colors = self.cm, self.colors
        for i, (curve, start, target), in enumerate(zip(trajs, starts, targets,)):
            ax.plot(curve[:,0], curve[:,1], curve[:,2], label=f"Drone {i+1}", c=cm[colors[i]])
            ax.scatter(start[0],start[1], start[2], c=cm[colors[i]], s=25, marker='o')
            ax.scatter(target[0],target[1], target[2], c=cm[colors[i]], s=25, marker='x')
            
        for cors in corridors:
            cor_axes = [np.nonzero(c[0,:]-c[1,:])[0] for c in cors]
            for (cor, cor_ax) in zip(cors, cor_axes):
                X,Y,Z = self.corridor_geo(cor,corridor_r)
                ax.plot_surface(X,Y,Z, color=cm['darkkhaki'], alpha=.3)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.view_init(view_angle[1], view_angle[0])

        ax.legend()
        set_axes_equal(ax)
        plt.tight_layout()
        
        if show_now:
            plt.show()
        
        return fig, ax

    
    
    def plot_animation(self, drones, corridors, view_angle, show_now = 1):
        if show_now:
            def updateLines(i):
                for j, curve, in enumerate(trajs):
                    lines[j].set_data(curve[:i+1,0], curve[:i+1,1])
                    lines[j].set_3d_properties(curve[:i+1,2])
                    objs[j].set_data(curve[i:i+1,0], curve[i:i+1,1])
                    objs[j].set_3d_properties(curve[i:i+1,2])
                titleTime.set_text(u"Time = {:.2f} s".format(t_all[i]))
                return lines, objs

            def ini_plot():
                for i in range(len(trajs)):
                    curve = trajs[i]
                    ax.plot(curve[:,0], curve[:,1], curve[:,2], c=cm[colors[i]], lw=1, linestyle = "--")

                for i in range(len(lines)):
                    lines[i].set_data(np.empty([1]), np.empty([1]))
                    lines[i].set_3d_properties(np.empty([1]))
                    objs[i].set_data(np.empty([1]), np.empty([1]))
                    objs[i].set_3d_properties(np.empty([1]))

            fig, ax = self.plot_final_frame(drones, corridors, view_angle, show_now=False)
            trajs = [d.traj_pt for d in drones]
            T = drones[0].T
            [ax.lines.pop(0) for _ in range(len(ax.lines))]

            cm, colors = self.cm, self.colors
            
            t_all = np.linspace(0,T,len(trajs[0]))
            lines = [ax.plot([], [], [], lw=2, label=f"Drone {i+1}", color=cm[colors[i]])[0] for i in range(len(trajs))]
            objs = [ax.plot([], [], [], c=cm[colors[i]], lw=5, marker='o',linestyle="")[0] for i in range(len(trajs))]
            titleTime = ax.text2D(0.1, 0.95, "", transform=ax.transAxes)
                    
            line_ani = animation.FuncAnimation(fig, updateLines, init_func=ini_plot, frames=len(trajs[0]), interval=(T/len(trajs[0])*1000), blit=False)
            line_ani.save('./archive/2dronesOpti2.gif', writer='PillowWriter', fps=len(trajs[0])/T)
            
            plt.show()
            return line_ani
    
    
    
def set_axes_equal(ax):
    """
    Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    """

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])
    