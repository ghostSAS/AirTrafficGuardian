from utils.pkg_func import *
from utils.Bezier import Bezier
import casadi as ca


class Planner():
    def __init__(self) -> None:
        self.cm = plt.cm.colors.CSS4_COLORS
        self.colors = ['cyan', 'green', 'navy', 'black', 'tomato', 'red', 'darkkhaki']
        self.optimized = False
        
        
    def corridor_geo(self, corridor, radius, axis):
        """
        geometry of corridor, corridor aligns with the given axis
        :param: axis (int)  
            0: x axis
            1: y axis
            2: z axis
        # TODO: finish function 
        """
        assert(axis in [0,1,2])
        
        height = np.ptp(corridor[:,axis])
        center = np.mean(corridor, axis=0)
        phi = np.linspace(0, 2 * np.pi, 25)
        z = np.linspace(-height/2, height/2, 2) + center[axis]

        Z, PHI = np.meshgrid(z, phi)
        X = radius * np.cos(PHI) + center[axis-2]
        Y = radius * np.sin(PHI) + center[axis-1]
        
        if axis == 0:
            return Z,X,Y
        elif axis == 1:
             return Y,Z,X
        else:
            return X,Y,Z
        
        
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
    

    def update_timer_map(self, timer_map):
        # TODO: impliment using optimization
        pt = [40, 12]



    
    def plot_final_frame(self, trajs, starts, targets, corridors, corridor_r, view_angle):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        cm, colors = self.cm, self.colors
        l, s, t = [], [], []
        for i, (curve, start, target), in enumerate(zip(trajs, starts, targets,)):
            l.append(ax.plot(curve[:,0], curve[:,1], curve[:,2], label=f"Drone {i+1}", c=cm[colors[i]]))
            s.append(ax.scatter(start[0],start[1], start[2], c=cm[colors[i]], s=25, marker='o'))
            t.append(ax.scatter(target[0],target[1], target[2], c=cm[colors[i]], s=25, marker='x'))
            
        for cors in corridors:
            cor_axes = [np.nonzero(c[0,:]-c[1,:])[0] for c in cors]
            for (cor, cor_ax) in zip(cors, cor_axes):
                X,Y,Z = self.corridor_geo(cor,corridor_r,cor_ax)
                ax.plot_surface(X,Y,Z, color=cm['darkkhaki'], alpha=.3)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.view_init(view_angle[1], view_angle[0])

        ax.legend()
        set_axes_equal(ax)
        plt.tight_layout()
        
        return fig, ax

    
    
    def plot_animation(self, trajs, starts, targets, corridors, corridor_r, T, view_angle):
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

        fig, ax = self.plot_final_frame(trajs, starts, targets, corridors, corridor_r, view_angle)
        [ax.lines.pop(0) for _ in range(len(ax.lines))]

        cm, colors = self.cm, self.colors
        
        t_all = np.linspace(0,T,len(trajs[0]))
        lines = [ax.plot([], [], [], lw=2, label=f"Drone {i+1}", color=cm[colors[i]])[0] for i in range(len(trajs))]
        objs = [ax.plot([], [], [], c=cm[colors[i]], lw=5, marker='o',linestyle="")[0] for i in range(len(trajs))]
        titleTime = ax.text2D(0.1, 0.95, "", transform=ax.transAxes)
                
        # 
        line_ani = animation.FuncAnimation(fig, updateLines, init_func=ini_plot, frames=len(trajs[0]), interval=(T/len(trajs[0])*1000), blit=False)
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
    