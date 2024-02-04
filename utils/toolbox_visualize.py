###############################################################################
# CONFIDENTIAL (C) Mitsubishi Electric Research Labs (MERL) 2023
# Author: Yunshen Huang
# Date: 01/10/2023
#
# TODO: Write description
###############################################################################
from utils.pkg_func import *

def palette():
    cm = plt.cm.colors.CSS4_COLORS
    colors = ['cyan', 'green', 'navy', 'black', 'tomato', 'red', 'darkkhaki']
    return cm, colors


def visualize_splines_1D(values, d, times, ax_list=None, lineWidth=2, \
    plotEndpoints=True, lineColor=None,figsize=(7,5), verbose={}):

    derivative_name = ['position', 'velocity', 'acceleration', 'jerk']
    ax_name = ['X','Y','Z']
    lc_list = ['b','g','r','c','m']
    # need values to be a fat matrix
    if values.shape[0] > values.shape[1]:
        values = values.T
    to_order = int(values.shape[0]/d)-1

    
    def plot_endpoints(tt, vals, lc, plotEndpoints,ax):
        if plotEndpoints:
            # print(tt[0], tt[-1])

            ax.plot(tt[0],   vals[0], color=lc, marker='o')
            ax.plot(tt[-1], vals[-1], color=lc, marker='x')

    def plot_1segment(times,values,last_plot=True):

        for dimension in range(d):
            ax = ax_list[dimension] if d>1 else ax_list

            line_labels = [] 

            for ii in range(to_order+1):

                if not lineColor:
                    lc = lc_list[ii]
                else:
                    lc = lineColor

                line_label, = ax.plot(times,values[dimension+d*ii,:], color=lc, linewidth=lineWidth)
                
                if last_plot: line_labels.append(line_label) 

                plot_endpoints(times,values[dimension+d*ii,:],lc,plotEndpoints, ax)
            
            if last_plot:
                ax.legend(line_labels,derivative_name[:to_order+1],loc=0)

                # ax.set_ylabel(r'$(\mathrm{d}^n/\mathrm{d}t^n)P(t) [s]$')
                ax.set_ylabel(ax_name[dimension])
                
            ax.grid()


    # if not len(ax_list): 
    #     fig, ax_list = plt.subplots(nrows=d, ncols=1, figsize=figsize)
    fig, ax_list = add_ax(ax_list, nrows=d)

    if len(values.shape) == 3:
        for ss in range(values.shape[2]):
            last_plot = True if ss == values.shape[2]-1 else False
                
            tt = times[:,ss]
            val = values[:,:,ss]
            plot_1segment(tt,val,last_plot)
    else:
        plot_1segment(times,values)
    
    ax_first = ax_list[0] if d>1 else ax_list
    ax_last = ax_list[-1] if d>1 else ax_list
    try:
        ax_first.set_title(f"Polynomials (id [{verbose['id'][0]}-{verbose['id'][-1]}]):{verbose['m']} splines with {verbose['n']} degrees")
    except:
        pass
    ax_last.set_xlabel('Time [s]')

    # fig.tight_layout(pad=2.0)
    ax_last.get_figure().tight_layout(pad=2.0)

    return fig, ax_list
            

def vis_bezier_regression(Dmat,Best,ax=None):
    fig, ax = add_ax_3d(figsize=[8,8])
    s1 = ax.scatter3D(Dmat[0,:],Dmat[1,:],Dmat[2,:],s=5,c='r')
    l2, = ax.plot3D(Best[0,:],Best[1,:],Best[2,:],linewidth=2,c='g')
    ax.set_title(f'Estimate error: {np.linalg.norm(Dmat-Best)/Dmat.shape[1]}')
    plt.legend([l2,s1],['Estimated curve','Sampled points'])
    
    return fig, ax


def corridor_geo(corridor, radius):
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

    
def drone_geo():
    """
    geometry of drone
    # TODO: finish function 
    """
    
    return 0

    
def plot_final_frame_3D(drones, corridor_info, view_angle, verbose):
    trajs       = [d.pos_pt for d in drones]
    targets     = [d.target for d in drones]
    starts      = [d.start for d in drones]
    corridors    = corridor_info['end_points']
    corridor_r  = corridor_info['radius']
            
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    fig, ax = add_ax_3d()

    cm, colors = palette()
    for i, (curve, start, target), in enumerate(zip(trajs, starts, targets,)):
        ax.plot(curve[:,0], curve[:,1], curve[:,2], label=f"Drone {i+1}", c=cm[colors[i]])
        ax.scatter(start[0],start[1], start[2], c=cm[colors[i]], s=25, marker='o')
        ax.scatter(target[0],target[1], target[2], c=cm[colors[i]], s=25, marker='x')
        
    for cors in corridors:
        # cor_axes = [np.nonzero(c[0,:]-c[1,:])[0] for c in cors]
        # for (cor, cor_ax) in zip(cors, cor_axes):
        #     X,Y,Z = corridor_geo(cor,corridor_r)
        #     ax.plot_surface(X,Y,Z, color=cm['darkkhaki'], alpha=.3)
        X,Y,Z = corridor_geo(cors[[0,-1]],corridor_r)
        ax.plot_surface(X,Y,Z, color=cm['darkkhaki'], alpha=.3)
            

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.view_init(view_angle[1], view_angle[0])

    ax.legend()
    set_axes_equal(ax)
    plt.tight_layout()
    
    if verbose['show']:
        plt.show()
    
    return fig, ax


def plot_animation_3D(drones, corridor_info, view_angle, verbose):
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

    
    fig, ax = plot_final_frame_3D(drones, corridor_info, view_angle, {'show':False})
    trajs = [d.pos_pt for d in drones]
    T = drones[0].T
    [ax.lines.pop(0) for _ in range(len(ax.lines))]

    cm, colors = palette()
    
    t_all = np.linspace(0,T,len(trajs[0]))
    lines = [ax.plot([], [], [], lw=2, label=f"Drone {i+1}", color=cm[colors[i]])[0] for i in range(len(trajs))]
    objs = [ax.plot([], [], [], c=cm[colors[i]], lw=5, marker='o',linestyle="")[0] for i in range(len(trajs))]
    titleTime = ax.text2D(0.1, 0.95, "", transform=ax.transAxes)
            
    line_ani = animation.FuncAnimation(fig, updateLines, init_func=ini_plot, frames=len(trajs[0]), interval=(T/len(trajs[0])*1000), blit=False)
    
    if verbose['save']:
        line_ani.save('./archive/2dronesOpti2.gif', writer='PillowWriter', fps=len(trajs[0])/T)
    if verbose['show']:
        plt.show()
    
    return line_ani
    

def add_ax_3d(fig=None, ax=None,figsize=[8,8]):
    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(projection='3d')
    return fig, ax
    
def add_ax(fig=None, ax=None, nrows=3, ncols=1, figsize=(10, 6)):
    if ax is None:
        fig, ax_list = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)
    return fig, ax_list
    
    
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
    

