from utils.pkg_func import *
import utils.toolbox_bezier as bzr
import casadi as ca



def get_A(spline, conti_order, T):
        d, n = spline.d, spline.n

        Ap = ca.MX(d*(n+1), d*(n+1))
        Am = ca.MX(d*(n+1), d*(n+1))
        for k in range(conti_order+1):
            Mk = bzr.get_M(spline.n, k, spline.d, T)
            Nm = np.c_[np.zeros((spline.d, Mk.shape[0]-spline.d)), np.eye(spline.d)]
            Np = np.c_[np.eye(spline.d), np.zeros((spline.d, Mk.shape[0]-spline.d))]
            # if not k:
            #     Ap = Np@Mk.copy()
            #     Am = Nm@Mk.copy()
            # else:
            #     Ap = np.r_[Ap, Np@Mk]
            #     Am = np.r_[Am, Nm@Mk]
            
            # Ap[np.arange(k*self.d,(k+1)*self.d),:] = Np@Mk
            # Am[np.arange(k*self.d,(k+1)*self.d),:] = Nm@Mk            
            Ap[k*spline.d:(k+1)*spline.d,:] = Np@Mk
            Am[k*spline.d:(k+1)*spline.d,:] = Nm@Mk   
        return Ap, Am   # weird order, due to its final implementation
    

def get_snap_cost(traj):
    """
    Compute the Hessian matrix forming the quadratic L2-cost function associated
    minimum snap optimization for the entire trajectory

    :return: Q matrix for the entire trajectory
    :rtype: numpy.array

    """
    H_full = ca.MX((traj.n+1)*traj.d*traj.m, (traj.n+1)*traj.d*traj.m)
    for index, spline in enumerate(traj.splines):

        Iv,Ix,Iy = traj.get_indices(index)

        # Add the minimum-snap cost (if applicable)
        H = spline.get_Q(spline.weights)
        H_full[Iv[0]:Iv[-1]+1, Iv[0]:Iv[-1]+1] = H
        
    return H_full

def get_EPC(traj, Tspan):
    """
    end points constraints
    Aeq = [[Aeq0]             ]
            [       [Aeq1]      ]
            [             ..    ]  
            [             [Aeqm]]
    beq follows the similar strcture
    """
    n, m, d, P0, PT, conti_order = traj.n,  traj.m, traj.d, \
                        traj.constraints['P0'], traj.constraints['PT'], traj.constraints['conti_order']

    # ----------- Construct the end point constraints --------------
    # At the initial point of each spline
    A_0_ep = ca.MX((n+1)*m*d, (n+1)*m*d)
    b_0_ep = np.zeros(((n+1)*m*d))*np.nan
    # At the end-point of each spline
    A_T_ep = ca.MX((n+1)*m*d, (n+1)*m*d)
    b_T_ep = np.zeros(((n+1)*m*d))*np.nan
    
    for spline in traj.splines:
        id = spline.id
        Iv,Ix,Iy = traj.get_indices(id)
        [A0, AT] = get_A(spline, spline.n, Tspan[id])
        # Form constraint matrices
        A_T_ep[Iv[0]:Iv[-1]+1, Iv[0]:Iv[-1]+1] = AT
        A_0_ep[Iv[0]:Iv[-1]+1, Iv[0]:Iv[-1]+1] = A0
        # Define constraints
        b_0_ep[Iv]   = P0[:,id]
        b_T_ep[Iv]   = PT[:,id]               
        
    # Remove all constraints that were not set
    A_0_ep = A_0_ep[np.where(~np.isnan(b_0_ep))[0],:]
    A_T_ep = A_T_ep[np.where(~np.isnan(b_T_ep))[0],:]
    b_0_ep = b_0_ep[~np.isnan(b_0_ep)]
    b_T_ep = b_T_ep[~np.isnan(b_T_ep)]
    
    A_ep = ca.vertcat(A_0_ep, A_T_ep)
    b_ep = ca.MX(np.r_[b_0_ep, b_T_ep])

    return A_ep,b_ep


def get_CC(traj, Tspan):   # continuous constraints
    n, m, d, CC = traj.n, traj.m, traj.d, traj.constraints['CC']
    A_con = ca.MX((n+1)*((m-1)*d), (n+1)*m*d)
    b_con = np.zeros((n+1)*((m-1)*d))

    for id in range(traj.m-1):
        Im,Imx,Imy = traj.get_indices(id)
        Ip,Ipx,Ipy = traj.get_indices(id+1)

        # manipulate the index for higher dimension
        [_, AT] = get_A(traj.get_spline(id), traj.constraints['conti_order'], Tspan[id])
        [A0, _] = get_A(traj.get_spline(id+1), traj.constraints['conti_order'], Tspan[id+1])

        A_con[Im[0]:Im[-1]+1, Im[0]:Im[-1]+1] = AT
        A_con[Im[0]:Im[-1]+1, Ip[0]:Ip[-1]+1] = -1*A0
        b_con[Im]       = CC[:,id]
        
    A_con = A_con[np.where(~np.isnan(b_con))[0],:]
    b_con = ca.MX(b_con[~np.isnan(b_con)])

    return A_con,b_con