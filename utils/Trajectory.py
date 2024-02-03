from utils.pkg_func import *

import utils.toolbox_bezier as bzr


############################# Spline super class ##############################

class SplineSegment():
    """
    Abstract spline segment class, contains shadow functions that are
    implemented in inheriting classes, and also contains various
    plotting functionality

    :param d: Trajectory dimension
    :type  d: int > 0 
    :param n: Trajectory order
    :type  n: int > 0 
    :param id: The identifier of the spline, mainly used for debugging
    :type  id: int > 0 

    :raises TypeError: When inputs are if incorrect type
    :raises ValueError: When inputs are if incorrect value
    """
    def __init__(self, d=1, n=6, id = 0, T=2.0):
        self.id          = is_type(id, int)              # Spline identifier
        self.d           = is_positive(is_type(d, int))  # Spline dimension
        self.n           = is_positive(is_type(n, int))  # Spline order
        self.isOptimized = False                         # Flag set to False no data has been written to {P,T}
        self.T           = T                          # Spline time
        self.P           = np.zeros((d*(n+1),1))         # Spline data
        self.poly_n      = self.d   

    def set_T(self, T):
        self.T = T
        # raise BrokenPipeError("Not yet implemented")

    def get_Q(self, cost):
        raise BrokenPipeError("Not yet implemented")

    def get_A(self):
        raise BrokenPipeError("Not yet implemented")

    def get_indices(self, dim):
        raise BrokenPipeError("Not yet implemented")
        
    def evaluate_in_time(self, time, derivative):
        raise BrokenPipeError("Not yet implemented")

    
############################ Bezier segment class #############################
class Bezier:
    def __init__(self, degree, dim, T, id=0) -> None:
        """
        self.P = [p0x p0y p0z | p1x p1y p1z | .. ].T
        """
        assert(degree >= 3)
        
        self.set_degree_dim(degree, dim)
        self.T = T
        self.id = id
        
    def set_degree_dim(self, degree, dim, idx_w=[4]):
        self.n = degree
        self.d = dim
        self.P = np.zeros((self.d*(self.n+1),1))
        self.idx_w = idx_w
        self.set_weights(self.idx_w)
        
        
    def set_weights(self, idx_w):
        derivative_w = np.zeros(self.n+1)
        derivative_w[idx_w] = np.ones(len(idx_w))
        self.weights = derivative_w
        
        
    def set_P(self, P : np.ndarray) -> None:
        """ 
        vectorized control points
        structure: [p0x p0y p0z | p1x p1y p1z | .. ].T
        """
        assert P.size == self.d*(self.n+1)
        
        if len(P.shape) == 1:
            self.P = P.reshape((-1,1))
        elif P.shape[1] == self.d:
            # given P: [[p0x p0y p0z]
            #           [p1x p1y p1z]
            #           [... ... ...]]
            self.P = P.reshape((-1,1))
        else:
            raise ValueError('Provided control potins shape: ', P.shape, "Wrong dimension of control points")
        
        # self.P_big = np.zeros(((self.degree+1)*self.dim, self.degree+1))
        
    
    
    def get_Pmat(self) -> np.ndarray:
        """ 
        matrix control points
        structure: [[p0x p0y p0z]
                    [p1x p1y p1z]
                    [... ... ...]]
        """
        return self.P.reshape((-1,self.d), order='F')
    
    
    def get_indices(self, dim):
        """
        Get the polynomial coefficients of a given dimension

        :param dim: Dimension to extract
        :type  dim: int in [0,d-1]

        :return: Polynomial coefficients in dimension "dim" (going from 0 to n)
        :rtype: ndarray in R(n+1)
        """
        return [ii for ii in range(dim, (self.n+1)*self.d, self.d)]

    
    def get_Q(self, cost):
        """ 
        matrix for snap cost P.T@Q@P
        """
        # Q = bzr.get_snap_cost_bezier(self.n, self.d, cost, self.T)
        # D = np.eye((self.n+1)*self.d) - \
        #         np.diag(np.ones((self.n*self.d)),self.d)
        # np.delete(D, D.shape[0]-np.arange(self.d)-1, 0)
        # return Q + self.dcost*(D.T@D)
        return bzr.get_snap_cost_bezier(self.n, self.d, cost, self.T)
            
        
    def evaluate_in_time(self, time, derivative=0):
        """
        evaluate bezier curve at given time points with derivative
        
        :param time: time points to be evaluated
        :param derivative: at -th order of derivative
        
        :return values 
        """
        # Input checks
        if isinstance(time, (np.ndarray, np.generic)):
            time = np.reshape(time, (len(time),))
        elif isinstance(time, (int, float)):
            time = np.array([time])
        else:
            raise TypeError("Input type not recognized")

        if isinstance(derivative,int):
            if derivative >= 0 and derivative <= self.n+1:
                derivative = int(derivative)
            else:
                raise ValueError("Derivative should be within [0,n+1]")
        else:
            raise TypeError("Derivative type should be int")

        # Get map from original points to the derivative
        M  = bzr.get_M(self.n, derivative, self.d, self.T)
        # Get points of the derivative bezier curve
        P  = M @ self.P
        # Convert points to a dx(n+1) dimensional array
        # CP = P.reshape([-1, self.d], order='F')
        CP = P.reshape([-1, self.d])
        
        # Define the time over which the spline is to be evaluated
        L  = time / self.T
        return bzr.evaluate_bezier(CP, L)
    
    def evaluate_in_space(self, time):
        return self.evaluate_in_time(time, derivative=0)
    
    def get_A_old(self, time, derivative=0):
        """
        at t \in time and derivative, find 
        Aeq@MP=beq  or Aun@MP<=bun 
        
        in the structure of 
          P0   P1  ... Pn 
        | .    .   ... .
      t0|  .    .  ...  .
       _|   .    . ...   .
        | .    .   ... .
      t1|  .    .  ...  .
       _|   .    . ...   .

        :param time \in R^nt
        :return A \in R^{len(time)*d,d*(n+1)}
        """
        assert self.T >= np.max(time)
        
        # Get map from original points to the derivative
        M  = bzr.get_M(self.n, derivative, self.d, self.T)
        # Define the time over which the spline is to be evaluated
        L  = time / self.T
        # get bernstein base matrix at derivative
        berMat = bzr.get_bernstein_mat(self.n-derivative, self.d, L)@M
        
        Nr,Nc = berMat.shape
        A = np.zeros((Nr*self.d, Nc))
        d, n = self.d, self.n+1
        for i in range(Nr):
            for j in range(0,Nc,d):
                A[i*d:(i+1)*d, j:j+d] = np.diag(berMat[i,j:j+d])
        
        return A
    
 
    def get_A(self, conti_order):
        Ap = np.zeros((self.d*(self.n+1),self.d*(self.n+1)))
        Am = np.zeros_like(Ap)
        for k in range(conti_order+1):
            Mk = bzr.get_M(self.n,k,self.d,self.T)
            Nm = np.c_[np.zeros((self.d, Mk.shape[0]-self.d)), np.eye(self.d)]
            Np = np.c_[np.eye(self.d), np.zeros((self.d, Mk.shape[0]-self.d))]
            # if not k:
            #     Ap = Np@Mk.copy()
            #     Am = Nm@Mk.copy()
            # else:
            #     Ap = np.r_[Ap, Np@Mk]
            #     Am = np.r_[Am, Nm@Mk]
            
            Ap[np.arange(k*self.d,(k+1)*self.d),:] = Np@Mk
            Am[np.arange(k*self.d,(k+1)*self.d),:] = Nm@Mk            

        return Ap, Am   # weird order, due to its final implementation
    
   
    def get_idx_constr(self, time):
        """
        get the index of constraints
        :return : equality idx (start, target), inequality idx
        """
        N_time = len(time)
        
        base = np.arange(self.d*N_time)
        eq_idx = np.r_[base[:self.d], base[-self.d:]]
        
        ineq_idx = np.setdiff1d(base, eq_idx)
        
        return eq_idx, ineq_idx

    def set_T(self, T):
        self.T = T
        
        
    def get_polynomial_representation(self):
        """
        Get a polynomial representation of the spline

        :return: A polynomial spline of the same dimensions as the Bezier spline
        :rtype: PolynomialSegment
        """
        polynomial_spline   = PolynomialSegment(self.d, self.n)
        polynomial_spline.P = bzr.convert_spline(self).reshape((-1,1), order='F')
        polynomial_spline.T = self.T
        return polynomial_spline
    


############################# Trajectory class ################################

class Trajectory():
    """
    A Trajectory class that contions multiple Bezier curves
    
    :param degree: Bezier curve order
    :type degree: int > 0
    :param dim: trajectory dimension
    :type dim: int > 0
    :param m: Number of Bezier used to define the Trajectory
    :type m: int > 0
    :param conti_order: up to the order that is ensured to be continuoued
    :type conti_order: int > 0
    
    self.P = [ |    |  ..  |
              B0.P B1.P.. Bm.P
               |    |  ..  |]
    Bi.P = [p0x p0y p0z | p1x p1y p1z | .. | pnx pny pnz ].T
    """
    def __init__(self, degree, dim, start, target, corridor, T, conti_order=4) -> None:
        self.n = degree
        self.d = dim
        self.m = sum([len(cor) for cor in corridor])+1
        self.T = T 
        
        self.constraints = {}
        self.set_waypoints(start, target, corridor)
        self.constraints['conti_order'] = conti_order  # Continuity in first kmax derivatives
        self.constraints['CC'] = np.empty(((self.n+1)*self.d,(self.m)))*np.nan
        self.constraints['CC'][:(conti_order+1)*self.d,:] = np.zeros(((conti_order+1)*self.d,self.m))

        self.constraints['useFEP'] = 1      # Use fixed-end-point constraints
        self.constraints['useCON'] = 1      # Use continuity constraints

        self.initial_spline(T)
        

    def set_waypoints(self, start, target, corridor):
        self.constraints['P0'] = np.empty(((self.n+1)*self.d,(self.m)))*np.nan
        self.constraints['PT'] = np.empty(((self.n+1)*self.d,(self.m)))*np.nan
        
        # set beginning point for each segment
        self.constraints['P0'][:self.d,0] = start        
        i = 1
        for cor in corridor:
            self.constraints['P0'][:self.d , i:i+len(cor)] = cor.T
            i += len(cor)
            
        # set ending point for the last segment (set to target)
        self.constraints['PT'][:self.d,-1] = target
        
        # start and end at 0 velocity and accerlation
        self.constraints['P0'][self.d:self.d*3,0] = np.zeros((2*self.d))
        self.constraints['PT'][self.d:self.d*3,-1] = np.zeros((2*self.d))
        

    def initial_spline(self, T):
        wp = np.c_[self.constraints['P0'][:self.d,:], self.constraints['PT'][:self.d,-1]]
        wp_dis = np.linalg.norm(np.diff(wp, axis=1), axis=0)
        total_dis = sum(wp_dis)
        
        self.splines = [Bezier(self.n, self.d, T=wp_dis[id]/total_dis*T, id=id) for id in range(self.m)]
        
        
    def set_T(self, T):
        for idx, spline in enumerate(self.splines):
            spline.set_T(T[idx])
            

    def get_T_ind(self):
        """
        Compute the domain of the trajectory as extracted from the splines

        :return: T list of times over which the individual spliens are defined
        :rtype: list of floats
        """
        return [self.splines[ii].T for ii in range(self.m)]
            

    def get_T_cum(self):
        """
        Compute the domain of the trajectory as extracted from the splines

        :return: T list of times over which the trajectory is defined, starting from 0
        :rtype: list of floats

        """
        return np.r_[0,np.cumsum(self.get_T_ind())]
        
    def get_indices(self, id):
        I = id*(self.n+1)*self.d + np.arange(self.d*(self.n+1))
        Ix,Iy = np.meshgrid(I,I,indexing='ij')
        return I,Ix,Iy

    def get_mesh_indices(self, Ix, Iy):
        Ix,Iy = np.meshgrid(Ix,Iy,indexing='ij')

        return Ix,Iy


    def get_spline(self,id):
        """
        Get the spline with the desired id

        : param id: id defined to find the spline
        : type id: int

        :return: spline with the specified id
        :rtype: SplineSegment
        """
        return self.splines[id]
    
        
    def set_P(self, P : np.ndarray) -> None:       
        P_mat = P.reshape((self.m,-1))
        for i in range(self.m):
            self.splines[i].set_P(P_mat[i])
        
    
    def set_degree_dim(self, degree, dim, idx_w=[4]):
        self.d = dim
        self.n = degree
        for sp in self.splines:
            sp.set_degree_dim(degree, dim, idx_w=idx_w)
            
    def get_snap_cost(self):
        """
        Compute the Hessian matrix forming the quadratic L2-cost function associated
        minimum snap optimization for the entire trajectory

        :return: Q matrix for the entire trajectory
        :rtype: numpy.array

        """
        H_full = np.zeros(((self.n+1)*self.d*self.m,(self.n+1)*self.d*self.m))
        for index, spline in enumerate(self.splines):

            Iv,Ix,Iy = self.get_indices(index)

            # Add the minimum-snap cost (if applicable)
            H = spline.get_Q(spline.weights)
            H_full[Ix,Iy] += H
            
        return H_full
            
    def get_EPC(self):
        """
        end points constraints
        Aeq = [[Aeq0]             ]
              [       [Aeq1]      ]
              [             ..    ]  
              [             [Aeqm]]
        beq follows the similar strcture
        """
        n, m, d, P0, PT, conti_order = self.n,  self.m, self.d, \
                            self.constraints['P0'], self.constraints['PT'], self.constraints['conti_order']

        # ----------- Construct the end point constraints --------------
        # At the initial point of each spline
        A_0_ep = np.zeros(((n+1)*m*d, (n+1)*m*d))
        b_0_ep = np.zeros(((n+1)*m*d))*np.nan
        # At the end-point of each spline
        A_T_ep = np.zeros(((n+1)*m*d, (n+1)*m*d))
        b_T_ep = np.zeros(((n+1)*m*d))*np.nan
        
        for spline in self.splines:
            id = spline.id
            I,Ix,Iy = self.get_indices(id)
            [A0, AT] = spline.get_A(self.n)
            # Form constraint matrices
            A_T_ep[Ix,Iy] = AT
            A_0_ep[Ix,Iy] = A0
            # Define constraints
            b_0_ep[I]   = P0[:,id]
            b_T_ep[I]   = PT[:,id]               
            
        # Remove all constraints that were not set
        A_0_ep = A_0_ep[~np.isnan(b_0_ep),:]
        A_T_ep = A_T_ep[~np.isnan(b_T_ep),:]
        b_0_ep = b_0_ep[~np.isnan(b_0_ep)]
        b_T_ep = b_T_ep[~np.isnan(b_T_ep)]
        
        A_ep = np.r_[A_0_ep, A_T_ep]
        b_ep = np.r_[b_0_ep, b_T_ep]

        return A_ep,b_ep
    
    def get_EPC(self):
        """
        end points constraints
        Aeq = [[Aeq0]             ]
              [       [Aeq1]      ]
              [             ..    ]  
              [             [Aeqm]]
        beq follows the similar strcture
        """
        n, m, d, P0, PT, conti_order = self.n,  self.m, self.d, \
                            self.constraints['P0'], self.constraints['PT'], self.constraints['conti_order']

        # ----------- Construct the end point constraints --------------
        # At the initial point of each spline
        A_0_ep = np.zeros(((n+1)*m*d, (n+1)*m*d))
        b_0_ep = np.zeros(((n+1)*m*d))*np.nan
        # At the end-point of each spline
        A_T_ep = np.zeros(((n+1)*m*d, (n+1)*m*d))
        b_T_ep = np.zeros(((n+1)*m*d))*np.nan
        
        # At the initial point of each spline
        A_0_unep = np.zeros(((n+1)*m*d, (n+1)*m*d))
        b_0_unep = np.zeros(((n+1)*m*d))*np.nan
        # At the end-point of each spline
        A_T_unep = np.zeros(((n+1)*m*d, (n+1)*m*d))
        b_T_unep = np.zeros(((n+1)*m*d))*np.nan
        
        for spline in self.splines:
            id = spline.id
            I,Ix,Iy = self.get_indices(id)
            [A0, AT] = spline.get_A(self.n*0)
            # only start and target points are equally constrained
            if id == 0:
                # Form constraint matrices
                A_0_ep[Ix,Iy] = A0
                # Define constraints
                b_0_ep[I]   = P0[:,id]
            elif id == m:
                A_T_ep[Ix,Iy] = AT
                b_T_ep[I]     = PT[:,id]
                A_0_unep[Ix,Iy] = A0
                b_0_unep[I]   = P0[:,id]
                
            # end points of tennels are bounded (unequally constrained)
            else:
                # Form constraint matrices
                A_T_unep[Ix,Iy] = AT
                A_0_unep[Ix,Iy] = A0
                # Define constraints
                b_0_unep[I]   = P0[:,id]
                b_T_unep[I]   = PT[:,id]
                
        # Remove all constraints that were not set
        A_0_ep = A_0_ep[~np.isnan(b_0_ep),:]
        A_T_ep = A_T_ep[~np.isnan(b_T_ep),:]
        b_0_ep = b_0_ep[~np.isnan(b_0_ep)]
        b_T_ep = b_T_ep[~np.isnan(b_T_ep)]
        
        A_0_unep = A_0_unep[~np.isnan(b_0_unep),:]
        A_T_unep = A_T_unep[~np.isnan(b_T_unep),:]
        b_0_unep = b_0_unep[~np.isnan(b_0_unep)]
        b_T_unep = b_T_unep[~np.isnan(b_T_unep)]
        
        A_ep = np.r_[A_0_ep, A_T_ep]
        b_ep = np.r_[b_0_ep, b_T_ep]
        Aun_ep = np.r_[A_0_unep, A_T_unep]
        bun_ep = np.r_[b_0_unep, b_T_unep]

        return A_ep,b_ep, Aun_ep, bun_ep
    
    
    def get_CC(self):   # continuous constraints
        n, m, d, CC = self.n, self.m, self.d, self.constraints['CC']
        A_con = np.zeros(((n+1)*((m-1)*d), (n+1)*m*d))
        b_con = np.zeros((n+1)*((m-1)*d))

        for id in range(self.m-1):
            Im,Imx,Imy = self.get_indices(id)
            Ip,Ipx,Ipy = self.get_indices(id+1)

            # manipulate the index for higher dimension
            [_, AT] = self.get_spline(id).get_A(self.constraints['conti_order'])
            [A0, _] = self.get_spline(id+1).get_A(self.constraints['conti_order'])

            A_con[Imx, Imy] = +AT
            A_con[Imx, Ipy] = -A0
            b_con[Im]       = CC[:,id]
            
        A_con = A_con[~np.isnan(b_con),:]
        b_con = b_con[~np.isnan(b_con)]

        return A_con,b_con
        
        
    def evaluate_in_time(self, time, derivative, resolution=50):
        """
        Evaluate the trajectory at a given time, pututting the terminal point
        if it is greater than the domain over which the trajectory is deifned
        and the initial point if the time is smaller than the domain over which
        the trajectory is defined.

        :param time: Time at which the trajectory is to be evaluated
        :type  time: float

        :return: Array of dimension (d,1)
        :rtype: Numpy array
        """
        T_cum = self.get_T_cum()        

        id_s = np.where(T_cum<=time[0])[0][-1] if time[0]>0 else 0
        id_e = np.where(T_cum>time[-1])[0][0]-1 if time[-1]<T_cum[-1] else len(T_cum)-1

        # T_cum = self.get_T_cum(eval=True)        

        # ensure the time duration for bezier is T \in [0,1]
        # if self.spline_type == 'bezier':
        #     T_cum = np.arange(self.m+1)
            
        values = np.zeros(((derivative+1)*self.d,(id_e-id_s)*resolution))
        t_span = np.zeros((id_e-id_s)*resolution)

        for id in range(-id_s+id_e):
            spline = self.get_spline(id+id_s)
            tt = np.linspace(0,spline.T,resolution)
            for dd in range(derivative+1):
                val = spline.evaluate_in_time(tt,dd).T
                values[dd*self.d:(dd+1)*self.d , id*resolution:(id+1)*resolution] = val
            t_span[id*resolution:(id+1)*resolution] = (tt + T_cum[id+id_s])

        return values, t_span