from utils.pkg_func import *

import utils.toolbox_bezier as bzr



class Bezier:
    def __init__(self, degree, dim, T) -> None:
        assert(degree >= 3)
        
        self.set_degree_dim(degree, dim)
        self.T = T
        
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
    
    def get_A(self, time, derivative=0):
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
    




