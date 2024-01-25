from utils.pkg_func import *


def bernstein_base(n,i,l):
    r"""
    Evaluates a basis associated with the `i`-th control point of an
    `n` degree Bezier curve at a locaiton :math:`\lambda\in[0,1]`.
    
    :param n: Degree of the considered Bezier curve
    :type  n: int greater than zero, typically in [1,10]
    :param i: Index of considered control point
    :type  i: int on the interval [0,n]
    :param l: Location at which the basis is evaluated
    :type  l: array on the interval [0,1]

    :returns:
        - a \in R^l.size() - array of evaluation of :math:`\alpha_n^i(\lambda)` 
    """
    return comb(n,i) * l**i * (1.0-l)**(n-i)
    
    # return (factorial(n)/(factorial(i)*factorial(n-i)))*(1-l)**(n-i)*l**i
    
    
def evaluate_bezier(P,lam):
    """
    Evaluates a bezier curve :math:`\mathbf{B}(\lambda; \mathcal{P})` of degree 
    :math:`n`  the control points :math:`\mathcal{P}` provided in a vector representation
    :math:`\mathbf{p}=\mathrm{mat}(\mathcal{P})\in\mathbb{R}^{d,(n+1)}` at the points
    :math:`\lambda\in[0,1]^m`
    
    :param P: Matrix representation :math:`\mathbf{v}=\mathrm{mat}(\mathcal{P})` of the control points
    :type  P: ndarray of dimension R(d,n+1)
    :param lam: Locations :math:`\lambda` at which the curve is to be evaluated
    :type  lam: ndarray of dimension R(m,)

    :returns:
        - B - ndarray matrix in R(d,M)
    """
    n, d = P.shape
    n -= 1
    B = np.zeros((len(lam), d))
    for i in range(n+1):
        B += bernstein_base(n,i,lam).reshape((-1,1))@P[i,:].reshape((1,-1))
    return B


def get_bernstein_mat(n,d,lam):
    """
    get bernstein matrix in the structure of 
         P0   P1  ... Pn
       | |||  ||| ... |||
       t |||  ||| ... |||
       | |||  ||| ... |||
       in R^{len(lam),d*(n+1)}
    """
    B = np.zeros((len(lam), d*(n+1)))
    for i in range(n+1):
        B[:,i*d:(i+1)*d] = np.tile(bernstein_base(n,i,lam),(d,1)).T
        # for j in range(d):
        #     B[:, i+(n+1)*j] = bernstein_base(n,i,lam)
    return B
    
  

def eveluate_point_DeCas(self, t:float):
    """
    evaluate the curve at a specific moment t using De Caseljau's algorithm
    """
    assert t>=0 and t<= 1.0

    P_new = self.P_big
    d = self.dim
    P_new[:,0] = self.get_Pvec().squeeze()
    
    for j in range(1, self.degree+1):   # j-th column
        for i in range(self.degree+1-j):  # i*dim:(i+1)*dim-th row
            P_new[i*d:(i+1)*d, j] = (1-t)*P_new[i*d:(i+1)*d,j-1] + t*P_new[(i+1)*d:(i+2)*d,j-1]
        
    return P_new[:d, -1]
    
    
def get_M1(n, T):
    """
    Compute the linear map :math:`\mathbf{M}_n^1` relating the one dimensional bezier curve
    :math:`\mathbf{B}(tT^{-1}; \mathcal{P})` with its stderivative :math:`\mathbf{B}(tT^{-1}; \mathcal{P}^{(1)})`
    in the sense that :math:`\mathrm{vec}(\mathcal{P}^{(1)}) = \mathbf{M}_n^1\mathrm{vec}(\mathcal{P})`

    :param n: Degree of the considered Bezier curve
    :type  n: int greater than zero, typically in [1,10]
    :param T: Time over which the spline is defined
    :type  T: float greater than zero, typically in [0.5,1.5] 

    :returns:
        - M - ndarray matrix in R(n,n+1)
    """
    return  (n / T) * (np.c_[np.zeros((n,1)),np.eye(n)]-
                       np.c_[np.eye(n),np.zeros((n,1))])

def get_M(n,k,d,T):
    """
    Recursively evaluate the map :math:`\mathbf{M}_n^k` relating a
    :math:`d`-dimensional bezier curve :math:`\mathbf{B}(tT^{-1}; \mathcal{P})` with its 
    :math:`k`-th derivative :math:`\mathbf{B}(tT^{-1}; \mathcal{P}^{(k)})` 
    in the sense that :math:`\mathrm{vec}(\mathcal{P}^{(k)}) = \mathbf{M}_n^k\mathrm{vec}(\mathcal{P})`

    :param n: Degree of the considered Bezier curve
    :type  n: int greater than zero, typically in [1,10]
    :param k: Derivative order
    :type  k: int in the interval [0,n]
    :param d: Dimension of the considered Bezier curve
    :type  d: int greater than zero, typically in [1,3]
    :param T: Time over which the spline is defined
    :type  T: float greater than zero, typically in [0.5,1.5] 

    :returns:
        - M - ndarray matrix in R(n+1-k,n+1)
    """
    M1D = np.eye(n+1)
    for ii in range(k):
        M1D = get_M1(n-ii,T)@M1D
    return np.kron(M1D, np.eye(d))



def get_c(i,j,m,n):
    r"""
    Evaluates the integral :math:`\int_0^1\alpha_n^i(\lambda)\alpha_m^j(\lambda) d\lambda` (*)

    :param n: Degree of the considered Bezier curve
    :type  n: int greater than zero, typically in [1,10]
    :param m: Degree of the considered Bezier curve
    :type  n: int greater than zero, typically in [1,10]
    :param i: Index of a control point
    :type  i: int on the interval [0,n]
    :param j: Index of a control point
    :type  j: int on the interval [0,m]

    :returns:
        - c - float evaluation of (*)
    """
    return (factorial(m)/(factorial(i)*factorial(m-i)))*\
           (factorial(n)/(factorial(j)*factorial(n-j)))*\
           gamma(i+j+1)*gamma(m+n-i-j+1)/gamma(m+n+2)

def get_C(m,n,d):
    """
    Evaluates a matrix with elements comprising of the integral exoressions
    computed in `get_c`, which is used in forming the varuous :math:`\mathcal{L}_2`
    norm objectives.

    :param n: Degree of the considered Bezier curve
    :type  n: int greater than zero, typically in [1,10]
    :param m: Degree of the considered Bezier curve
    :type  n: int greater than zero, typically in [1,10]
    :param d: Dimension of the considered Bezier curve
    :type  d: int greater than zero, typically in [1,3]

    :returns:
        - C - ndarray of dimension R(m+1,n+1)
    """
    C1D = np.zeros((m+1,n+1))
    for ii in range(m+1):
        for jj in range(n+1):
            C1D[ii,jj] = get_c(ii,jj,m,n)
    return np.kron(C1D, np.eye(d))

def get_snap_cost_bezier(n,d,c,T):
    """
    Decompose a squared :math:`\mathcal{L}_2` cost related to total variation
    of a Bezier curve :math:`\mathbf{B}(tT^{-1}; \mathcal{P})` of degree `n` and
    and its derivatives. The cost function of interest is

    :math:`\sum\limits_{k=0}^n c_k\| (d^k/dt^k)\mathbf{B}(tT^{-1}; \mathcal{P})\|_{\mathcal{L}_2([0,T])}^2`

    which can be written as a quadratic in :math:`\mathbf{p}=\mathrm{vec}(\mathcal{P})`, as

    :math:`J(\mathbf{p}) = \mathbf{p}^T\mathbf{Q}\mathbf{p}`

    Note, this function is tested and validated in `testBezierSegment.py`.

    :param n: Degree of the considered Bezier curve
    :type  n: int greater than zero, typically in [1,10]
    :param d: Dimension of the considered Bezier curve
    :type  d: int greater than zero, typically in [1,3]
    :param c: Derivative weights in the cost function
    :type  c: ndarray of dimension R(n+1,)
    :param T: Time over which the spline is defined
    :type  T: float greater than zero, typically in [0.5,1.5] 

    :returns:
        - Q - Positive definite ndarray matrix in R(d(N+1),d(N+1))
    """
    Q = get_C(n,n,d)* c[0] * T
    for k in range(1,n+1):
        Mnkd = get_M(n,k,d, T)
        Qn   = get_C(n-k,n-k,d)
        Q    += T * (Mnkd.T @ Qn @ Mnkd) * c[k]
    return Q



def get_snap_cost_bezier_derivatives(n,d,c,T,P):
    nP = (n+1)*d
    grad_N_P = np.zeros((nP,nP))
    grad_N_T = np.zeros((nP,nP))
    hess_N_T = np.zeros((nP,nP))
    for k in range(0,n+1):
        Mnkd      = get_M(n,k,d,1.0)
        Qn        = get_C(n-k,n-k,d)
        Qtmp      = c[k] * (Mnkd.T@Qn@Mnkd)
        grad_N_P += Qtmp * (T ** (1.0-2.0*k))
        grad_N_T += Qtmp * (T ** (0.0-2.0*k)) * (1.0-2.0*k)
        if k > 0:
            hess_N_T += (Qtmp * (T ** (-1.0-2.0*k))) * ((1.0-2.0*k) * (-2.0*k))
    
    gP, gT = P.T@grad_N_T@P, 2*grad_N_P@P
    hPP, hPT, hTT = 2*grad_N_P, 2*grad_N_T@P, P.T@hess_N_T@P
    return gP, gT, hPP, hPT, hTT


def get_wind_cost_bezier_derivatives(N,d,w,q,T,V,P):
    HP = np.zeros(((N+1)*d,(N+1)*d))
    HV = np.zeros(((N+1)*d,(N+1)*d))
    f  = np.zeros((1,(N+1)*d))

    nP = (N+1)*d
    gP = np.zeros((nP, 1))
    gT = np.zeros(( 1, 1))
    for a in range(0,N+1):
        MNad = get_M(N,a,d,T)
        for b in range(0,N+1):
            MNbd = get_M(N,b,d,T)
            Cmnd = get_C(N-a,N-b,d)
            Cbar = (MNad.T @ Cmnd @ MNbd) * T
            HP  += Cbar * w[a] * w[b]
            f   += V.T @ Cbar * w[b] * q[a] + V.T @ Cbar.T * w[a] * q[b]
            HV  += Cbar * q[a] * q[b]
    gP = T*(2*HP@P + f.T)
    gT = 1.0
    return gP, gT

def get_wind_cost_bezier(N,d,w,q,T,V):
    """
    Decompose the special squared :math:`\mathcal{L}_2` cost related to the wind
    formed in two Bezier curves :math:`\mathbf{B}(tT^{-1}; \mathcal{P})` and
    :math:`\mathbf{B}(tT^{-1}; \mathcal{V})` of the same degree, here :math:`N`,
    and is defined with a weight vector, here :math:`\mathbf{w}=(w_0,...,w_N)`.
    The cost function of interest is

    :math:`\Bigg\| \sum\limits_{k=0}^N (d^k/dt^k)[w_k \mathbf{B}(tT^{-1}; \mathcal{P}) + q_k\mathbf{B}(tT^{-1}; \mathcal{V})] \Bigg\|_{\mathcal{L}_2([0,T])}^2`

    which, rather amazingly can be written as a quadratic in :math:`\mathbf{p}=\mathrm{vec}(\mathcal{P})`, as

    :math:`J(\mathbf{p}; \mathcal{V}) = \mathbf{p}^T\mathbf{H}\mathbf{p} + \mathbf{f}^T\mathbf{p} + c`

    Note, this function is tested and validated in `testBezierSegment.py`.

    :param N: Degree of the considered Bezier curve
    :type  N: int greater than zero, typically in [1,10]
    :param d: Dimension of the considered Bezier curve
    :type  d: int greater than zero, typically in [1,3]
    :param w: Derivative weights in the thrust cost function
    :type  w: ndarray of dimension R(N+1,)
    :param q: Derivative weights in the thrust cost function
    :type  q: ndarray of dimension R(N+1,)
    :param T: Time over which the spline is defined
    :type  T: float greater than zero, typically in [0.5,1.5] 
    :param V: Vector representation :math:`\mathbf{v}=\mathrm{vec}(\mathcal{V})` of the control points defining the wind speeds
    :type  V: ndarray of dimension R(d(N+1),1)

    :returns:
        - H - Positive semidefinite ndarray matrix in R(d(N+1),d(N+1))
        - f - Ndarray in R(1,d(N+1))
        - c - Constant term in (**) mainly used to sanity check implementation
    """
    HP = np.zeros(((N+1)*d,(N+1)*d))
    HV = np.zeros(((N+1)*d,(N+1)*d))
    f  = np.zeros((1,(N+1)*d))

    for a in range(0,N+1):
        MNad = get_M(N,a,d,T)
        for b in range(0,N+1):
            MNbd = get_M(N,b,d,T)
            Cmnd = get_C(N-a,N-b,d)
            Cbar = (MNad.T @ Cmnd @ MNbd)
            HP  += Cbar * w[a] * w[b]
            f   += V.T @ Cbar * w[b] * q[a] + V.T @ Cbar.T * w[a] * q[b]
            HV  += Cbar * q[a] * q[b]
    H  = HP * T
    f *= T
    c  = (V.T @ HV @ V) * T
    return H, f, c

def convert_spline(spline, resolution=100, method='regression'):
    """
    Convert a time-scaled bezier curve to a standard polynomial curve,
    either using regression or by analytical expressions.

    TODO: Implement and test the analytical version
    TODO: Migrate to the bezier curve class instead

    :param spline: A spline object of type Bezier
    :type  BezierSpline: A bezier spline object
    :param resolution: The resolution used in the regression (if applicable)
    :type  resolution: int greater than zero
    :param method: The method used for conversion, regression is generally more stable for higher degrees
    :type  method: string, either "regression" or "analytical

    :returns:
        - PolynomialSpline - A polynomial spline object
    """
    P = spline.P.reshape([-1,spline.d],order='F')
    
    if method == 'regression':
        l = np.linspace(0,1,resolution)
        C = np.zeros((P.shape[0],P.shape[1]))
        n = P.shape[1]-1
        Bval = evaluate_bezier(P,l)
        for idx in range(spline.d):
            C[idx, :] = np.polynomial.polynomial.Polynomial.fit(spline.T*l, Bval[idx,:], n).convert().coef
    else:
        SyntaxError("Not yet implemented. Exists in Matlab.")
    return C.reshape((-1), order='F')


def get_regression_cost_bezier(lvec, Dmat, n):
    """
    TODO: Add docstring

    :param lvec: 
    :type  lvec: 
    :param Dmat: 
    :type  Dmat: 
    :param n: Degree of the considered Bezier curve
    :type  n: int greater than zero, typically in [1,10]

    :returns:
        - H - Positive definite ndarray matrix in R(d(N+1),d(n+1))
        - f - Ndarray in R(1,d(n+1))
        - c - Constant term in (**) mainly used to sanity check implementation
    """
    
    d = Dmat.shape[0]
    L = Dmat.shape[1]

    # Build total variation cost maps
    delta = np.array([lvec[1]-lvec[0]])
    delta = np.append(delta, (lvec[2:]-lvec[:-2])/2)
    delta = np.append(delta, lvec[-1]-lvec[-2])

    N1D = np.diag(delta)
    A1D = np.zeros((L,n+1))

    for k in range(L):
        for i in range(n+1):
            A1D[k,i] = bernstein_base(n,i,lvec[k])

    N = np.kron(N1D, np.eye(d))
    A = np.kron(A1D, np.eye(d))
    D = np.reshape(Dmat, [d*L], order='F')
    H =    A.T@N@A
    f = -2*A.T@N@D
    c =    D.T@N@D
    
    return H, f, c