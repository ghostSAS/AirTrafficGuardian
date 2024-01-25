import sys
sys.path.append('.')

from utils.pkg_func import *
from utils.Bezier import Bezier



start = np.array([0,0,0])
target = np.array([5,5,5])
corridor = np.array([[2,2,2],[2,2,3.5]])
corridor_r = 0.15

degree = 10
dim = 3
T = 2
dt = .08
m = int(T/dt)+1

# examine n-th other of derivative
derivative = 0

# ----------------- create bezier -------------
P = np.array([[-1, 0],
              [0,1],
              [1,0],
              [1, -.5],
              [1, -1],
              [1, -1.5],
              [1, -2],
              [2,-3],
              [3, -2]])
Bstep1 = Bezier(P.shape[0]-1, P.shape[1], T)
Bstep1.set_P(P)

t_sample = np.linspace(0,T,20)
A = Bstep1.get_A(t_sample, derivative=derivative)

curve = Bstep1.evaluate_in_time(np.linspace(0,T,50), derivative=derivative)
curve2 = (A@Bstep1.P).reshape((-1, P.shape[1]))

v_curve = (np.diff(curve.T)*(P.shape[0]-1)/T).T
plt.plot(v_curve[:,0], v_curve[:,1], label='Vel on difference')

plt.plot(curve[:,0], curve[:,1], label='Direct evaluation')

plt.plot(curve2[:,0], curve2[:,1], label='From A matrix')
plt.legend()
plt.scatter(P[:,0], P[:,1])
plt.show()