import sys
sys.path.append('.')

from utils.pkg_func import *

from utils.Trajectory import Bezier



P = np.array([(0.09374999999999989, 0.15297619047619054),
         (0.549107142857143, 0.1648809523809524),
         (0.7083333333333335, 0.6142857142857144),
         (0.5282738095238095, 0.8940476190476193),
         (0.24404761904761907, 0.8776785714285716),
         (0.15327380952380942, 0.6321428571428573),
         (0.580357142857143, 0.08303571428571432),
         (0.8839285714285716, 0.28988095238095246)])

P = np.array([[-1, 0],
              [0,1],
              [1,0],
              [1, -.5],
              [1, -1],
              [1, -1.5],
              [1, -2],
              [2,-3],
              [3, -2]])


B = Bezier(degree=8, dim=2, T=4)

B.set_P(P)

ts = time.time()
curve = B.evaluate_in_time(np.linspace(0,4,50), derivative=0)
print(f'{time.time()-ts}')


plt.scatter(B.get_Pmat()[:,0], B.get_Pmat()[:,1])

plt.plot(curve[:,0], curve[:,1])
plt.axis('equal')
plt.show()