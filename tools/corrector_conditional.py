import numpy as np
import pylab as plt
import math
import sys

bones = bool(int(sys.argv[1]))
if bones:
	rho, lamda, mu = 1.904,  5.891, 0.982
else:
	rho, lamda, mu = 0.916,  1.415, 0.236
K = lamda + 2 * mu
c1 = math.sqrt(K / rho)
c2 = math.sqrt(mu / rho)

def BOmega(x):
	y = np.zeros(x.shape)
	z = np.sqrt(1 - x**2 - y**2)
	zero = np.zeros(x.shape)
	tmp = np.array([
		[    K * x * c2, mu * y * c1, mu * z * c1],
		[lamda * y * c2, mu * x * c1,        zero],
		[lamda * z * c2,        zero, mu * x * c1]
	]) / (c1 * c2)
	return tmp.transpose(2, 0, 1)

def acsBOmega(x):
	return np.array([[x]]).transpose()

def plotBOmega(fBOmega, color, name):
	step = 0.001
	x = np.arange(step / 4, 1 + step / 4, step)
	BO = fBOmega(x)
	abs_det = np.abs(np.linalg.det(BO))
	max_det = abs_det[-1]
	plt.plot(x, abs_det / max_det, color + '-', label = '$|det| / |det|_{max}$ - ' + name)
	plt.plot(x, np.linalg.cond(BO) / 1000, color + ':', label = '$\\mu / 1000$ - ' + name)
	#plt.plot(x, 1 / np.linalg.norm(BO, axis = (1, 2), ord = -2) / 1000, color + ':', label = '$norm_{-2} / 1000$ - ' + name)


plotBOmega(BOmega, 'r', 'elastic-3d')
plotBOmega(acsBOmega, 'b', 'acoustic')

plt.xlabel('$(\\vec{n} \\cdot \\vec{l})$', fontsize="xx-large")
plt.grid()
plt.xlim(xmin = 0)
plt.xlim(xmax = 1)
plt.ylim(ymin = 0)
plt.ylim(ymax = 1)
plt.legend(fontsize="xx-large")
plt.show()

