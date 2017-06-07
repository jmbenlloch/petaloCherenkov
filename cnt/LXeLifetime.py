"""
Parameterizes the scintillation lifetime of LXe
JJGC February 2016

Scintillation in LXe has two origins: 

1) exciton trapping
2) Recombination

In both cases, the excited molecular states are: singles (tau1 = 2.2 ns)
and triple (tau2 = 27 ns).

In the case of recombination, the intensity of the scintillation light can be parameterized
as:

Ir(t) = Ar1* Ir1(t) + Ar2 * Ir2(t)

where Ar1 = 0.44, Ar2 = 0.56

I1r(t)  = (1 + Tr^-1 * t)^-2 * (exp^{-t/tau1})/tau1
Tr = 15 ns, tau1 = 2.2 ns

Ir2(t) = (exp^{-t/tau2})/tau2 
tau2 = 27 ns

In the case of exciton trapping:

Ie(t) = Ae1* I1(t) + Ae2 * I2(t)

where Ae1 = 0.17, Ae2 = 0.83

I1e(t)  =  (exp^{-t/tau1})/tau1
 tau1 = 2.2 ns

Ie2(t) = (exp^{-t/tau2})/tau2 
tau2 = 27 ns

And: I(t) = 0.7 *Ir(t) + 0.3 Ie(t)
"""

from Util import *
from scipy.optimize import curve_fit


Tr = 15*ns # characteristic recombiantion time in LXe
tau1 = 2.2*ns
tau2 = 27*ns
Ar1 = 0.44
Ar2 = 0.22
Ae1 = 0.17
Ae2 = 0.83
  

def photon_generator(tau, taumin, taumax, size):
	x = np.linspace(taumin, taumax, size)
	y = np.random.exponential(scale=tau, size=size)
	return x,y
	


if __name__ == '__main__':
	nphotons = 30000
	tau=tau1
	taumin = 0*ns
	taumax = 30*ns

	x,y = photon_generator(tau, taumin, taumax, 10000)
	print("x = {0}, y ={1}".format(x,y))
	plt.plot(x, y)
	plt.show()

	fig = plt.figure()
	ax = fig.add_subplot(111)
	# the histogram of the data
	n, bins, patches = ax.hist(y, 50, normed=1, facecolor='green', alpha=0.75)
	ax.set_xlabel('t')
	ax.set_ylabel('Probability')
	ax.set_title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
	ax.set_xlim(0, 30)
	#ax.set_ylim(0, 0.03)
	ax.grid(True)

	plt.show()

	