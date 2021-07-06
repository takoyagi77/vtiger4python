import math
import control
from control import matlab
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fmin
import cmath
import copy
import pyswarms as ps


class Vtiger:
	def __init__(self, **kwargs):
		if 'freq' in kwargs:
			kwargs = kwargs['freq']
		if 'ts' in kwargs:
			ts = kwargs['ts']
		else:
			ts = 1
		z = matlab.tf([1, 0], [0, 1], ts)
		p = (1 - 1 / z) / ts
		self.count = 0
		self.f = 0
		if 'y00' in kwargs:
			self.y0jw = self.fft4step(kwargs['y00'])[0]
		if 'u00' in kwargs:
			self.u0jw = self.fft4step(kwargs['u00'])[0]
		if 'r00' in kwargs:
			self.r0jw = self.fft4step(kwargs['r00'])[0]
		self.p = self.fft4tf(p, len(self.u0jw))
		if 'r' in kwargs:
			self.r = kwargs['r']
		if 'th0' in kwargs:
			self.th0 = kwargs['th0']
		if 'wST' in kwargs:
			self.wST = kwargs['wST']
		else:
			self.wST = 0.02
		if 'OVr' in kwargs:
			self.OVr = kwargs['OVr']
		else:
			self.OVr = 2
		if 'GMr' in kwargs:
			self.GMr = kwargs['GMr']
		else:
			self.GMr = 3
		if 'PMr' in kwargs:
			self.PMr = kwargs['PMr']
		else:
			self.PMr = 20
		if 'optimize' in kwargs:
			temp = kwargs['optimize']
			try:
				for i in temp.keys():
					self.optimization = i
				self.optimizer = temp[self.optimization]
			except Exception:
				self.optimization = 'fmin'
		else:
			self.optimization = 'fmin'

	def fft4step(self, y00):
		'''
		Convert steady-state signals, such as step responses, into periodic signals.

		Exsamples
		---------
		>>> import numpy as np\n
		>>> from matplotlib import pyplot as plt\n
		>>> import pyfftw\n
		>>> u00 = np.ones([100, 1])\n
		>>> u00[0] = 0 # IMPORTANT\n
		>>> u0jw, w_ts = fft4step(u00)\n
		>>> plt.plot(w_ts, u0jw); plt.show()\n

		Parameters
		----------
		y00 : np.array
			Step response like date. y[0] and y[-1] must be stational.
		
		Returns
		-------
		y0jw : complex
			The Discrete Fourier Transform (DFT) of y00 using a Fast Fourier Transform (FFT) algorithm.
		w_ts : np.array
			Vector of the frequency domain.
		
		See Also
		--------
		fft4step.m : The MATLAB version.
		'''
		if (np.linalg.norm(y00 - y00[0] * np.ones([len(y00), 1]), ord=2)) == 0:
			print('Error in fft4step: Data in which all elements are the same cannot be cycled.')
			print('\tSet initial value different from others. Ex. u00(0)=0, r00(0)=0;')
			return
		n = min(20, round(len(y00) / 10 + 1))
		y00ave = sum(y00[-n:]) / n + y00[0]
		y0 = np.array(y00.tolist() + list(map(lambda x: x + y00ave, -1 * y00)))
		out = np.asarray(y0)
		y0jw = np.fft.fft(out.reshape(-1))
		N = len(y0jw)
		w_ts = np.append(1e-10, 2 * math.pi / (N / np.arange(1, N, 1)))
		return y0jw, w_ts

	def fft4tf(self, G, N):
		'''
		Convert the transfer function to a frequency signal using a Fast Fourier Transform (FFT) algorithm.

		Exsamples
		---------
		>>> import numpy as np\n
		>>> import control\n
		>>> from control import matlab\n
		>>> from matplotlib import pyplot as plt\n
		>>> import pyfftw\n
		>>> u = np.concatenate([np.ones([10, 1]), np.zeros([10, 1])], 0)\n
		>>> a = pyfftw.empty_aligned(len(u), dtype='complex128'); b = pyfftw.empty_aligned(len(u), dtype='complex128'); c = pyfftw.empty_aligned(len(u), dtype='complex128')\n
		>>> fft_object = pyfftw.FFTW(a, b)\n
		>>> c[:] = list(u)\n
		>>> ujw = fft_object(c)\n
		>>> s = matlab.tf('s')\n
		>>> G = matlab.c2d(1 / (s + 1), 1)\n
		>>> h = fft4tf(G, len(u))\n
		>>> yf = h * ujw\n
		>>> y = np.fft.ifft(yf)\n
		>>> plt.plot(np.abs(y)); plt.show()\n

		Parameters
		----------
		G : control.xferfcn.TransferFunction(Discreate)
			The transfer function you want to convert into the frequency domain.
		N : int
			Length of input signal.
		
		Returns
		-------
		h : np.array
			The frequency response of G(z).
		
		See Also
		--------
		fft4tf.m : The MATLAB version.
		'''
		w = np.append(1e-10, 2 * math.pi / (N / np.arange(1, N / 2 + 1, 1) * G.dt))
		g, p, _ = matlab.bode(G, w.astype(np.float64), plot=0)  # g:Gain p:Phase[rad]
		gl, pl, _ = control.freqresp(G, [w[-1]])
		g = np.append(g, gl); p = np.append(p, pl)
		h = g * np.exp(1j * p)
		h = np.array(h.tolist() + np.flipud(h[1:len(h) - 1].conjugate()).tolist())
		return h

	def freq2yu(self, *args): # args[0]=freq,args[1]=K,args[2]=ts
		'''
		Convert the Frequency response in a plant date with Controller to a Time response.

		Parameters
		---------
		freq : dict
			MUST key		1. y0jw : Frequency response of cyclic output data of the plant.
						2. u0jw : Frequency response of cyclic input data of the plant.
			NOT MUST key		1. r0jw : Frequency response of cyclic reference signal into the feedback system with K.
						2. d0jw : Frequency response of cyclic output disturbance into the feedback system with K.
						3. du0jw : Frequency response of cyclic input disturbance into the feedback system with K.
		K : control.xferfcn.TransferFunction(Discreate) or list
			TransferFunction : Controller to be evaluated.
			list : PID gains. K[0] = kp, K[1] = ki, K[2] = kd
		ts : float
			Sampling time [s]
		
		Returns
		--------
		y : np.array
			Predicted time response of the closed loop output in which K is inserted.
		u : np.array
			Predicted time response of the closed loop system input in which K is inserted.
		r2 : np.array
			r(t) in which set rjw = 0 when u0jw = 0.
		
		See Also
		---------
		freq2yu.m : The MATLAB version.
		'''
		if hasattr(self, 'u0jw'):
			u0jw = self.u0jw
		else:
			u0jw = args[0]['u0jw']
		if hasattr(self, 'y0jw'):
			y0jw = self.y0jw
		else:
			y0jw = args[0]['y0jw']

		if len(args) != 3:
			K = args[1]
			pass
		else:
			K = args[1] + 0 * matlab.tf([1, 0], [0, 1], args[2])
		
		if isinstance(K[0], float) and isinstance(K[1], float) and isinstance(K[2], float):
			Kjw = K[0] + K[1] / self.p + K[2] * self.p
			invKjw = 1 / Kjw
		else:
			invKjw = self.fft4tf(1 / K, len(u0jw))
		
		if hasattr(self, 'r0jw'):
			r0jw = self.r0jw
		elif 'r0jw' in args[0]:
			r0jw = args[0]['r0jw']
		else:
			r0jw = y0jw * 0
		if hasattr(self, 'd0jw'):
			d0jw = self.d0jw
		elif 'd0jw' in args[0]:
			d0jw = args[0]['d0jw']
		else:
			d0jw = y0jw * 0
		if hasattr(self, 'du0jw'):
			du0jw = self.du0jw
		elif 'du0jw' in args[0]:
			du0jw = args[0]['du0jw']
		else:
			du0jw = u0jw * 0
		
		r1jw = invKjw * u0jw + y0jw
		r1jw[r1jw == 0] = 1e-10
		
		
		yjw = (r0jw * y0jw + invKjw * (d0jw * u0jw + du0jw * y0jw)) / r1jw
		ujw = (r0jw - yjw) / invKjw
		
		yjw[r1jw == 0] = 0
		ujw[r1jw == 0] = 0
		ujw[invKjw == 0] = 0
		
		y = np.fft.ifft(yjw); y = y[0:round(len(y) / 2)]
		u = np.fft.ifft(ujw); u = u[0:round(len(u) / 2)]
		r2 = np.fft.ifft(r0jw); r2 = r2[0:round(len(r2) / 2)]

		if np.linalg.norm(y[:].imag) != 0:
			tmp = np.linalg.norm(y.imag) / np.linalg.norm(y.real)
			if tmp > 1e-5:
				print('Warning in V-Tiger: An imaginary part(im / re = 1e-5) exists at y. imag / real = ', tmp)
			y = y.real
			u = u.real
			r2 = r2.real
			if not hasattr(self, 'r'):
				self.r = r2

		r00 = np.fft.ifft(r0jw); r00 = r00[0:round(len(r00) / 2)]
		if np.linalg.norm(r00[:].real - r2) and self.count != 1:
			self.count = 1
			print('Error in V-Tiger: Removed the frequency component of r that has a it not in u.')
		return y, u, r2

	def stepinfo(self, y, T, *args, **kwargs):
		"""
		Computes the step-response characteristics from an array of step-response data.\n
		This is alpha version.

		Parameters
		-----------
		y : np.array
			An array of step-response data.
		T : np.array
			Corresponding time vector with y.
		args : float
			args[0]( = yfinal) : The steady-state value of y.
		kwargs : dict
			SettlingTimeThreshold : The threshold in the definition of settling time. The default value is 0.02.
			RiseTimeLimits : NOT already. The lower and upper threshold in the definition of rise time. The default range is 10% ~ 90%.

		Returns
		-------
		si : dict
			Peak : Peak absolute value of y.
			PeakTime : Time at which the peak value occur.
			Overshoot : Percentage overshoot, relative to yfinal.
			SettlingTime : Time it takes for the error |y - yfinal| between the response y and the steady-state response yfinal to fall to within 2% of yfinal.
			RiseTime : NOT already. Time it takes for the response to rise from 10% to 90% of the steady-state response.
		
		See Also
		---------
		stepinfo.m : The MATLAB version.
		"""

		if 'SettlingTimeThreshold' in kwargs:
			SettlingTimeThreshold = kwargs['SettlingTimeThreshold']
		else:
			SettlingTimeThreshold = 0.02
		# Calculate steady-state value
		ylast = sum(y[-20:]) / len(y[-20:])
		if args == ():
			yfin = ylast
		else:
			yfin = args[0]
		ts = max(T) / len(T)
		Peak = max(y.real)
		PeakTime = np.argmax(y) * ts
		if Peak - yfin > 0:
			Overshoot = (max(y.real) - yfin) / yfin * 100
		else:
			Overshoot = (max(y.real) - yfin) / yfin * 0
		
		if (ylast <= yfin * (1.0 + SettlingTimeThreshold)) and (ylast >= yfin * (1.0 - SettlingTimeThreshold)):
			# 定常値が収束してる場合
			tmp = copy.copy(y)
			tmp[yfin * (1 + SettlingTimeThreshold) < tmp] = 0
			tmp[yfin * (1 - SettlingTimeThreshold) > tmp] = 0
			tmp[tmp != 0] = 1
			tmp = np.flipud(tmp[:])
			num = 0
			for i in tmp:
				if i == 0:
					break
				else:
					num += 1
			SettlingTime = T[-num]
		else:
			# 定常値が収束していない場合
			SettlingTime = np.inf
		
		si = {'Peak': Peak, 'PeakTime': PeakTime, 'Overshoot': Overshoot, 'SettlingTime': SettlingTime}

		return si

	def constraints(self, th):
		global gl
		if hasattr(self, 'r'):
			r = self.r
		else:
			r = gl['r']
		p2 = self.p2

		if math.isnan(th[0]) or math.isnan(th[1]) or math.isnan(th[2]):
			c = 1e99*np.ones([4,1])
		else:
			Kp = th[0]
			Ki = th[1]
			Kd = th[2]
			c = []
			Kjw = Kp + Ki / p2 + Kd * p2
			GKjw = self.Gjw * Kjw
			GKjw = GKjw[0:round(len(GKjw) / 2)]
			phi = np.array(list(map((lambda x: math.degrees(cmath.phase(x))), GKjw)))
			Gm, Pm, _, _ = matlab.margin(abs(GKjw), phi, self.w[0:round(len(self.w) / 2)])
			c.append(self.GMr - 20 * math.log10(Gm))
			c.append(self.PMr - Pm)
			if math.isnan(Kd + Kp + Ki):
				c.append(np.inf)
				c.append(np.inf)
			else:
				c.append(np.max(matlab.pole(matlab.tf([0,0,1],[Kd,Kp,Ki])).real) + 0.001)
				y, _, _ = self.freq2yu(gl, th)
				k = np.arange(1, len(y) + 1, 1)
				si = self.stepinfo(y, k, sum(r[len(r) - 20:len(r)]) / 20)
				c.append((si['Overshoot'] - self.OVr).tolist()[0])
		c = list(map((lambda x: float(x)), c))
		ceq=[]
		return c, ceq

	def J_cost(self, th, *f):
		global gl
		if hasattr(self, 'r'):
			r = self.r
		else:
			r = gl['r']
		if math.isnan(th[0]) or math.isnan(th[1]) or math.isnan(th[2]):
			J = 1e99
		else:
			y, _, _ = self.freq2yu(gl, th)
			k = np.arange(1, len(y) + 1, 1)
			y[y < 0] = 0
			si = self.stepinfo(y[0:len(y) - 2], k[0:len(k) - 2], sum(r[len(r) - 20:len(r)]) / 20, SettlingTimeThreshold=self.wST)
			J = si['SettlingTime']
		if f:
			plt.plot(np.ones(len(y), 1) * sum(r[-20:]) / 20 * [1 + self.wST, 1 - self.wST, 1 + self.OVr / 100])
			plt.show()
			return
		if math.isnan(J):
			J = 1e99
		print('J=' + str(J))
		return J

	def J_cost2(self, th):
		if self.optimization == 'PSO':
			J = []
			for i in th:
				c, ceq = self.constraints(i)
				ctmp = c + ceq
				ctmp.append(0)
				Ji = self.J_cost(i) + 1e10*max(ctmp)
				if self.f:
					print(Ji)
				J.append(Ji)
		else:
			c, ceq = self.constraints(th)
			ctmp = c + ceq
			ctmp.append(0)
			J = self.J_cost(th) + 1e10*max(ctmp)
			if self.f:
				print(J)
		return J

	def vtigerPID(self, **kwargs):
		global gl
		if 'freq' in kwargs:
			gl = kwargs['freq']
		else:
			gl = ''
		if hasattr(self, 'y0jw'):
			y0jw = self.y0jw
		else:
			y0jw = gl['y0jw']
		if hasattr(self, 'u0jw'):
			u0jw = self.u0jw
		else:
			u0jw = gl['u0jw']
		if hasattr(self, 'th0'):
			th0 = self.th0
		else:
			th0 = gl['th0']
		G = np.divide(y0jw, u0jw, where=abs(u0jw) != 0)
		N = len(u0jw)
		ts = 1
		self.w = np.append(1e-10, 2 * math.pi / (N / np.arange(1, N, 1) * ts))
		ignor = np.divide(sum(abs(u0jw)), N) * 1e-10
		G = G[~(abs(u0jw) < ignor)]
		self.w = self.w[~(abs(u0jw) < ignor)]
		if 'freq' in kwargs:
			self.p2 = gl['p']
		else:
			self.p2 = self.p
		self.p2 = self.p2[~(abs(u0jw) < ignor)]
		self.Gjw = G

		if 'f' in kwargs:
			if kwargs['f'] == 0:
				pass
			else:
				self.J_cost(th0, f)
				return
		else:
			pass
		
		if self.optimization == 'PSO':
			cost, th = self.optimizer.optimize(self.J_cost2, iters=500)
		else:
			th = fmin(self.J_cost2, th0)
		
		return th



if __name__ == '__main__':
	s = matlab.tf('s')
	count = 0
	ts = 0.01
	num = [0, 0, 5]
	den = [0.01, 0.2, 10]
	Gs = matlab.tf(num, den)

	num1, den1 = matlab.pade(0.1, 5)
	delay = matlab.tf(num1, den1)
	G = matlab.c2d(Gs * delay, ts)
	z = matlab.tf([1, 0], [0, 1], ts)
	p = (1 - 1 / z) / ts

	Ku, Pm, Wu, Wcp = matlab.margin(G)
	Tu = 1 / (Wu / 2 / math.pi)
	kp0 = 0.6 * Ku; ki0 = kp0 / (0.5 * Tu); kd0 = kp0 * 0.125 * Tu
	K0 = kp0 + ki0 / p + kd0 * p
	th0 = [kp0, ki0, kd0]
	th1 = copy.copy(th0)


	N = 3000
	u00 = np.ones([N, 1])
	u00[0] = 0
	t = np.arange(0, N * ts, ts)
	(y00, t00, a) = matlab.lsim(Gs, U=u00, T=t)
	plt.figure()
	plt.plot(t00, y00)
	r = u00
	y00 = np.array(y00)
	y00 = y00.reshape((len(y00), 1))

	options = {'c1': 0.5, 'c2': 0.3, 'w': 0.9}
	optimizer = ps.single.GlobalBestPSO(n_particles=20, dimensions=3, options=options)

	V = Vtiger(y00=y00, u00=u00, r00=r, r=r, ts=ts, th0=th0, optimize={'PSO': optimizer})
	th = V.vtigerPID()
	print(th)

	K = th[0] + th[1] / s + th[2] * s
	K0 = th1[0] + th1[1] / s + th1[2] * s
	G = Gs * K / (1 + Gs * K)
	G0 = Gs * K0 / (1 + Gs * K0)
	N = 3000
	u = np.ones([N, 1])
	u[0] = 0
	t = np.arange(0, N * ts, ts)
	(y, t, _) = matlab.lsim(G, U=u, T=t)
	(y0, t0, _) = matlab.lsim(G0, U=u, T=t)
	plt.figure()
	plt.plot(t, y, label='vtiger')
	plt.plot(t0, y0, label='init')
	plt.legend()
	plt.show()
	plt.close('all')
