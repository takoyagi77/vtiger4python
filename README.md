# vtiger4python

'vtiger4python' allows you to design controllers in Python using the '[V-Tiger](https://github.com/kosaka3/vtiger_matlab)(This is MATLAB program.)' methodology.



# Features

* A controller can be designed from ONLY a set of step responses which are input and output datas from control target.
* Overshoot and settling time can be set as constraint conditions.
* No reference model required for optimization.

# Requirement

* Python 3.8.6 or above
* control 0.8.4
* numpy 1.20.1
* matplotlib 3.3.4
* scipy 1.6.1
* pyswarms 1.3.0

# Installation

vtiger_class.py is not need aditional install.
Please install above packages by the bottom command.

```bash
pip install -r requirement.txt
```

# Usage

'vtiger_class.py' can not only be excuted, but also declared as a 'class'.

## Basic usage

In the default settings, the following transfer function is optimized by [fmin](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fmin.html).

<img src='https://latex.codecogs.com/png.latex?G\left(s\right)&space;=&space;{\frac{5}{0.01s^2&plus;0.2s&plus;10}}' />

Plese type following command.


```bash
git clone https://github.com/takoyagi77/vtiger4python vtiger_demo_python
cd vtiger_demo_python
pip install -r requirement.txt # <- if necessary
python vtiger_class.py
```

Optimization is finished, result is visualized by [matplotlib](https://matplotlib.org)


## Advanced usage

If you want to optimize more transfer functions or insert V-Tiger into another program, try the following code.

```python
from vtiger_class import Vtiger
...
# setup the Particle Swarm Optimization
options = {'c1': 0.5, 'c2': 0.3, 'w': 0.9}
optimizer = ps.single.GlobalBestPSO(n_particles=20, dimensions=3, options=options)

# setup the V-Tiger
# The smallest argument setting
V = Vtiger(u00=u00, y00=y00, r00=r, r=r, ts=ts, th0=th0)
# The largest argument setting
V = Vtiger(u00=u00, y00=y00, r00=r, r=r, ts=ts, th0=th0, wST=0.02, OVr=2, GMr=3, PMr=20, optimize={'PSO': optimizer, 'iters':500})
# Optimize the controller
th = V.VtigerPID()

# Calculated final result
y, u, r = V.freq2yu({0: 0}, th)
```

### Parameter

* u00, y00 : Input/Output data list(type : ndarray).
* r00, r   : Reference data list for feedback system(type : ndarray).
* ts       : Sampling time \[s\](type : float).
* th0      : Init controller parameters(type : list). In case of PID controller, first one is proportional gain, second one is integral gain and last one is differential gain(all type : float).
* wST      : Settling time threshold setting for evaluation step response(type : float).
* OVr      : Allowable overshoot amount \[%\](type : float or int).
* GMr      : Gain margin \[dB\](type : float or int).
* PMr      : Phase margin \[deg\](type : float or int).
* optimize : Optimization method. Currently, it have 'fmin' and 'PSO'. If you want to optimize by PSO, attach the PSO optimizer and iters number(type : dict or str).


## Including method

* fft4step
* fft4tf
* stepinfo
* vtigerPID

# Note

I don't test environments under Linux and Mac.

# License

'vtiger4python' is under [MIT license](https://en.wikipedia.org/wiki/MIT_License).