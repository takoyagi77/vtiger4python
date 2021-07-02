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
# if necessary
pip install -r requirement.txt
python vtiger_class.py
```

Optimization is finished, result is visualized by [matplotlib](https://matplotlib.org)


## Advanced usage

If you want to optimize more transfer functions or insert V-Tiger into another program, try the following code.

```python
import vtiger_class
...
# setup the V-Tiger
# The smallest argument setting
V = Vtiger(u00=u00, y00=y00, r00=r, r=r, ts=ts, th0=th0)
# The largest argument setting
V = Vtiger(u00=u00, y00=y00, r00=r, r=r, ts=ts, th0=th0, wST=0.02, OVr=2, GMr=3, PMr=20, optimize='PSO')
# Optimize the controller
th = V.VtigerPID()

# PID controller design
K = matlab.tf([th[2], th[0], th[1]], [0, 1, 0])
```

### Parameter

* u00, y00 : Input/Output data list(type : ndarray).
* r00, r   : Reference data list for feedback system(type : ndarray).
* ts       : Sampling time(type : float).
* th0      : Init controller parameters(type : list). In case of PID controller, first one is proportional gain, second one is integral gain and last one is differential gain(all type : float).
* wST      : Settling time threshold setting for evaluation step response(type : float).
* OVr      : Allowable overshoot amount\[%\](type : float or int)



# Note

注意点などがあれば書く


# License

'vtiger4python' is under [MIT license](https://en.wikipedia.org/wiki/MIT_License).