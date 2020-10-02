# Some code done in C++ to experiment
This sections contains various different subjects.  

1. [Economic and population model based on Malthus theory](https://github.com/mintDan/Cppcode#economic-and-population-model-based-on-malthus-theory)
2. [RK4 on system of equations with change of variables](https://github.com/mintDan/Cppcode#rk4-on-system-of-equations-with-change-of-variables)
3. [Adaptive RK4 for Baseball](https://github.com/mintDan/Cppcode#adaptive-rk4-for-baseball)


## Economic and population model based on Malthus theory
Thomas Malthus reasoned that technological advancement does not necessarily lead to improved quality of life, as measured by proxy GDP/capita.
The reason is that a new technology might increase output, GDP, but, with increased output follows an increase in population, and thus GDP/capita can be constant.
Any one person in this model will thus not feel major improvements in quality of life. The break of this model for a country, that is, when GDP/capita starts to significantly increase,
is called the demographic transition, which is when a country becomes a more modern society.

To describe the model we use Y = GDP, L = population, y = GDP/population, n is a measure for birthrates, d is a measure for deathrates(assumed constant here),
A is a measure for productivity in terms of technological capabilities.

![n](https://github.com/mintDan/Cppcode/blob/master/figs/n.png)

The rate of change in population w.r.t time is given below. This can be used to forecast the model with different outcomes for different parameters.

![L](https://github.com/mintDan/Cppcode/blob/master/figs/L.png)

The script predicts the evolution of GDP and population L and writes the data to files, which are plotted with Python.
As can be seen, even though GDP increases, the population also increases, thus GDP/capita approaches an asymptote.

![Malthus](https://github.com/mintDan/Cppcode/blob/master/figs/Malthus.png)


## RK4 on system of equations with change of variables
Just an experiment to try RungeKutta4 on systems of equations expressed in cartesian and polar coordinates. RK4 and other methods are good because they can reduce systematic errors of gradients among other things.
Also as an alternative to reducing the timestep dt, since in the end, continually reducing dt would lead to float errors and misrepresentation, so in the end, one would go to better integration techniques anyway.

The system in cartesian and polar coordinates is seen below. Note these are all functions of t and derivatives w.r.t variable t.

![xysystem](https://github.com/mintDan/Cppcode/blob/master/figs/xysystem.png)

Solving them shows that RK4 gives pretty much equilevant answers when expressed.

![pp](https://github.com/mintDan/Cppcode/blob/master/figs/PP.png)

## Adaptive RK4 for Baseball
The forces modelled here are Magnus force, gravity and drag force. RK4 has been simply modified to include an adaptive step size.  
The Magnus force can bend the path of the baseball in funny ways.

![BB](https://github.com/mintDan/Cppcode/blob/master/figs/Slider.png)


