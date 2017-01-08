#Economic/population model based on Malthus' theory
Thomas Malthus reasoned that technological advancement does not necessarily lead to improved quality of life, as measured by proxy GDP/capita.
The reason is that a new technology might increase output, GDP, but, with increased output follows an increase in population, and thus GDP/capita can be constant.
Any one person in this model will thus not feel major improvements in quality of life. The break of this model for a country, that is, when GDP/capita starts to significantly increase,
is called the demographic transition, which is when a country becomes a more modern society.

To describe the model we use Y = GDP, L = population, y = GDP/population, n is a measure for birthrates, d is a measure for deathrates(assumed constant here),
A is a measure for productivity in terms of technological capabilities.

![n](https://github.com/mintDan/Cppcode/blob/master/figs/n.png

The rate of change in population w.r.t time is given below. This can be used to forecast the model with different outcomes for different parameters.

![L](https://github.com/mintDan/Cppcode/blob/master/figs/L.png

The script predicts the evolution of GDP and population L and writes the data to files, which are plotted with Python.
As can be seen, even though GDP increases, the population also increases, thus GDP/capita approaches an asymptote.

![Malthus](https://github.com/mintDan/Cppcode/blob/master/figs/Malthus.png
 



