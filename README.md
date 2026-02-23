# Meiotic-Exit-Code
Code implementation of the mathematical model describing meiotic exit during meiosis II in budding yeast.

This project analyzes and implements a minimal model of the **meiosis process in budding yeast**, a type of cell division essential for sexual reproduction.

This repository contains numerical implementations of the model presented in the associated manuscript. The system of nonlinear ordinary differential equations (ODEs) was simulated in **Python** and independently corroborated using **Mathematica (Wolfram Research, Version 13.3)** through the `NDSolve` function with a stiff solver.

 For the Python simulations we used interactive visualizations created using **Bokeh** and **Panel**. The plotting and modeling style is inspired by [Caltech's Bi/BE 150 course](http://be150.caltech.edu/2020/content/lessons/00_setting_up_python_computing_environment.html).
> It is **strongly recommended** to run the code in a **Jupyter Notebook or JupyterLab environment** to fully interact with the plots.


## Installation

Make sure you have Python 3.7+ installed. Use the package maneger pip to install:

```bash
pip install numpy
pip install scipy
pip install bokeh
pip install colorcet
pip install panel
pip install holoviews
pip install sympy
                 
