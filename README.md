# Casus 1: Dynamic Modeling (ODEs)

**Authors**: Joris Deelstra (j.deelstra.2@st.hanze.nl), Emiel Bosma (e.l.bosma@st.hanze.nl)

**Professor**: Tsjerk Wagenaar (t.a.wassenaar@pl.hanze.nl)

_Hanze University of Applied Sciences Groningen_

## Table of contents

1. [Goal of this repository](#1-goal-of-this-repository);
2. [Content of this repository](#2-content-of-this-repository);
3. [Background information](#3-background-information);
4. [Installation & Requirements](#4-installation--requirements);
5. [Models](#5-models);
6. [Methods](#6-methods);
7. [Optimization](#7-optimization);
8. [Model-Comparing](#8-model-comparing);
9. [References](#9-references).

## 1. Goal of this repository

This repository is for the school assignment Bio-Informatics - modeling cancer - Casus A.

The objectives are:
* Making a python class;
* Adding an [SIR-model](https://en.wikipedia.org/wiki/Compartmental_models_(epidemiology)) for infectious diseases;
* Adding ordinary differential equations (ODE) solving methods for cancer growth;
* Optimising the ODE solving methods by means of the [Mean Squared Error](https://en.wikipedia.org/wiki/Mean_squared_error) (MSE) and the [Monte Carlo method](https://en.wikipedia.org/wiki/Monte_Carlo_method);
* Comparing the ODE solving methods by means of the [Bayes Information Criterion](https://en.wikipedia.org/wiki/Bayesian_information_criterion) (BIC) and the [Akaike Information Criterion](https://en.wikipedia.org/wiki/Akaike_information_criterion) (AIC), with also the corrected method (AICc).


## 2. Content of this repository

```
C:.
|   .gitignore                          # What git ignores, e.g. pycharm .idea/ or visual studie code vscode/
|   main.py                             # The python class for dynamic moddeling of ODEs
|   README.md                           # This file, an explanation of the project and installation guide
|   requirements.txt                    # The required python libraries for this project
|   Usecase_example.ipynb               # An example of how to handle the class
|             
+---generated_data                      # Example data
|       peerj-07-6983-s001.csv          # Cancer growth data from an research on mouse
|       sir_dummy_dataset.csv           # Generated infection data with outliers
|       tumorgroei_dummy_dataset.csv    # Generated cancer growth data with outliers
|       
\---logboeken                           # The logs of this project, written in Dutch
    |   Logboek_Emiel.md
    |   Logboek_Joris.md
    |   
    \---images                          # Images used in the logs
            Euler_plot_sirdata_infected_per_region.png
            optimalization_linear.jpg
            optimization_exp_1.png
```

## 3. Background information

This casus is inspired by the student Kenneth Chan of the Pine Crest School, Claremont Mckenna College. They wrote a short paper called "[Treatment Optimization for Tumor Growth by Ordinary Differential Equations](https://doi.org/10.47611/jsrhs.v12i4.5202)".
In this paper applies the author among other things several different models of tumor growth to various existing datasets in order to identify the best-fitting model.

Mathematical models, i.e. ordinary differential equations (ODEs), can be used to model tumor growth and study treatment of cancer [[1]](#9-references). 
This makes it essential for bioinformatics to understand and be able to apply such models and equations on patient data. 

## 4. Installation & Requirements

How to install this repository correctly:

1. To install this git repository, go to the folder you want to pase this git repo in;
2. Open the console in that folder;
3. Run the console command `git clone https://github.com/emielbosma/case1_modelleren.git`;

This project requires python version 3 or higher.
If not done before, download python from the [official website](https://www.python.org/downloads/), or if you use microsoft go to microsoft store and search for python.

This project uses the python libraries matplotlib, numpy and pandas.

4. Run in the same console as before the console command `pip install -r requirements.txt`.

## 5. Models

The repository implements several mathematical models to 
describe tumor growth and infectious disease spread:

Parameters:

(V): Tumor volume

(t): Time

(c): Growth-rate constant

(d): Growth exponent

(Vmax): Maximum tumor volume

### a. SIR-model

The SIR-model describes the spread of infectious diseases by 
dividing a population into three groups: susceptible (S), infected (I) and 
recovered (R).

### b. Linear growth

Constant growth rate, independant of the current tumor volume.

$$
\frac{dV}{dt}=c
$$

### c. Exponential growth

Assumes that growth rate is proportional to the current tumor volume.

$$
\frac{dV}{dt}=c\cdot V
$$

### d. Mendelsohn growth

Describes tumor growth using a power-law relationship.

$$
\frac{dV}{dt}=c\cdot V^d
$$

### e. Exponential flattening growth

Assumes that tumor growth slows down as the tumor volume approaches 
a maximum.

$$
\frac{dV}{dt}=c\cdot(V_{\max}-V)
$$

### f. Logistic growth 

A logistic function is a common S-shaped curve with the following equation:

$$ f(x) = \frac{L}{1 + e^{-k(x - x_0)}} $$

When used for biology, i.e. for population growth, it typically gets rebuilt as:

$$ \frac{dN}{dt} = rN \cdot \left(1 - \frac{N}{K}\right) $$

Where `N` is the population size, or in our case volume, `t` is time, `r` is the growth rate and `K` is the carrying capacity. [[7]](#9-references)

We've rewritten it as:

$$ \frac{dV}{dt} = c \cdot V \cdot (V_{max} - V) $$

Where `V` is the volume of the tumor, `Vmax` is the maximum volume (`K`) and `c = r / Vmax`. The tumor first grows exponentially, then slows down and levels off at `Vmax`. The growth is fastest at `V = Vmax / 2`.

### g. Montroll

The Montroll model is a generalized logistic model with an extra shape parameter `d`:

$$ \frac{dV}{dt} = c \cdot V \cdot (V_{max}^d - V^d) $$

With `d = 1` this is exactly the logistic model. The value of `d` determines at which volume the growth is fastest: with a small `d` that point lies lower, with a large `d` it lies closer to `Vmax`. This makes the S-curve asymmetric. [8](#9-references)

### h. Allee effect

The Allee effect describes populations that need a minimum size to be able to grow. In biology it is often written as:

$$ \frac{dN}{dt} = rN \cdot \left(\frac{N}{A} - 1\right) \cdot \left(1 - \frac{N}{K}\right) $$ [[9]](#9-references)

We've rewritten it as:

$$ \frac{dV}{dt} = c \cdot (V - V_{min}) \cdot (V_{max} - V) $$

Where `Vmin` is the threshold volume. Below `Vmin` the tumor shrinks, above `Vmin` it grows towards `Vmax`. For tumors this can represent cells needing each other (cooperation) to survive.

### i. Lineair gelimiteerd

In this model the growth starts exponentially and then turns into linear growth:

$$ \frac{dV}{dt} = c \cdot \frac{V}{V + d} $$

For small `V` this is about `(c/d) · V` (exponential); for large `V` it is about `c` (linear). `d` is the volume at which the growth rate is half of the maximum `c`. The form is similar to Michaelis-Menten kinetics.

### j. Oppervlakte-gelimiteerd

Here the growth of a large tumor is limited by its surface, because nutrients and oxygen enter through the outside:

$$ \frac{dV}{dt} = c \cdot \frac{V}{(V + d)^{1/3}} $$

For small `V` the growth is exponential; for large `V` it is about `c · V^(2/3)`, which is proportional to the surface area of a sphere. `d` determines where the switch happens.

### k. Von Bertalanffy

The Von Bertalanffy model sees growth as the difference between building up and breaking down:

$$ \frac{dV}{dt} = c \cdot V^{2/3} - d \cdot V $$

The first term is the build-up (anabolism), which goes through the surface (`V^(2/3)`). The second term is the breakdown (catabolism), which is proportional to the volume. The tumor stops growing when both are equal, at:

$$ V_{max} = \left(\frac{c}{d}\right)^3 $$ [[10]](#9-references)

### l. Gompertz

The Gompertz model is one of the most used models for tumor growth. The relative growth rate decreases as the tumor gets bigger:

$$ \frac{dV}{dt} = c \cdot V \cdot \ln\left(\frac{V_{max}}{V}\right) $$

With the solution:

$$ V(t) = V_{max} \cdot e^{\ln(V_0 / V_{max}) \cdot e^{-c t}} $$

Like the logistic model this gives an S-curve, but it is asymmetric: the growth is fastest at `V = Vmax / e` (about 37% of `Vmax`) instead of at the halfway point. [[11]](#9-references)

## 6. Methods

### a. Euler method

Euler's method calculates the derivative at the start of an step, and makes on basis of that the next step.

The Euler method is:

$$ y_n+1 = hf(t_n,y_n) $$ [[4]](#9-references)

We've written it as:

```python
def euler_method(data, y_name, h=1):
    f, ts, ys = prepare(data, y_name, h)
    for i in range(len(ts) - 1):
        ys[i + 1] = ys[i] + h * f(ts[i], ys[i])
```

Were `f` is one of the mathematical models described in `5. Models`. \
`ts` is the (time)frame. \
`h` is how large the step is over `ts`. \
`ys` is in our case how big the formula thinks the tumor gets.

### b. Heun's method

Euler's method is used as the foundation for Heun's method. 
Instead of calculating the derivative based on one point, it uses the following mathematical formulas to calculate it:

$$ Slope_{\text{left}} = f(x_i,y_i) $$
$$ Slope_{\text{right}} = f(x_i + h, y_i + hf(x_i, y_i)) $$
$$ Slope_{\text{ideal}} = 1/2 (Slope_left + Slope_right) $$

[[5]](#9-references)

We've rewritten it as:

```python
def heun_method(data, y_name, h=1):
    f, ts, ys = prepare(data, y_name, h)
    for i in range(len(ts) - 1):
        k1 = f(ts[i], ys[i])
        k2 = f(ts[i] + h, ys[i] + h * k1)
        ys[i + 1] = ys[i] + h / 2 * (k1 + k2)
    return ts, ys
```

### c. Runge-Kutta methods

With the Runge-Kutta method we mean the Fourt-order Runge-Kutta method (RK4).
Instead of one or two points to calculate the deritative from, RK4 uses as the nume suggest 4 points.
The formulas are as follows:

$$ k_1 = f(t_n, y_n) $$
$$ k_2 = f(t_n + h/2, y_n + k_1 \cdot h/2) $$
$$ k_3 = f(t_n + h/2, y_n + k_2 \cdot h/2) $$
$$ k_4 = f(t_n + h, y_n + hk3) $$

[[6]](#9-references)

We've rewritten it as:

```python
def rungekutta_method(data, y_name, h=1):
    f, ts, ys = prepare(data, y_name, h)
    for i in range(len(ts) - 1):
        k1 = f(ts[i], ys[i])
        k2 = f(ts[i] + h / 2, ys[i] + h / 2 * k1)
        k3 = f(ts[i] + h / 2, ys[i] + h / 2 * k2)
        k4 = f(ts[i] + h, ys[i] + h * k3)
        ys[i + 1] = ys[i] + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return ts, ys
```

As you can see in the code, we don't weigh the points the same as we do with Heun's method. \
We calculate ys[i + 1] as:

$$ ys_{i+1} = k1 \cdot 1/6 + k2 \cdot 2/6 + k3 \cdot 2/6 + k4 \cdot 1/6 $$

## 7. Optimization

The `optimization()` method estimates model parameters using a Random 
Search algortihm. It adjusts parameters to minimize the Mean Squared 
Error (MSE) between the observed tumor volumes and the model predictions. 

Optimization is performed independently for each treatment group. 
The method supports different growth models and ODE-solving methods. 

The results, with an example run in the `usecase.ipynb`, include:

* Optimized model parameters
* MSE values 
* Visualizations

## 8. Model-comparing

To compare how effective the optimization is, and see which mathematical method is the best case for our data, we calculate the quality of the methods using `AICc`. For educative purposes we also added `BIC` and `AIC`, but we use `AICc` to sort the results from best to worst.

`BIC` formula is:

$$ BIC = n \cdot ln(MSE) + ln(n) \cdot k $$

The `AIC` formula is:

$$ AIC = n \cdot ln(MSE) + 2k $$

The corrected `AIC` formula is:

$$ AIC_c = n \cdot ln(MSE) + \frac{2k \cdot n}{n - k - 1} $$

We've implemented the formulas as following:

```python
BIC = n * math.log(mse) + k * math.log(n)
AIC = n * math.log(mse) + 2 * k
AICc = AIC + (2 * k * (k + 1)) / (n - k - 1) if n - k - 1 > 0 else np.inf
```

We use the if-statement for an overflow-guard.

Expanding the code's version of the corrected AIC we get the following:

$$ 2k + \frac{2k(k + 1)}{n - k - 1} = \frac{2k(n - k - 1) + 2k(k + 1)}{n - k - 1} = \frac{2kn}{n - k - 1}$$

So it's the same as the above shown version.

## 9. References
* [1] "Treatment Optimization for Tumor Growth by Ordinary Differential Equations", Kenneth Chan, https://www.researchgate.net/publication/378286328_Treatment_Optimization_for_Tumor_Growth_by_Ordinary_Differential_Equations, DOI: https://doi.org/10.47611/jsrhs.v12i4.5202
* [2] Workshops provided by the Hanze, professor Dave R.M. Langers & professor Tsjerk Wagenaar, https://brightspace.hanze.nl/d2l/le/lessons/71129/lessons/1040456 (only accessible for Hanze Bioinformatic-students, for more information approach prof. Wagenaar)
* [3] "Compartmental models (epidemiology)", wikipedia, https://en.wikipedia.org/wiki/Compartmental_models_(epidemiology)
* [4] "Euler method", wikipedia, https://en.wikipedia.org/wiki/Euler_method
* [5] "Heun's method", wikipedia, https://en.wikipedia.org/wiki/Heun%27s_method
* [6] "Runge-Kutta methods", wikipedia, https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
* [7] "Logistic function", wikipedia, https://en.wikipedia.org/wiki/Logistic_function
* [8] "Mathematical modeling of Gompertz cancer cell dynamics using real data-driven physics informed neural networks", Biswadip Pal et al., https://www.tandfonline.com/doi/full/10.1080/29937574.2026.2702897#d1e510
* [9] "Allee effect", wikipedia, https://en.wikipedia.org/wiki/Allee_effect
* [10] "Best fitting tumor growth models of the von Bertalanffy-PütterType", Manfred Kühleitner et al., https://pmc.ncbi.nlm.nih.gov/articles/PMC6624893/
* [11] "Gompertz function", wikipedia, https://en.wikipedia.org/wiki/Gompertz_function