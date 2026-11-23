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

TODO: toevoegen dmv command `tree /f /a` aan het einde

```

```

## 3. Background information

This casus is inspired by the student Kenneth Chan of the Pine Crest School, Claremont Mckenna College. They wrote a short paper called "[Treatment Optimization for Tumor Growth by Ordinary Differential Equations](https://doi.org/10.47611/jsrhs.v12i4.5202)".
In this paper applies the author among other things several different models of tumor growth to various existing datasets in order to identify the best-fitting model.

Mathematical models, i.e. ordinary differential equations (ODEs), can be used to model tumor growth and study treatment of cancer [[1]](#9-references). 
This makes it essential for bioinformatics to understand and be able to apply such models and equations on patient data. 

## 4. Installation & Requirements

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

Describes growth that initially accelerates but gradually 
slows down as the tumor approaches its maximum volume. 

$$
\frac{dV}{dt}=c\cdot V\cdot(V_{\max}-V)
$$

## 6. Methods

### a. Euler method

### b. Heun's method

### c. Runge-Kutta methods

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

## 9. References
* [1] "Treatment Optimization for Tumor Growth by Ordinary Differential Equations", Kenneth Chan, https://www.researchgate.net/publication/378286328_Treatment_Optimization_for_Tumor_Growth_by_Ordinary_Differential_Equations, DOI: https://doi.org/10.47611/jsrhs.v12i4.5202
* [2] Workshops provided by the Hanze, professor Dave R.M. Langers & professor Tsjerk Wagenaar, https://brightspace.hanze.nl/d2l/le/lessons/71129/lessons/1040456 (only accessible for Hanze Bioinformatic-students, for more information approach prof. Wagenaar)
* [3] "Compartmental models (epidemiology)", wikipedia, https://en.wikipedia.org/wiki/Compartmental_models_(epidemiology)
* [4] "Euler method", wikipedia, https://en.wikipedia.org/wiki/Euler_method
* [5] "Heun's method", wikipedia, https://en.wikipedia.org/wiki/Heun%27s_method
* [6] "Runge-Kutta methods", wikipedia, https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods