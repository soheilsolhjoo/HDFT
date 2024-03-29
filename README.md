# Hot Deformation Fitting Tool (HDFT) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5226216.svg)](https://doi.org/10.5281/zenodo.5226216)

## Table of Contents
1. [Introduction](#introduction)
2. [Copyright and license](#copyright-and-license)
3. [How to cite](#how-to-cite)
4. [Start guide](#start-guide)
5. [Data file formatting](#data-file-formatting)
6. [Temperature flag](#temperature-flag)

## Introduction
HDFT code (written in MATLAB) is developed to study hot deformation's stress-strain data sets. Various models are available, including the hyperbolic sine model of Sellars and Tegart [1], its revised version [2], and its limiting exponential and power-law models. Moreover, a code for training Artificial Neural Networks is available as well.

[1] C. M. Sellars and W. McTegart, “On the mechanism of hot deformation”, Acta Metallurgica 1966, 14(9), 1136–1138.
<br />[2] S. Solhjoo, "Revisiting the common practice of Sellars and Tegart’s hyperbolic sine constitutive model: a critical note", Modelling 2022, 3(3), 359-373, doi:[10.3390/modelling3030023](https://doi.org/10.3390/modelling3030023).

## Copyright and license
The HDFT code is licensed under the [GNU GPLv3 license](https://choosealicense.com/licenses/gpl-3.0/). The details of this license can be found in the files COPYRIGHT.txt. The code is written by Dr. Soheil Solhjoo.

## How to cite
The essentials of this code have been outlined in:
- S. Solhjoo, "Revisiting the common practice of Sellars and Tegart’s hyperbolic sine constitutive model: a critical note", Modelling 2022, 3(3), 359-373, doi:[10.3390/modelling3030023](https://doi.org/10.3390/modelling3030023).
```
@article{solhjoo2022modelling,
  title={Revisiting the Common Practice of Sellars and Tegart’s Hyperbolic Sine Constitutive Model},
  author={Solhjoo, Soheil},
  journal={Modelling},
  volume={3},
  number={3},
  pages={359--373},
  year={2022},
  publisher={MDPI}
}
```
Moreover, upon using the HDFT code, please cite the code's DOI and its GitHub page as well:
- https://github.com/soheilsolhjoo/HDFT
- https://doi.org/10.5281/zenodo.5226216

## Start guide
The HDFT code is written in MATLAB 2021a and tested on a machine running on Windows 10.
To run the HDFT code, the stress-strain data should be prepared as described in [Data file formatting](#Data-file-formatting). Then, various functions can be called to perform the fitting procedure. File "HD_runner.m" is a sample file for calling functions, and each has a detailed description on how it works.

## Data file formatting
The **input data** must be stored in files with the name **#.data**, where "#" is a number; e.g., "1.data", "2.data", etc. Each file must have the following format:
<br />line 1: *Temperature*
<br />line 2: *Strain-rate*
<br />line 3: *strain   stress*
<br />line 4: *strain   stress*
<br />  ...

For example:
<br />100
<br />0.1
<br />0	100
<br />0.01	110
<br />  ...
<br />0.1	200

The main **output** is saved in **HDFT_fit.mat** containing various structures. The general one is "input" that is the processed input data. Moreover, each function generates a separate structure for itself:
<br />PL# (HD_power_law)
<br />EX# (HD_exponential)
<br />ST## (HD_sinh_conventional)
<br />STR### (HD_sinh_revisited)
<br />ANN#_# (HD_ANN)

***NOTE**: HDFT does not save two structures with the same name, and the newer structure overwrites the previous one, if existed.*

## Temperature flag
Among the controlling input arguments, ***CK*** is the first one for all the functions, which determines the temperature unit: *Celsius -> CK=1* and *Kelvin -> CK=0*.
