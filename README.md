# HotDeformation [![DOI](https://zenodo.org/badge/339745520.svg)](https://zenodo.org/badge/latestdoi/339745520)
The code is published under [GNU GPLv3 license](https://choosealicense.com/licenses/gpl-3.0/).


This code is developed to study hot deformation's stress-strain data sets. Various models are available, including the hyperbolic sine model of Sellars and Tegart [1], its revised version [2], and its limiting exponential and power-law models. Moreover, a code for training Artificial Neural Networks is available as well.

[1] C. M. Sellars and W. McTegart, “On the mechanism of hot deformation”, Acta Metallurgica 14.9 (1966), pp. 1136–1138.
<br />[2] S. Solhjoo. "Revisiting the common practice of Sellars and Tegart’shyperbolic sine constitutive model: a critical note", to be submitted.

## File formatting
The input data must be stored in files with the name **#.data**, where "#" is a number; e.g., "1.data", "2.data", etc. Each file must have the following format:
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


File "HD_runner.m" is a sample file for calling functions, and each has a detailed descrptoin on how it works. Among the controlling input arguments, ***CK*** is crucial, which should be decided by the user based on the temperature unit: *Celcius -> CK=1* and *Kelvin -> CK=0*.
