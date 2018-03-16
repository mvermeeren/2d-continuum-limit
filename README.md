# 2d-continuum-limit
A SageMath implementation of continuum limits of 2-dimensional lattice equations.

This is an implementation of the continuum limit for two-dimensional multidimensionally consistent lattice equations and their pluri-Lagrangian structure.

The continuum limit procedure is discussed in

	Mats Vermeeren, “Continuum limits of pluri-Lagrangian systems” 
	arxiv.org/abs/1706.06830

This package consists of the following files:

	* 2d.sage (Main)
	* 2d_auxiliaries.sage
	* 2d_clean_lagrangian.sage
	* 2d_equation_limit.sage
	* 2d_lagrangian_limit.sage
	* 2d_variational_calculus.sage
	* output.sage
	* weierstrass.py

This software is written by

	Mats Vermeeren
	TU Berlin
	vermeeren@math.tu-berlin.de
	
and published under the MIT License (see 'LICENCE' for more details)


RUNNING THE PROGRAM

Execute the file '2d.sage' in SageMath, preferably version 7.5.1 or higher

The input has to be given within the '### DATA INPUT ###' section of '2d.sage'.
Most importantly the lattice equation needs to be specified by setting 'switch'
and the dimension of the multi-time in which it is embedded needs to be specified 
with 'numvars'

Additional options are described within the file '2d.sage'


OUTPUT

Aside from printin in the console, the program will create log files in the 
subfolder 'datadump'. If the parameter 'viewpdf' is set to True, a pdf of the
output will be compiled.

If a continuous pluri-Lagrangian structure is successfully computed, it will be 
written to a file in the subfolder 'lagrangians'


ADDING MORE EQUATIONS OR LAGRANGIANS

Custom equations and Lagrangians can be added in the files 
'2d_equation_limit.sage' and  '2d_lagrangian_limit.sage' respectively. 
A unique 'switch' should be introduced to select them.
