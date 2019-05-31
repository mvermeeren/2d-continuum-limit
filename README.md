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
	* 2d_input.sage (Specify input here)
	* 2d_lagrangian_limit.sage
	* 2d_variational_calculus.sage
	* weierstrass.py

This software is written by

	Mats Vermeeren
	TU Berlin
	vermeeren@math.tu-berlin.de
	
and published under the MIT License (see 'LICENCE' for more details)


RUNNING THE PROGRAM

	This code was developed in Sage 7.5.1. 
	
	In Sage 8, there is a so far unindentified bug in the series expansion of very large expressions.
	
Execute the file '2d.sage' in SageMath.

All input is given in the file '2d_input.sage'.

Most importantly the lattice equation needs to be specified by setting 'switch'
and the dimension of the multi-time in which it is embedded needs to be specified 
with 'numvars'.

Additional options are described in the comments of the file '2d_input.sage'.


OUTPUT

Aside from printin in the console, the program will create log files in the 
subfolder 'datadump'. If the parameter 'viewpdf' is set to True, a pdf of the
output will be compiled.

If a continuous pluri-Lagrangian structure is successfully computed, it will be 
written to a file in the subfolder 'lagrangians'


ADDING MORE EQUATIONS OR LAGRANGIANS

Custom equations and Lagrangians can be added in the method get_equation(switch)
of the file '2d_input.sage'. A unique 'switch' should be introduced to select them.
