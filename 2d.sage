#!/usr/bin/env sage

### DATA INPUT ###

### LATTICE EQUATION

switch = 'H1'
#OPTIONS - with Lagrangian: 
#lin,
#Q1, Q2, Q3, Q4w
#H1, H3sin (=H3tan),
#BSQ, GD4skew
#OPTIONS - equation only:
#H3
#Q4rat

var('delta', latex_name='\\delta')
delta = 0 # Specify parameter (commonly 0 or 1)

### DIMENSION OF MULTI-TIME
numvars = 3
includeeven = True

### LAGRANGIANS
onlyequation = False # Forget about the Lagrangians?
autosimplify = True # Simplify Lagrangians?

### CHECKS
double_eqn_exp = True # Calculate and check full array of equations?
# Not necessary if multidimensional consistency is known
# Time-intensive as code has not been optimized for speed

negdepth = 3 #How deep to check for negative order terms

elcheckdepth = numvars #To which extent to verify the EL equations
#-1 for none
#0 for quick
#[1..numvars+1] to specify order of multi-indeces to check
#Note: this does not check the 'corner equations' involving 3 coefficients.

checkcomm = True #Explicitly verify commutativity?

### REQUESTED OUTPUT
viewpdf = True # Let Tex built output pdf
vieweqnarray = False # Output the double series expansion of the quad equation
viewldisc = False # Output the series expansion of the Lagrangian
viewlpluri = False # Output the Lagrangian 2-form before simplification
viewELeqs = False # Print the first few EL equations even if they are verified

### END OF DATA INPUT ###

#########################

### import custom class for weierstrass elliptic functions
from weierstrass import *

### initialize variables
w = walltime()
var('z')
#var('g2')
warning = False

miwaconst = -2
order = 2*numvars
secondorder = False
components = 1
lagnumvars = numvars
if switch == 'BSQ':
	order = 2*numvars+2
	secondorder = True
if switch == 'GD4' or switch == 'GD4skew' or switch == 'BSQ3':
	miwaconst = 1
	components = 3
	order = 2*numvars+1 #check
	secondorder = True
	#lagnumvars = int( (numvars+1)/2 )
	#lagnumvars = int( numvars/2+1 )
	lagnumvars = numvars-1
	
### Initialize output (to console, to pdf, and to text file)
load('output.sage')
textadd(switch + ('' if delta==0 else ' with parameter ' + str(delta) ) )

### Load custom methods
load('2d_auxiliaries.sage')
load('2d_variational_calculus.sage')

### Compute limit of equations
load('2d_equation_limit.sage')
print walltime(w)

if not(onlyequation):
	### Compute limit of Lagrangians
	load('2d_lagrangian_limit.sage')
	print walltime(w)

	### Eliminate alien derivatives
	load('2d_clean_lagrangian.sage')
	print walltime(w)

### Finalize output
if warning:
    textadd('!!! WARNING(S) GENERATED !!!')
    
plaindoc.close()
latexdoc.close()

if viewpdf:
	view(output)

### If Lagrangians calculated and verified, write to file
if (not(onlyequation) and (not(warning) and (elcheckdepth >= numvars))):
	lagfile = open('lagrangians/' + filename + '-plain','w')
	lagfile.write(str(pde))
	lagfile.write('\n')
	lagfile.write(str([[utriang(cleanlag)[i,j] for j in [0..lagnumvars-1]] for i in [0..lagnumvars-1]]))
	lagfile.close()
	lagfile = open('lagrangians/' + filename + '-tex','w')
	lagfile.write(latex(pde))
	lagfile.write('\n')
	lagfile.write(latex(utriang(cleanlag)))
	lagfile.close()

