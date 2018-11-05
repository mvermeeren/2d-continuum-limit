#!/usr/bin/env sage

### Import custom class for weierstrass elliptic functions
from weierstrass import *
import os 

### Load input
load('2d_input.sage')

### Initialize variables
var('z')
warning = False

### Doubly infinite hierarchy?
doubletime = bool(switch == 'H3double')
if doubletime:
	dt = 2
	weights = 2*[1..numvars]
	autosimplify = False #Simplification algorithm not implemented for double hierarchy
else:
	dt = 1
	weights = [1..numvars]
dtnumvars = dt*numvars

# Set default values
components = 1
order = 2*numvars
secondorder = 0
lagnumvars = numvars
miwaconst = -2
squarearray = True

# Override default values if set in input file
if switch in multicompontent_list:
	components = multicompontent_list[switch]
if switch in order_list:
	order = order_list[switch]
if switch in secondorder_list:
	secondorder = secondorder_list[switch]
if switch in lagnumvars_list:
	lagnumvars = lagnumvars_list[switch]
if switch in miwa_list:
	miwaconst = miwa_list[switch]
if switch in squarearray_list:
	squarearray = squarearray_list[switch]
	
### Initizalize variables and load methods
load('2d_auxiliaries.sage')
load('2d_variational_calculus.sage')

### Impose equations listed in input file, if any
constraints = constraint_list(switch)
imposedpdes = pde_list(switch)
	
### Open output files
if not os.path.exists('datadump/'):
    os.makedirs('datadump/')
filename = switch
if not(delta == 0):
	filename += '-delta'
filename += '-' + str(numvars)
if includeeven:
	filename += 'full'
plaindoc = open('datadump/' + filename + '-plain','w')
latexdoc = open('datadump/' + filename + '-tex','w')


### Customized latex function to use to export tex file
### For immediate display, the built in latex() should be used
def cleanlatex(eqn):
	string = latex(eqn)
	for power in [1..9]:
		string = string.replace('^{'+str(power)+'}','^'+str(power))
	string = string.replace('{v_{','v_{') 
	string = string.replace('}}','}') 
	string = string.replace('\, ','') 
	string = string.replace('\left(v_{}\right)','(v)') 
	return string

### print to latex and plaintext files
### if all==True, print to console and pdf as well
output = []	

def latexadd(eqn,all=True):
	global output
	print ""
	print eqn
	print ""
	if all:
		output += ['\\scalebox{.1}{$' + latex(eqn) + '$}'] #['\\tiny ' + latex(eqn)]
	if not(latexdoc.closed):
		latexdoc.write(cleanlatex(eqn) + '\n\n')
	if not(plaindoc.closed):
		plaindoc.write(str(eqn) + '\n\n')
		
def textadd(string,all=True):
	global output
	print "%.2f" % (walltime(w)) + "s: " + string
	if all:
		output += [latex('') + '\\scalebox{.1}{' + string + '}']
	if not(latexdoc.closed):
		latexdoc.write(string + '\n\n')
	if not(plaindoc.closed):
		plaindoc.write(string + '\n\n')
		
### start timing
w = walltime()

### Print mission statement
textadd('EQUATION ' + switch + ('' if delta==0 else ' with parameter ' + str(delta) ) )

### Compute limit of equations
load('2d_equation_limit.sage')

if not(onlyequation):
	### Compute limit of Lagrangians
	load('2d_lagrangian_limit.sage')
	### Eliminate alien derivatives
	load('2d_clean_lagrangian.sage')

### Finalize output
if warning:
    textadd('!!! WARNING(S) GENERATED !!!')
    
plaindoc.close()
latexdoc.close()

if viewpdf:
	view(output)

### If Lagrangians calculated and verified, write to file
if (not(onlyequation) and not(str(L)=='0') and (not(warning) and (elcheckdepth >= numvars))):
	if not os.path.exists('lagrangians/'):
		os.makedirs('lagrangians/')
	lagfile = open('lagrangians/' + filename + '-plain','w')
	lagfile.write(str(pde))
	lagfile.write('\n\n')
	lagfile.write(str(c))
	lagfile.write('\n\n')
	if doubletime:
		lagfile.write(str(cc))
		lagfile.write('\n\n')
		lagfile.write(str([[cleanlag[i,j] for j in [0..lagnumvars-1]] for i in [0..lagnumvars-1]]))
	else:
		lagfile.write(str([[utriang(cleanlag)[i,j] for j in [0..lagnumvars-1]] for i in [0..lagnumvars-1]]))
	lagfile.close()
	
	lagfile = open('lagrangians/' + filename + '-tex','w')
	lagfile.write(cleanlatex(pde))
	lagfile.write('\n\n')
	lagfile.write(cleanlatex(c))
	lagfile.write('\n\n')
	if doubletime:
		lagfile.write(cleanlatex(cc))
		lagfile.write('\n\n')
		lagfile.write(cleanlatex(cleanlag))
	else:
		lagfile.write(cleanlatex(utriang(cleanlag)))
	lagfile.close()

