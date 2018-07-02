#!/usr/bin/env sage

textadd("LAGRANGIAN")

scaling = 1/miwaconst * diagonal_matrix([-(-1)^i*i for i in [1..lagnumvars]]) 

### First series expansion
if doubletime:
	series = L.series(b,lagnumvars+1).truncate()
else:
	series = L.series(b,lagnumvars).truncate()

### Check for nonpositive order terms
locwarning = False
for i in [1-negdepth..0]:
    test = expand(dertovar(series.coefficient(b,i).series(a,lagnumvars).truncate()).simplify_full())
    if not(test == 0):
        locwarning = True
        textadd("Term at order " + str(i) + "? - Will check again later")
        latexadd(test)
if ((locwarning == False) & (negdepth > 0)):
    textadd("No nonpositive order terms detected")

### Second series expansion
@parallel
def expanda(i):
	series2 = expand( series.coefficient(b,i).series(a,lagnumvars+1).truncate() )
	if doubletime:
		return [expand(dertovar( series2.coefficient(a,j) )) for j in [1..lagnumvars]]
	else:
		return [0 for j in [1..i]] + [expand(dertovar( series2.coefficient(a,j) )) for j in [i+1..lagnumvars]]

cgen = expanda([i for i in [1..lagnumvars + dt -2]])

Lcoeff = [[0 for j in [1..lagnumvars]] for i in [1..lagnumvars]]
for i in cgen:
    Lcoeff[i[0][0][0]-1] = i[1]

textadd("L disc = ",viewldisc)
latexadd(matrix(Lcoeff),viewldisc)

if doubletime:
	Lcoeff = matrix(Lcoeff)
else:
	Lcoeff = matrix(Lcoeff) - matrix(Lcoeff).transpose()

### Euler-Maclaurin correction ###

#vertical EM operator
@parallel
def diffEMcell(row,col,array):
	out = 0*a
	for i in [1..row]:
		if includeeven or is_odd(i):
			out += (-1)^(i+1) * miwaconst/i * vdiff(array[row-i][col], eval('t'+str(i + (dt-1)*numvars )))
	return out

def EMdiff(array):
	out = [[0 for i in [1..lagnumvars]] for j in [1..lagnumvars]]
	cellgen = diffEMcell(flatten([[(i,j,array) for i in [0..lagnumvars-1]] for j in [0..lagnumvars-1]],max_level=1))
	for i in cellgen:
		out[i[0][0][0]][i[0][0][1]] = i[1]
	return matrix(out)

#powers of vertical EM operator
EMdiffs = [matrix(lagnumvars,lagnumvars) for k in [0..lagnumvars]]
EMdiffs[0] = Lcoeff
for i in [1..lagnumvars]:
    EMdiffs[i] = EMdiff(EMdiffs[i-1])

#vertical correction
lagarray = matrix(EMdiffs[0]) - 1/2*EMdiffs[1] + sum(bernoulli(2*i)/factorial(2*i)*EMdiffs[2*i] for i in [1..lagnumvars/2])

#powers of horizontal EM operator
@parallel
def diffEMcell(row,col,array):
	out = 0*a
	for i in [1..row]:
		if includeeven or is_odd(i):
			out += (-1)^(i+1) * miwaconst/i * vdiff(array[row-i][col], eval('t'+str(i)))
	return out
		
EMdiffs = [matrix(lagnumvars,lagnumvars) for k in [0..lagnumvars]]
EMdiffs[0] = lagarray
for i in [1..lagnumvars]:
    EMdiffs[i] = EMdiff(EMdiffs[i-1].transpose()).transpose()

#horizontal correction
lagarray =  matrix(EMdiffs[0]) - 1/2*EMdiffs[1] + sum(bernoulli(2*i)/factorial(2*i)*EMdiffs[2*i] for i in [1..lagnumvars/2])
lagarray = scaling * lagarray * scaling

triang = utriang(lagarray)

### Ouptut
textadd("L pluri = ",viewlpluri)
if doubletime:
	latexadd(lagarray,viewlpluri)
else:
	latexadd(triang,viewlpluri)
