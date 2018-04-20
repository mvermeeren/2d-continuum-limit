#!/usr/bin/env sage

### Check if a quantity contains only native derivatives
@parallel
def is_native(times,func):
	times = [1] + times
	out = True
	for index in indices(order,numvars):
		alienindex = copy(index)
		for i in times:
			alienindex[i-1] = 0
		if not(alienindex == [0 for i in [1..numvars]]):
			if not(diff(func,field(index)) == 0):
				out = False
	return out

### Isolate terms with alien derivatives
#ONLY NEEDED FOR TROUBLESHOOTING?
def alien_terms(times,func):
	terms = summands(func)
	gen = is_native([(times,term) for term in terms])
	out = 0
	for i in gen:
		if not(i[1]):
			out += i[0][0][1]
	return out
	
### Eliminate products of time derivatives
# Adds a term that as a double zero on solutions - does not affect EL eqns
relpde = [] #relevant equations
fv = [] #create flow variables
for e in allpde:
	#if weight(fieldtoindex(e.lhs())) <= 2*numvars-3: # might be too strict
	if weight(fieldtoindex(e.lhs())) <= 2*numvars-2: # might be too strict
		relpde += [e]
		var('r' + str(e.lhs()))
		fv += [ [eval('r' + str(e.lhs())) , weight(fieldtoindex(e.lhs())) ] ]

@parallel
def double0(func,i,j):
	if func in RR:
		return func
		
	#replace to flow vars    
	temp = expand(func)
	for k in [0..len(relpde)-1]:
		e = relpde[k]
		temp = expand(temp.subs({e.lhs(): fv[k][0] + e.rhs()}))
	
	terms = summands(temp)
	out = 0
	
	#substitute double zeros
	if secondorder > 0:
		cutoff = 2*(numvars)
	else:
		cutoff = 2*(numvars)-1
	for t in terms:
		include = True
		for i in [1..len(fv)]:
			var1 = fv[i-1]
			for var2 in fv[i-1:]:
				if var1[1] + var2[1] <= cutoff:
					if (t.has(var1[0])) and (expand(t/var1[0]).has(var2[0])):
						include = False
		if include:
			out += t

	#replace to original vars    
	for k in [0..len(relpde)-1]:
		e = relpde[k]
		out = expand(out.subs({fv[k][0]: e.lhs() - e.rhs()}))
	
	return out
	
### Check again for nonpositive order terms
for i in [1-negdepth..0]:
	test = expand(double0(dertovar(series.coefficient(b,i).series(a,lagnumvars).truncate()),0,0))
	if not(test == 0):
		warning = True
		textadd("Term at order " + str(i) + "!")
		latexadd(test.simplify_full())
if ((warning == False) & (negdepth > 0)):
	textadd("No nonpositive order terms detected")

### Simplify by double zeroes
cleantriang = copy(triang)
cleangen = double0(flatten([[(triang[i-1,j-1],i,j) for j in [1..lagnumvars]] for i in [1..lagnumvars]],max_level=1))
for i in cleangen:
	cleantriang[i[0][0][1]-1,i[0][0][2]-1] = i[1]

### Find suitable 1-form c1 dt1 + c2 dt2 + ...
c = [0 for i in [1..lagnumvars]]
if autosimplify:
	for i in [1..lagnumvars]:
		integral = vintegrate(cleantriang[0,i-1])
		if integral == 0:
			c[i-1] = 0
		else:
			c[i-1] = -integral.combine()
	
	textadd("Simplifying 1-form:")
	latexadd(c)

### calculate d( c1 dt1 + c2 dt2 + ... )
corr = 0*copy(triang)
for i in range(lagnumvars):
	for j in [1..lagnumvars]:
		if i < j:
			corr[i-1,j-1] += vdiff(c[j-1],eval('t'+str(i)))
		if i > j:
			corr[j-1,i-1] += -vdiff(c[j-1],eval('t'+str(i)))

### Add d( c1 dt1 + c2 dt2 + ... )
for i in [1..lagnumvars]:
	for j in [1..lagnumvars]:
		cleantriang[i-1,j-1] = cleantrig( corr[i-1,j-1] + cleantriang[i-1,j-1] )

### In case the previous step introduces new producs of time derivatives
cleangen = double0(flatten([[(cleantriang[i-1,j-1],i,j) for j in [1..lagnumvars]] for i in [1..lagnumvars]],max_level=1))
for i in cleangen:
	func = cleantrig(i[1])
	if func == 0:
		cleantriang[i[0][0][1]-1,i[0][0][2]-1] = 0
	else:
		cleantriang[i[0][0][1]-1,i[0][0][2]-1] = func.combine()

### Output
textadd('Simplified Lagrangian:')
latexadd(cleantriang)
cleanlag = cleantriang - transpose(cleantriang)
elcheck(cleanlag)
