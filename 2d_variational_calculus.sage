#!/usr/bin/env sage

### integrate away terms with time derivatives as far as possible, and
### integrate away excessively high space derivatives (e.g. v_111*v1 to -d(v_11*v))
def vintegrate(func):
	out = 0
	for j in (range(numvars)[1:] + ['']):
		if j == '':
			addorder = numvars
		else:
			addorder = numvars-j
		for i in reversed([1..addorder]): 
			for component in [1..components]:
				for t in summands(func):
					if components == 1:
						vstring = 'v_'
					else:
						vstring = 'v' + str(component) + '_'
					new = 0
					dhighest = diff(t, eval(vstring+i*'1'+str(j)))
					intcandidate = dhighest * eval(vstring+(i-1)*'1'+str(j))
					#check that it is highest derivative
					if (i+1)*'1'+str(j) in str(t) or bool(intcandidate == 0):
						new = 0
					#catch logaritmic derivatives
					elif not( vstring in str(intcandidate) ):
						new = intcandidate*log(eval(vstring+(i-1)*'1'+str(j)))
					elif components == 1 and not( vstring in str( expand(intcandidate*sin(v_)/cos(v_)/v_) ) ):
						new = intcandidate*sin(v_)/cos(v_)/v_ * log(sin(eval(vstring+(i-1)*'1'+str(j))))
					elif components == 1 and 'v_*wzeta(2*v_)' in str(intcandidate):
						new = intcandidate/(v_*wzeta(2*v_)) * (1/2)*log(wsigma(2*v_)) #NOT ROBUST
					#leave other elliptic terms be
					#elif 'wzeta(2*v_)' in str(intcandidate):
					#	new = 1
					#highest derivative has to occur linearly
					elif not(diff(dhighest, eval(vstring+i*'1'+str(j))) == 0):
						new = 0
					#else:
					#	intcandidate = diff(t, eval('v_'+i*'1'+str(j))) * eval('v_'+(i-1)*'1'+str(j))
					#handle rational derivatives
					else:
						newterms = summands(intcandidate)
						for term in newterms:
							if not(term == 1):
								coeff = diff(term, eval(vstring+(i-1)*'1'+str(j))) / term*eval(vstring+(i-1)*'1'+str(j))
								if not(coeff == 0):
									new += term/coeff
					#update function, add term to integral
					func = expand(func - vdiff(new,t1))
					out += new
	return out

#deletes powers of parameters higher than truncorder
def killpowers(f,truncorder=numvars):
	return f.subs([a^(truncorder+i)==0 for i in [1..truncorder]]).subs([b^(truncorder+i)==0 for i in [1..truncorder]])

### variational derivative in terms of v-variables
var('target')
def varder(dir,func,var,component=1):
    result = 0
    freeorder = order - weight(var)
    for index in indices(freeorder,numvars):
        if sum(index) == sum(index[d-1] for d in dir):
            realindex = [index[i-1]+var[i-1] for i in [1..numvars]]
            a = diff(func.subs(field(realindex,component)==target),target)
            a = vdiff(a.subs(target==field(realindex,component)),deri(index,numvars))
            result += (-1)^sum(index) * a
    return result
    

#@parallel
#def pcheck(i,j,func):
#	out = expand(replace(func))
#	if out == 0:
#		return 0
#	else:
#		return cleantrig(out).simplify_rational()
#	    #return trig_to_exp(out).simplify_full()

def check(matrix,message="All entries = 0 mod eqns",verbose=True):
	size = matrix.dimensions()[1]
	#checkgen = pcheck(flatten([[(i,j,matrix[i-1,j-1]) for j in [1..size]] for i in [1..size]],max_level=1))
    
	remainder = 0*copy(matrix)
	for i in [1..size]:
		for j in [1..size]:
			out = expand(replace(matrix[i-1,j-1]))
			if out == 0:
				remainder[i-1,j-1] = 0
			else:
				remainder[i-1,j-1] = cleantrig(out).simplify_rational()
	if verbose:
		if remainder == 0*remainder:
			textadd(message)
		else:
			textadd("Nonvanishing terms!")
			latexadd(remainder)
	return remainder == 0*remainder


### variational derivatives
@parallel
def pvarder(i,j,func,index,component=1):
	if func == 0:
		return 0
	else:
	    return varder([i,j],func,index,component)

def varders(matrix,index,edge=False,component=1):
	eindex = [copy(index) for k in [1..lagnumvars]]
	for i in [1..lagnumvars]:
		eindex[i-1][i-1] = eindex[i-1][i-1]+1
	if edge:
	    vgen = pvarder(flatten([[(i,j,matrix[i-1,j-1],eindex[i-1],component) for i in range(lagnumvars)] for j in [1..lagnumvars]],max_level=1))
	else:
	    vgen = pvarder(flatten([[(i,j,matrix[i-1,j-1],index,component) for i in [1..lagnumvars]] for j in [1..lagnumvars]],max_level=1))

	out = 0*copy(matrix)
	for i in vgen:
		out[i[0][0][0]-1,i[0][0][1]-1] = i[1]
	return out
### take difference of rows
def diff_in_columns(matrix):
	out = copy(matrix)
	for i in range(lagnumvars):
		for j in range(lagnumvars):
			if not(i==j):
				if j==1:
					out[i-1,j-1] = matrix[range(lagnumvars)[1]-1,j-1] - matrix[i-1,j-1]
				else:
					out[i-1,j-1] = matrix[0,j-1] - matrix[i-1,j-1]
	return out

### calculate multi-time EL expressions of 1st kind
@parallel
def elplane(matrix,index,component=1):
	global warning
	relevantmatrix = copy(matrix)
	for i in [1..lagnumvars]: # eliminate irrelevant terms
		if index[i-1] > 0:
			for j in [1..lagnumvars]:
				relevantmatrix[i-1,j-1] = 0
				relevantmatrix[j-1,i-1] = 0
	vararray = varders(relevantmatrix,index,False,component)
	if not(check(vararray,verbose=False)):
		warning = True
		textadd(str(index) + " - " + str(component) + " - plane")
		check(vararray,verbose=True)
	#return warning

### calculate multi-time EL expressions of 2nd kind
@parallel
def eledge(matrix,index,component=1):
	global warning
	relevantmatrix = copy(matrix)
	for i in [1..lagnumvars]: # eliminate irrelevant terms
		if index[i-1] > 0:
			for j in [1..lagnumvars]:
				if not(i==j):
					relevantmatrix[j-1,i-1] = 0
	vararray = varders(relevantmatrix,index,True,component)
	diffvararray = diff_in_columns(vararray)
	if not(check(diffvararray,verbose=False)):
		warning = True
		textadd(str(index) + " - " + str(component) + " - edge")
		latexadd(vararray)
		latexadd(diffvararray)


### check multi-time EL equations for given matrix of coefficients
def elcheck(matrix):
	triang = utriang(matrix)
	global warning
	if elcheckdepth > -1:
		for component in [1..components]:
			vararray = varders(triang,[0 for i in [1..numvars]],False,component)
			if viewELeqs:
				textadd("Multi-time EL equations; - " + str(component) + " - plane:")
				latexadd(vararray)
			if not(check(vararray)):
				warning = True
		
		for component in [1..components]:
			vararray = varders(matrix,[0 for i in [1..numvars]],True,component)
			if viewELeqs:
				textadd("Multi-time EL equations; - " + str(component) + " - edge:")
				latexadd(vararray)
			diffvararray = diff_in_columns(vararray)
			if not(check(diffvararray,"Columns agree mod eqns")):
				warning = True
	
	if elcheckdepth > 0:
		textadd("Deep EL check...")
		#if elcheckdepth > 2:
		#	print("elcheckdepth = " + str(elcheckdepth) + ". This could take a while...")
		for index in indices(elcheckdepth,numvars):
			for component in [1..components]:
				elplane(triang,index,component)
		#print("Half way there...")
		for index in indices(elcheckdepth-1,numvars):
			for component in [1..components]:
				eledge(matrix,index,component)
		if not(warning):
			textadd("Deep EL check successful.")
