#!/usr/bin/env sage

### Initialize time variables, parameters, and field
def range(t):
    if includeeven:
        return [1..t]
    else:
        return [1,3..t]

var('a', latex_name='\\alpha')
var('b', latex_name='\\beta')

for i in [1..numvars]:
    var('t' + str(i))
    
function('u',nargs=numvars+1)

### the discrete field u_lattice
def ul(n,m,component=1):
    times = [eval('t' + str(i)) for i in [1..numvars]]
    for i in [1..numvars]:
        if includeeven | is_odd(i):
            times[i-1] += miwaconst * (-1)^(i+1)* (n/i * a^i + m/i * b^i)
    out = 'u('
    for i in [1..numvars]:
        out += 'times[' + str(i-1) + '],'
    out = eval(out + str(component) + ')')
    return out
    
### Create admissible multi-indices [n_1,...n_numvars]
def weight(multiindex):
    return sum([i*multiindex[i-1] for i in [1..len(multiindex)]]) 

def indices_sub(wgt,len):
    if len == 1:
        yield [wgt]
    else:
        if includeeven | is_odd(len):
            for i in [0..int(wgt/len)]:
                for index in indices_sub(wgt-i*len,len-1):
                    yield index+[i]
        else:
            for index in indices_sub(wgt,len-1):
                yield index+[0]
                
def indices(wgt,len):
    for w in [0..wgt]:
        for index in indices_sub(w,len):
            yield index
            
### Initialize v-variables
for index in indices(order,numvars):
	for j in [1..components]:
		if components == 1:
			string = 'v_'
		else:
			string = 'v' + str(j) + '_'
		for i in [1..numvars]:
			string += str(i)*index[i-1]
		var(string)

### returns v-variable corresponding to given multiindex (and component)
def field(index,i=1):
	if components ==1:
		string = 'v_'
	else:
		string = 'v' + str(i) + '_'
	for i in [1..len(index)]:
		string += str(i)*index[i-1]
	return eval(string)

### returns multiindex corresponding to given v-variable
def fieldtoindex(func):
	for index in indices(order,numvars):
		for i in [1..components]:
			if field(index,i) == func:
				return index
				
### returns derivative array [t_i,t_j,...] corresponding to multiindex    
def deri(index,n=numvars):
    deri = []
    for i in [1..n]:
        deri += index[i-1]*[eval('t' + str(i))]
    return deri

### converts a function of partial derivatives of u into a function of corresponding v-variables
def dertovar(func):
	for index in indices(order,numvars):
		for i in [1..components]:
			func = func.subs(ul(0,0,i).diff(deri(index,numvars))==field(index,i))
	return func
	
### converts a function of v-variables into a function of corresponding partial derivatives of u
def vartoder(func):
	for index in indices(order,numvars):
		for i in [1..components]:
			func = func.subs(field(index,i)==ul(0,0,i).diff(deri(index,numvars)))
	return func
    
### takes time-derivatives on the level of v-variables
def vdiff(func,*var):
    func = vartoder(func)
    func = diff(func,*var)
    func = dertovar(func)
    return func

### return terms of a sum
def summands(func):
	if func == 0:
		return []
	elif type(func) == type(1+a) and func.operator() == (1+a).operator():
		return func.operands()
	else:
		return [func]
		
### returns upper triangular part of a matrix
def utriang(matrix):
	triang = 0*matrix
	for i in [0..matrix.dimensions()[1]-1]:
		triang[i,i:] = matrix[i,i:]
	return triang
