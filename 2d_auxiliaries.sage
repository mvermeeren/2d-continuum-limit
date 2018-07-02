#!/usr/bin/env sage

### Initialize time variables, parameters, and field
def range(t,double=False):
	if includeeven:
		if double:
			return [1..2*t]
		else:
			return [1..t]
	else:
		if double:
			return [1,3..t] + [t+1,t+3..2*t]
		else:
			return [1,3..t]

var('a', latex_name='\\alpha')
var('b', latex_name='\\beta')

for i in [1..dtnumvars]:
	var('t' + str(i))

if doubletime:
	function('u',nargs=2*numvars+1)
else:
	function('u',nargs=numvars+1)

### the discrete field u_lattice
def ul(n,m,component=1):
	times = [eval('t' + str(i)) for i in [1..dtnumvars]]
	for i in range(numvars,doubletime):
		if doubletime:
			if i <= numvars:
				times[i-1] += miwaconst * (-1)^(weights[i-1]+1)* (n/weights[i-1] * a^weights[i-1])
			else: #only applies if doubletime
				times[i-1] += miwaconst * (-1)^(weights[i-1]+1)* (m/weights[i-1] * b^weights[i-1])
		else:
			times[i-1] += miwaconst * (-1)^(i+1)* (n/i * a^i + m/i * b^i)
	out = 'u('
	for i in [1..dtnumvars]:
		out += 'times[' + str(i-1) + '],'
	out = eval(out + str(component) + ')')
	return out
    
### Create admissible multi-indices [n_1,...n_numvars]
def weight(multiindex):
	return sum([multiindex[i-1]*weights[i-1] for i in [1..len(multiindex)]]) 

def indices_sub(wgt,length):
	if length == 1:
		yield [int(wgt)]
	else:
		if length in range(numvars,doubletime):
			for i in [0..int(wgt/weights[length-1])]:
				for index in indices_sub(wgt-i*weights[length-1],length-1):
					yield index+[int(i)]
		else:
			for index in indices_sub(wgt,length-1):
				yield index+[0]
                
def indices(wgt,length):
    for w in [0..wgt]:
        for index in indices_sub(w,length):
            yield index
            
### Initialize v-variables
if doubletime:
	length = 2*numvars
else:
	length = numvars
for index in indices(order,length):
	for j in [1..components]:
		if components == 1:
			string = 'v_'
		else:
			string = 'v' + str(j) + '_'
		for i in [1..length]:
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
	for index in indices(order,dtnumvars):
		for i in [1..components]:
			if field(index,i) == func:
				return index
				
### returns derivative array [t_i,t_j,...] corresponding to multiindex    
def deri(index,n=numvars):
    deri = []
    for i in [1..dt*n]:
        deri += index[i-1]*[eval('t' + str(i))]
    return deri

### converts a function of partial derivatives of u into a function of corresponding v-variables
def dertovar(func):
	for index in indices(order,dtnumvars):
		for i in [1..components]:
			func = func.subs(ul(0,0,i).diff(deri(index))==field(index,i))
	return func
	
### converts a function of v-variables into a function of corresponding partial derivatives of u
def vartoder(func):
	if doubletime:
		length = 2*numvars
	else:
		length = numvars
	for index in indices(order,length):
		for i in [1..components]:
			func = func.subs(field(index,i)==ul(0,0,i).diff(deri(index)))
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
	
### function to simplify trig and ellipitc expressions (for Q3-1 resp Q4)
def cleantrig(func):
	old = 0
	new = expand(func)
	if components == 1:
		while True:
			if hash(new) == hash(old):
				return old
			else:
				old = new
				for k in reversed([2..numvars]):
					new = new.subs( cos(v_)^k == (1 - sin(v_)^2)*cos(v_)^(k-2) )
					new = new.subs(wp2p(v_)^k == (4*wp(2*v_)+8*wp(v_))*wpp(v_)^2*wp2p(v_)^(k-2) ) 
				new = new.subs(wp2p(2*v_) == diff(wp2p(v_)/wpp(v_)*(wp(v_)-wp(2*v_)) - wpp(v_),v_)/2 )
				new = new.subs( wpp(2*v_) == wp2p(v_)/wpp(v_)*(wp(v_)-wp(2*v_)) - wpp(v_) )
				new = expand(new)
	else:
		return new
