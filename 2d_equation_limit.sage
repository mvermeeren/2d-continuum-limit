#!/usr/bin/env sage

### EQUATION LIST ###

squararray = True

quad = 0*z

if switch == 'lin':
    quad = (1/b + 1/a)*(ul(0,1) - ul(1,0)) + (1/b - 1/a)*(ul(1,1) - ul(0,0))
    
elif switch == 'Q1':
    quad = 1/b^2*(ul(0,0)-ul(0,1))*(ul(1,0)-ul(1,1)) - 1/a^2*(ul(0,0)-ul(1,0))*(ul(0,1)-ul(1,1)) + delta^2*(a^2-b^2)
elif switch == 'Q2':
    quad = 1/b^2*(ul(0,0)^2-ul(0,1)^2)*(ul(1,0)^2-ul(1,1)^2) - 1/a^2*(ul(0,0)^2-ul(1,0)^2)*(ul(0,1)^2-ul(1,1)^2) + (a^2-b^2)*(ul(0,0)^2+ul(1,0)^2+ul(0,1)^2+ul(1,1)^2) - (a^2-b^2)*(a^4-a^2*b^2+b^4)
    
elif switch == 'Q3':
	if delta == 0:
		quad = -sin(a^2-b^2)*(exp(I*ul(0,0)+I*ul(1,1))+exp(I*ul(1,0)+I*ul(0,1))) + sin(a^2)*(exp(I*ul(0,0)+I*ul(1,0))+exp(I*ul(0,1)+I*ul(1,1))) - sin(b^2)*(exp(I*ul(0,0)+I*ul(0,1))+exp(I*ul(1,0)+I*ul(1,1)))
	else:
		quad = -sin(a^2-b^2)*( cos(ul(0,0))*cos(ul(1,1)) + cos(ul(1,0))*cos(ul(0,1)) ) + sin(a^2)*( cos(ul(0,0))*cos(ul(1,0)) + cos(ul(0,1))*cos(ul(1,1)) ) - sin(b^2)*( cos(ul(0,0))*cos(ul(0,1)) + cos(ul(1,0))*cos(ul(1,1)) ) + delta*sin(a^2-b^2)*sin(a^2)*sin(b^2)
	quad = quad/(a*b)^2

elif switch == 'Q4rat':
	var('g2','g3')
	var('z')
	wpcoeff = [g2/20,g3/28]
	for k in [4..numvars]:
		wpcoeff += [ 3/(2*k+1)/(k-3) * sum([ wpcoeff[m-2]*wpcoeff[k-m-2] for m in [2..k-2]]) ]
	def wpseries(func,prime=False):
		out = z^(-2) + sum([wpcoeff[k-2]*z^(2*k-2) for k in [2..numvars] ])
		if prime==False:
			return out.subs(z==func)
		else:
			return diff(out,z).subs(z==func)
	quad = wpseries(a^2,True) * ( (ul(0,0) - wpseries(b^2))*(ul(0,1) - wpseries(b^2)) - (wpseries(a^2)-wpseries(b^2))*(wpseries(b^2-a^2)-wpseries(b^2)) ) * ( (ul(1,0) - wpseries(b^2))*(ul(1,1) - wpseries(b^2)) - (wpseries(a^2)-wpseries(b^2))*(wpseries(b^2-a^2)-wpseries(b^2)) )
	quad += wpseries(b^2,True) * ( (ul(0,0) - wpseries(a^2))*(ul(1,0) - wpseries(a^2)) - (wpseries(b^2)-wpseries(a^2))*(wpseries(b^2-a^2)-wpseries(a^2)) ) * ( (ul(0,1) - wpseries(a^2))*(ul(1,1) - wpseries(a^2)) - (wpseries(b^2)-wpseries(a^2))*(wpseries(b^2-a^2)-wpseries(a^2)) )
	quad += - wpseries(a^2,True)*wpseries(b^2,True)*wpseries(b^2-a^2,True) * (wpseries(a^2)-wpseries(b^2))
	quad = quad*a^6*b^6
	
	squararray = False
    
elif switch == 'Q4w':
	g2 = 12*wp(v_)^2 - 2*wp2p(v_)
	#g2 = 12*wp(v_)^2 - 2*(wpp(2*v_) + wpp(v_))*wpp(v_)/(wp(v_) - wp(2*v_))
	g3 = expand(4*wp(v_)^3 - wpp(v_)^2 - g2*wp(v_))
	wpcoeff = [g2/20,g3/28]
	truncorder = 3 # CHECK ORDER
	for k in [4..truncorder]:
		c += [ expand( 3/(2*k+1)/(k-3) * sum([ wpcoeff[m-2]*wpcoeff[k-m-2] for m in [2..k-2]]) ) ]
	def wpseries(func,prime=False):
		out = z^(-2) + sum([wpcoeff[k-2]*z^(2*k-2) for k in [2..truncorder] ])
		if prime==False:
			return out.subs(z==func)
		else:
			return diff(out,z).subs(z==func)
	bb = wpseries(b^2)
	aa = wpseries(a^2)
	cc = wpseries(b^2 - a^2)
	
	quad =  b^6 * a^6 * wpseries(a^2,True) * ( (wp(ul(0,0)) - bb)*(wp(ul(0,1)) - bb) - (aa - bb)*(cc - bb) ) * ( (wp(ul(1,0)) - bb)*(wp(ul(1,1)) - bb) - (aa - bb)*(cc - bb) )
	quad += a^6 * b^6 * wpseries(b^2,True) * ( (wp(ul(0,0)) - aa)*(wp(ul(1,0)) - aa) - (bb - aa)*(cc - aa) ) * ( (wp(ul(0,1)) - aa)*(wp(ul(1,1)) - aa) - (bb - aa)*(cc - aa) )
	quad += - a^6 * b^6 * wpseries(a^2,True)*wpseries(b^2,True) * wpseries(b^2-a^2,True) * (aa - bb)
	
	squararray = False
    
elif switch == 'H1':
	quad = (1/b + 1/a + ul(1,1) - ul(0,0))*(1/b - 1/a + ul(0,1) - ul(1,0)) - (1/b^2 - 1/a^2)
elif switch == 'H3':
	if delta == 0:
		quad = 1/b*(ul(0,0)*ul(1,0) - ul(0,1)*ul(1,1)) - 1/a*(ul(0,0)*ul(0,1) - ul(1,0)*ul(1,1)) + delta*(a/b-b/a)
elif switch == 'H3sin':
	if delta == 0:
		quad = 1/b*sin(1/4*(ul(0,0)+ul(1,0)-ul(0,1)-ul(1,1))) - 1/a*sin(1/4*(ul(0,0)+ul(0,1)-ul(1,0)-ul(1,1)))
elif switch == 'H3tan':
	if delta == 0:
		quad = (1/b-1/a)*tan((ul(1,1)-ul(0,0))/4) - (1/b+1/a)*tan((ul(1,0)-ul(0,1))/4) + delta*(a/b-b/a)/(2*I)*exp(-I*(ul(0,0)+ul(1,0)+ul(0,1)+ul(1,1))/4)

elif switch == 'BSQ':
	quad = (1/a^3-1/b^3)/(1/a-1/b+ul(1,1)-ul(2,0)) - (1/a^3-1/b^3)/(1/a-1/b+ul(0,2)-ul(1,1)) - ul(0,1)*ul(1,2) + ul(1,0)*ul(2,1) + ul(2,2)*(1/a-1/b+ul(1,2)-ul(2,1)) + ul(0,0)*(1/a-1/b+ul(0,1)-ul(1,0)) - (2/a+1/b)*(ul(1,0)+ul(1,2)) + (1/a+2/b)*(ul(0,1)+ul(2,1))
	quad = quad/a/b
	
elif switch == 'BSQ3':
    #quad1 = ul(0,1,2) - ul(1,0,2) - (1/a - 1/b + ul(0,1,1) - ul(1,0,1))*ul(1,1,1) + 1/a*ul(0,1,1) - 1/b*ul(1,0,1)
    #quad2 = ul(0,1,3) - ul(1,0,3) + (1/a - 1/b + ul(0,1,1) - ul(1,0,1))*ul(0,0,1) + 1/b*ul(0,1,1) - 1/a*ul(1,0,1)
    #quad3 = (1/a - 1/b + ul(0,1,1) - ul(1,0,1))*(ul(1,1,2) - ul(0,0,3)) - (1/a + 1/b + ul(0,0,1)) * ( (1/a - 1/b + ul(0,1,1) - ul(1,0,1))*ul(1,1,1) - 1/a*ul(0,1,1) + 1/b*ul(1,0,1) ) - 1/a^2*(ul(1,0,1) - ul(0,0,1)) + 1/b^2*(ul(0,1,1) - ul(0,0,1))
    quad1 = ul(1,0,3) - ul(0,0)*ul(1,0) + ul(0,0,2)
    quad2 = ul(0,1,3) - ul(0,0)*ul(0,1) + ul(0,0,2)
    quad3 = ul(0,0)*ul(1,1) - ul(1,1,2) - ul(0,0,3) - (1/a-1/b)/(ul(1,0)-ul(0,1))
    quad = [quad1,quad2,quad3]

elif switch == 'GD4':
    quad1 = ul(0,1,2) - ul(1,0,2) - (1/a - 1/b + ul(0,1) - ul(1,0))*ul(1,1) + 1/a*ul(0,1) - 1/b*ul(1,0)
    quad2 = ul(0,1,3) - ul(1,0,3) + (1/a - 1/b + ul(0,1) - ul(1,0))*ul(0,0) + 1/b*ul(0,1) - 1/a*ul(1,0)

    quad3 = (1/a - 1/b + ul(0,1) - ul(1,0))*ul(1,1,2) - 1/a*ul(0,1,2) + 1/b*ul(1,0,2)
    quad3 += (1/a - 1/b + ul(-1,0) - ul(0,-1))*ul(-1,-1,3) + 1/b*ul(-1,0,3) - 1/a*ul(0,-1,3)
    quad3 += (1/a^4-1/b^4)/(1/a - 1/b + ul(-1,1) - ul(0,0)) - (1/a^4-1/b^4)/(1/a - 1/b + ul(0,0) - ul(1,-1))
    quad3 += (1/a+1/b)*(ul(0,1)*ul(-1,0) - ul(1,0)*ul(0,-1))
    quad3 += -( ul(0,1,2)*ul(-1,0) - ul(1,0,2)*ul(0,-1) + ul(0,1)*ul(-1,0,3) - ul(1,0)*ul(0,-1,3) )
    #quad3 += (1/a^2+1/a/b+1/b^2)*(ul(0,1)-ul(1,0)-ul(-1,0)+ul(0,-1)) - (1/a+1/b)*(ul(0,1,2) - ul(1,0,2) - ul(-1,0,3) + ul(0,-1,3))
    quad3 += (1/a^2+1/a/b+1/b^2)*(ul(0,1)-ul(1,0)-ul(-1,0)+ul(0,-1)) - (1/a+1/b)*((1/a - 1/b + ul(0,1) - ul(1,0))*ul(1,1) - 1/a*ul(0,1) + 1/b*ul(1,0)) - (1/a+1/b)*((1/a - 1/b + ul(-1,0) - ul(0,-1))*ul(-1,-1) + 1/b*ul(-1,0) - 1/a*ul(0,-1))
    
    quad = [quad2,quad1,quad3/a/b] #v2_-v3_-v1_^2+2*v1_1
    
elif switch == 'GD4skew':
    quad1 = ul(0,1,2) - ul(1,0,2) - (1/a - 1/b + ul(0,1) - ul(1,0))*(ul(1,1)-ul(0,0)) + (1/a+1/b)*(ul(0,1) - ul(1,0))
    quad2 = ul(0,1,3) - ul(1,0,3) - (1/a - 1/b + ul(0,1) - ul(1,0))*(ul(1,1)+ul(0,0)) + (1/a-1/b)*(ul(0,1) + ul(1,0))

    quad3 = (1/a - 1/b + ul(0,1) - ul(1,0))*(ul(1,1,2)+ul(1,1,3)) - 1/a*ul(0,1,2) + 1/b*ul(1,0,2) - 1/a*ul(0,1,3) + 1/b*ul(1,0,3)
    quad3 += (1/a - 1/b + ul(-1,0) - ul(0,-1))*(ul(-1,-1,2)-ul(-1,-1,3)) + 1/b*ul(-1,0,2) - 1/a*ul(0,-1,2) - 1/b*ul(-1,0,3) + 1/a*ul(0,-1,3)
    quad3 += 2*(1/a^4-1/b^4)/(1/a - 1/b + ul(-1,1) - ul(0,0)) - 2*(1/a^4-1/b^4)/(1/a - 1/b + ul(0,0) - ul(1,-1))
    quad3 += 2*(1/a+1/b)*(ul(0,1)*ul(-1,0) - ul(1,0)*ul(0,-1))
    quad3 += -( ul(0,1,3)*ul(-1,0) - ul(1,0,3)*ul(0,-1) + ul(0,1)*ul(-1,0,2) - ul(1,0)*ul(0,-1,2) ) -( ul(0,1,2)*ul(-1,0) - ul(1,0,2)*ul(0,-1) - ul(0,1)*ul(-1,0,3) + ul(1,0)*ul(0,-1,3) )
    quad3 += 2*(1/a^2+1/a/b+1/b^2)*(ul(0,1)-ul(1,0)-ul(-1,0)+ul(0,-1)) - (1/a+1/b)*(ul(0,1,3) - ul(1,0,3) - ul(-1,0,2) + ul(0,-1,2)) - (1/a+1/b)*(ul(0,1,2) - ul(1,0,2) + ul(-1,0,3) - ul(0,-1,3))
    
    quad = [quad1,quad3,quad2] #+ 2*quad2*(1/a^2/b+1/a/b^2)
    squararray = False
	
### END OF EQUATION LIST ###

############################

if not(type(quad) == type([])):
	quad = [quad]

### First series expansion
if double_eqn_exp:
	quadseries = [quad[i-1].series(b,numvars).truncate() for i in [1..components]]
else:
	quadseries = [quad[i-1].series(b,1).truncate() for i in [1..components]]
	
print(walltime(w))

### Check for negative order terms
for i in [-negdepth..-1]:
	for j in [1..components]:
		if not(expand(dertovar(quadseries[j-1].coefficient(b,i).series(a,numvars+1).truncate()) == 0)):
			textadd("Order " + str(i) + "? - Will check again later")
			latexadd(expand(dertovar(quadseries[j-1].coefficient(b,i).series(a,numvars+1).truncate())))

### Second series expansion
@parallel
def pexpand(i,component=1):
	if squararray:
		truncorder = numvars
	else:
		truncorder = numvars-i
	out = expand(dertovar( quadseries[component-1].coefficient(b,i).series(a,truncorder).truncate() ))
	return out

### Create array of coefficients of the double power series
textadd('Double series expansion:',vieweqnarray)

scaling = 1/miwaconst * diagonal_matrix([-(-1)^i*i for i in [1..numvars]]) 
eqnarray = [ matrix(0) for k in [0..components-1]]
matr = [[0 for j in [0..numvars-1]] for i in [0..numvars-1]]
#eqnarray = [[[0 for j in [0..numvars-1]] for i in [0..numvars-1]] for k in [0..components-1]]
for component in [1..components]:
	if double_eqn_exp:
		eqngen = pexpand([ (i,component) for i in [0..numvars-1]])
	else:
		eqngen = pexpand([ (0,component) ])
	for i in eqngen:
		#eqnarray[component-1][i[0][0][0]] = [i[1].coefficient(a,j) for j in [0..numvars-1]]
		matr[i[0][0][0]] = [i[1].coefficient(a,j) for j in [0..numvars-1]]
	eqnarray[component-1] = scaling * matrix(matr) * scaling

	latexadd(eqnarray[component-1],vieweqnarray)

### function to simplify trig and ellipitc expressions (for Q3-1 resp Q4)
def cleantrig(func, iterations=numvars/2):
	out = expand(func)
	if components == 1:
		out = out.subs( cos(v_)^2 == 1 - sin(v_)^2 )
		out = out.subs( cos(v_)^3 == (1 - sin(v_)^2)*cos(v_) )
		out = out.subs( cos(v_)^4 == (1 - sin(v_)^2)^2 )
		out = out.subs( cos(v_)^5 == (1 - sin(v_)^2)^2*cos(v_) )
		out = out.subs( cos(v_)^6 == (1 - sin(v_)^2)^3 )
		for i in [1..iterations]:
			out = out.subs(wp2p(2*v_) == diff(wp2p(v_)/wpp(v_)*(wp(v_)-wp(2*v_)) - wpp(v_),v_)/2 )
			out = out.subs(wpp(2*v_) == wp2p(v_)/wpp(v_)*(wp(v_)-wp(2*v_)) - wpp(v_) )
			out = out.subs(wp2p(v_)^5 == (4*wp(2*v_)+8*wp(v_))*wpp(v_)^2*wp2p(v_)^3 ) 
			out = out.subs(wp2p(v_)^4 == (4*wp(2*v_)+8*wp(v_))*wpp(v_)^2*wp2p(v_)^2 ) 
			out = out.subs(wp2p(v_)^3 == (4*wp(2*v_)+8*wp(v_))*wpp(v_)^2*wp2p(v_) ) 
			out = out.subs(wp2p(v_)^2 == (4*wp(2*v_)+8*wp(v_))*wpp(v_)^2 ) 
			#out = out.subs(wpp(2*v_) == wp2p(v_)/wpp(v_)*(wp(v_)-wp(2*v_)) - wpp(v_) )
			out = expand(out)
	return out

### Compile hierarchy of PDEs
pde = [{} for i in [1..components]]
for j in range(numvars)[1:]:
	for component in [1..components]:
		error = False
		if components == 1:
			derivative = 'v_' 
		else:
			derivative = 'v' + str(component) + '_' 
		if switch == 'BSQ': # and is_even(j):
			derivative += str(2) + str(j)
		elif switch == 'GD4' and component==3:
			derivative += str(2) + str(j)
		elif switch == 'GD4skew' and component==2:
			derivative += str(j)
		#elif switch == 'GD4' and j==2 and component==2:
		#	derivative += str(1)
		else:
			derivative += str(j)
		if eqnarray[component-1][0,j-1] == 0:
			sols = []
		else:
			sols = solve(eqnarray[component-1][0,j-1] == 0, eval(derivative))
	
		if eqnarray[component-1][0,j-1] == 0:
			error = True
		elif len(sols) == 0:
			error = True
	
		if error:
			print("No PDE for t" + str(j) + " (" + str(component) + ")")
			pde[component-1][j] = eval(derivative) == eval(derivative)
		else:
			pde[component-1][j] = expand(sols[0])
		
		pde[component-1][j] = cleantrig(pde[component-1][j])
#if not(components == 1):
#	print pde[1]
#	pde[1].update({2: v2_ == v3_+v1_^2-2*v1_1})
#	print pde[1]
#	latexadd(pde)

### Constraints
constraints = []
if switch == 'GD4':
	pde[2].update({2: v3_22 == 1/2*(pde[2][2].rhs() + v3_22)}) #botch
	constraints = [v2_ == v3_+v1_^2-2*v1_1]
if switch == 'GD4skew':
	constraints = [v3_ == v1_^2 + v1_1]
	var('k,l')
	pde[1].update({2: v2_22 == - 4*v1_11*v2_1 - 8*v1_1*v2_11 - v2_1111}) #botch
	if numvars >=3:
		pde[1].update({3: v2_3 == -1/2*v2_111 - 3*v1_1*v2_1}) #botch
	if numvars >=4:
		pde[1].update({4: v2_4 == 0}) #botch
	#if numvars >=5:
	#	pde[1].update({5: v2_5 == 0}) #botch
	
### List all differential consequences of hierarchy
allpde = []
for e in constraints:
	for index in indices(order,numvars):
		allpde += [vdiff(e,deri(index))]
if switch == 'GD4skew' and numvars >=3:
	allpde += [vdiff(pde[1][3],t2)]
	allpde += [vdiff(pde[1][3],t1,t2)]
for index in indices(order,numvars):
	for component in [1..components]:
		add = True
		### Determine lowest equation we can use to replace
		lowest = 1
		for i in reversed([2..numvars]): 
			if not(index[i-1] == 0):
				lowest = i
		if switch == 'GD4skew' and component == 2 and lowest == 2 and weight(index) > order-2:
				lowest = 1 #abort if order is too high (due to v2_22 = ...)
		if switch == 'BSQ' and weight(index) > order-2:
				lowest = 1 #abort if order is too high (due to v2_22 = ...)
		### Do replacement
		if not(lowest == 1):
			lowindex = copy(index)
			lowindex[lowest-1] += -1
			rhs = vdiff(pde[component-1][lowest].rhs(),deri(lowindex))
			lhs = vdiff(pde[component-1][lowest].lhs(),deri(lowindex))
			for eqn in allpde:
				rhs = expand(rhs.subs(eqn))
				if eqn.lhs() == lhs: #avoid multiple substitution
					add = False
			if add:
				allpde += [lhs == cleantrig(rhs).expand()]
#else:
#	if numvars >= 4:
#		pde[0][4] = (v_4 == 4/3*v_112 - 4*v_1*v_2)
#	if numvars >= 5:
#		pde[0][5] = (v_5 == 20/3*v_11^2 - 5/2*v_2^2 - 11/9*v_11111 + 5/12*v_122 - 20/3*v_1^3 + 10*v_1*v_111)
#	for index in indices(order,numvars):
#		### Determine lowest equation we can use to replace
#		lowest = 1
#		for i in reversed([3..numvars]):
#			if not(index[i-1] == 0):
#				lowest = i
#		if index[2-1] >= 2:
#			lowest = 2
#		### Do replacement
#		if not(lowest == 1):
#			lowindex = copy(index)
#			lowindex[lowest-1] += -1
#			if lowest == 2:
#				lowindex[2-1] += -1
#			lhs = vdiff(pde[0][lowest].lhs(),deri(lowindex))
#			rhs = vdiff(pde[0][lowest].rhs(),deri(lowindex))
#			for eqn in allpde:
#				rhs = expand(rhs.subs(eqn))
#			allpde += [lhs == rhs.simplify_full()]

### Replace time derivatives
def replace(func,iterations=components):
	out = func
	for i in [1..iterations]:
		out = out.subs(allpde)
	return expand(out)
			
### Simplify the hierarchy
for component in [1..components]:
	for i in pde[component-1]:
		pde[component-1][i] = pde[component-1][i].lhs() == cleantrig(replace(pde[component-1][i].rhs()))

### Output hierarchy
textadd('Simplified system of PDEs:')
for component in [1..components]:
	for i in pde[component-1]:
		latexadd(pde[component-1][i])

### Check again for negative order terms
for i in [-negdepth..-1]:
	for j in [1..components]:
		if not(replace(dertovar(quadseries[j-1].coefficient(b,i).series(a,numvars+1).truncate()) == 0)):
			textadd("ORDER " + str(i) + "!")
			latexadd(replace(dertovar(quadseries[j-1].coefficient(b,i).series(a,numvars+1).truncate())))
			
### Check if all coefficients of the double power series vanish on solutions
if double_eqn_exp:
	for component in [1..components]:
		if not(check(eqnarray[component-1])):
			warning = True

### Check commutativity
if checkcomm and numvars >= 5 and components == 1:
	for i in range(numvars)[1:-1]:
		for j in range(numvars)[2:]:
			if i < j:
				if not( i == 2 and switch == "BSQ" ):
					test = cleantrig(replace(vdiff(pde[i].rhs(), eval("t"+str(j)) ))) - cleantrig(replace(vdiff(pde[j].rhs(), eval("t"+str(i)) )))
				else:
					test = cleantrig(replace(vdiff(pde[2].rhs(), eval("t"+str(j)) ))) - cleantrig(replace(vdiff(pde[j].rhs(), t2,2 )))
				if not( test == 0 ):
					warning = True
					textadd("Nonzero commutator! " + str(i) + str(j))
					latexadd(expand( test ))
				else:
					textadd("PDEs for t" + str(i) + " and t" + str(j) + " commute.")
		
