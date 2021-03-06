#!/usr/bin/env sage

### CONSTRUCT SOME BUILDING BLOCKS ###

def logseries(f,order,truncorder=lagnumvars):
	g = expand(f-1)
	powers1 = [1]
	for k in [1..order]:
		pow = powers1[k-1]*g
		pow = killpowers(pow,truncorder)
		powers1 += [pow] 
	out = expand(sum( (-1)^(i+1)/i * powers1[i] for i in [1..order]))
	return out

def dilogseries(func,truncorder): # expansion of Li2(exp(func))
	f0 = func.series(a,1).truncate().series(b,1).truncate()
	if f0 == 0:
		c = [dilog(1),1-log(-func)]
		for coeff in ( (-log((1-exp(v_))/v_)).series(v_,truncorder+1).coefficients() )[1:]:
			c += [coeff[0]*factorial(coeff[1])]
	else:
		c= [dilog(exp(f0))] + [diff(target/2 - log(sin(target/2/I)),target,i).subs(target == f0) for i in [0..truncorder]]
	out = 0
	powers = [1]
	for i in [1..truncorder]:
		powers += [killpowers( powers[i-1]*(func-f0), numvars).simplify_rational()]
	out = killpowers( expand( sum( c[i]/factorial(i)*powers[i] for i in [0..truncorder]) ),truncorder)
	return out
	
### Ui-U and Ui+U
diff1 = dertovar( (ul(1,0)-ul(0,0)).series(a,lagnumvars+1).truncate() )
diff2 = dertovar( (ul(0,1)-ul(0,0)).series(b,lagnumvars+1).truncate() )
sum1 = dertovar( (ul(1,0)+ul(0,0)).series(a,lagnumvars+1).truncate() )
sum2 = dertovar( (ul(0,1)+ul(0,0)).series(b,lagnumvars+1).truncate() )

### U1-U2
diff12 = expand( dertovar( ul(1,0).series(a,lagnumvars+1).truncate() - ul(0,1).series(b,lagnumvars+1).truncate() )) 

### (U1-U2)/(a-b)
quot12 = expand( dertovar( (ul(1,0) - ul(0,0)).series(a,numvars+3).truncate() ) )
var("a2")
for e in reversed([1..numvars+2]):
    quot12 = quot12.subs(a^e == sum( a2^i*b^(e-i-1) for i in [0..e-1] ) )
quot12 = expand((quot12.subs(a2 == a)))

### U1+U2
sum12 = expand( dertovar( (ul(1,0) + ul(0,1)).series(a,lagnumvars+2).truncate().series(b,lagnumvars+2).truncate() ) )

############################

### List constraints in allpde
allpde = []

def replace(func,iterations=components*dt):
	out = func
	for i in [1..iterations]:
		out = out.subs(allpde)
	return expand(out)
	
for e in constraints:
	for index in indices(order - dt,dtnumvars):
		ediff = vdiff(e,deri(index))
		if not( ediff.lhs() in [eqn.lhs() for eqn in allpde] ):  #avoid multiple substitution
			allpde += [replace(ediff)]
	
### Get equation (and Lagrangian) from input file
if onlyequation:
	quad = get_equation(switch,False)[0]
else:
	pair = get_equation(switch,True)
	quad = pair[0]
	L = pair[1]

if not(type(quad) == type([])):
	quad = [quad]

### First series expansion
if double_eqn_exp:
	quadseries = [quad[i-1].series(b,numvars).truncate() for i in [1..components]]
else:
	quadseries = [quad[i-1].series(b,1).truncate() for i in [1..components]]

### Check for negative order terms
for i in [-negdepth..-1]:
	for j in [1..components]:
		if not(replace(dertovar(quadseries[j-1].coefficient(b,i).series(a,numvars+1).truncate()) == 0)):
			textadd("Order " + str(i) + "? - Will check again later")
			latexadd(replace(dertovar(quadseries[j-1].coefficient(b,i).series(a,numvars+1).truncate())))

### Second series expansion
@parallel
def pexpand(i,component=1):
	if squarearray:
		truncorder = numvars
	else:
		truncorder = numvars-i
	out = replace(dertovar( quadseries[component-1].coefficient(b,i).series(a,truncorder).truncate() ))
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
		matr[i[0][0][0]] = [replace(i[1].coefficient(a,j)) for j in [0..numvars-1]]
	eqnarray[component-1] = scaling * matrix(matr) * scaling

	latexadd(eqnarray[component-1],vieweqnarray)

### Compile hierarchy of PDEs
pde = [{} for i in [1..components]]
	
for component in [1..components]:
	if not(doubletime):
		for j in range(numvars)[1:]:
			error = False
			if components == 1:
				derivative = 'v_' 
			else:
				derivative = 'v' + str(component) + '_' 
			if component == secondorder and j == 2:
				derivative += str(2) + str(j)
			else:
				derivative += str(j)
			if eqnarray[component-1][0,j-1] == 0:
				sols = []
			else:
				sols = solve(eqnarray[component-1][0,j-1] == 0, eval(derivative))
			if len(sols) == 0:
				error = True
	
			if error:
				#textadd("No PDE for t" + str(j) + " (" + str(component) + ")")
				pde[component-1][j] = eval(derivative) == eval(derivative)
			else:
				pde[component-1][j] = ( sols[0].lhs() == cleantrig(replace(sols[0].rhs())) )
	
	pde[component-1].update(imposedpdes[component-1])

	### List all differential consequences of hierarchy
	for index in indices(order,dtnumvars):
	#for component in [1..components]:
		add = True
		### Determine lowest equation we can use to replace
		lowest = 1
		for i in reversed([3..numvars] + [numvars+2..dtnumvars]): 
			if index[i-1] > 0:
				lowest = i
		if component == secondorder:
			if index[1] > 1:
				lowest = 2
		else:
			if index[1] > 0:
				lowest = 2
		### Do replacement
		if not(lowest == 1 or lowest == numvars + 1):
			lowindex = copy(index)
			if component == secondorder and lowest == 2:
				lowindex[lowest-1] += -2
			else:
				lowindex[lowest-1] += -1
			lhs = vdiff(pde[component-1][lowest].lhs(),deri(lowindex))
			if not( lhs in [eqn.lhs() for eqn in allpde] ):  #avoid multiple substitution
				rhs = replace(vdiff(pde[component-1][lowest].rhs(),deri(lowindex)))
				allpde += [lhs == expand(cleantrig(rhs))]
	
	allpde = [ e.lhs() == replace(e.rhs()) for e in allpde] # catch incomplete replacements
	
	if vieweqnarray_iterate:
		eqnarray = [replace(i) for i in eqnarray]
		textadd("Iteration " + str(component))
		for i in eqnarray:
			latexadd(i,vieweqnarray_iterate)

### Simplify the hierarchy
for component in [1..components]:
	for i in pde[component-1]:
		pde[component-1][i] = pde[component-1][i].lhs() == cleantrig(replace(pde[component-1][i].rhs()))

### Output hierarchy
textadd('Simplified system of PDEs:')
for c in constraints:
	latexadd(c)
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
if not(doubletime) and checkcomm and components == 1:
	for component in [1..components]:
		for i in range(numvars)[1:-1]:
			for j in range(numvars)[2:]:
				if i < j:
					if not( component == secondorder and i == 2):
						test = cleantrig(replace(vdiff(pde[component-1][i].rhs(), eval("t"+str(j)) ))) - cleantrig(replace(vdiff(pde[0][j].rhs(), eval("t"+str(i)) )))
					else:
						test = cleantrig(replace(vdiff(pde[component-1][2].rhs(), eval("t"+str(j)) ))) - cleantrig(replace(vdiff(pde[0][j].rhs(), t2,2 )))
					if not( test == 0 ):
						warning = True
						textadd("Nonzero commutator! " + str(i) + str(j))
						latexadd(expand( test ))
					else:
						textadd("PDEs for t" + str(i) + " and t" + str(j) + " commute.")
		
