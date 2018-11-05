#!/usr/bin/env sage

### SELECT LATTICE EQUATION ###

switch = 'H1'
#OPTIONS 
### with Lagrangian: 
# lin, Q1, Q2, Q3, Q4, H1, H3tan, BSQ, GD4skew
# for doubly infinite hierarchy: combine H3double, H3tan
### equation only:
# H3, Q4rat

# Specify parameter:
delta = 1
# Only relevant for Q1 and Q3
# Usually 0 or 1
# Comment out to keep as variable $\delta$


### DIMENSION OF MULTI-TIME

numvars = 3
# Recommend to start at 3. Computing time increases quicky for higher values.

includeeven = False
# For many lattice equation there will be trivial equations in the even times.
# Ignoring these speeds up computation.

### LAGRANGIANS

# Forget about the Lagrangians?
onlyequation = False

# Automatically find an exact form to symplify the Lagrangian 2-form?
autosimplify = True
# If False the 1-form given in c_list(switch) is used (if available).
# Not available for H3double

# Try to group terms rather than expand?
combineterms = False


### CHECKS

# Calculate and check full array of equations?
double_eqn_exp = True
# Not necessary if multidimensional consistency is known.
# Can be time-intensive.

#How deep to check for negative order terms:
negdepth = 3 

#To which extent to verify the EL equations:
elcheckdepth = numvars 
#-1 for none, 0 for quick, [1..numvars+1] to specify order of multi-indeces to check
#Note: The 'corner equations' involving 3 coefficients are not checked.

checkcomm = True #Explicitly verify commutativity?


### REQUESTED OUTPUT

# Let Tex built output pdf?
viewpdf = True

# Output the double series expansion of the quad equation?
vieweqnarray = False
vieweqnarray_iterate = False

# Output the series expansion of the Lagrangian?
viewldisc = False

# Output the Lagrangian 2-form before simplification?
viewlpluri = False

# Print the first few EL equations even if they are satisfied?
viewELeqs = False

############################


### NUMBER OF COMPONENTS ###
# Default: 1
multicompontent_list = {'GD4': 3, 'GD4skew': 3}

### REQUIRED ORDER ###
# Default: 2*numvars
order_list = {'BSQ': max(7,2*numvars),'GD4': 2*numvars+1, 'GD4skew': 2*numvars-1}

### IN WHICH COMPONENT DO WE EXPECT A SECOND ORDER LEADING EQUATION? ###
# Default: None
secondorder_list = {'BSQ': 1,'GD4': 3, 'GD4skew': 2}

### HIGHEST COMPONENTS IN LAGRANGIAN TWO-FORM ###
# Default: numvars
lagnumvars_list = {'GD4': numvars-1, 'GD4skew': numvars-1}

### MIWA CONSTANT ###
# Default: -2
miwa_list = {'BSQ': -3, 'GD4': 1, 'GD4skew': 1}

### SQUARE ARRAY or TRIANGLE ###
# Expand equations in square array or just a triangle?
# Default: square array
squarearray_list = {'Q4rat': False, 'Q4': False, 'GD4skew': False}

### SUGGEST EQUATIONS FOR CONTINUUM LIMIT ###
# Default: None
# Equations to be used from the beginning:
def constraint_list(switch):
	if switch == 'H3double': return [v_1_1 == sin(v__)]#2*sin(v_/2)*cos(v_/2)]
	elif switch == 'GD4': return [v2_ == v3_ + v1_^2 - 2*v1_1], 
	elif switch == 'GD4skew': return [v3_ == v1_^2 + v1_1]#, v2_22 == - 4*v1_11*v2_1 - 8*v1_1*v2_11 - v2_1111]
	else: return []
# Equations to replace the pde found at the corresponding order:
var('v_4,v_5,v_6,v_5_,v__5') #need strategy to avoid un-initialized variables
def pde_list(switch):
	if switch == 'H3double': 
		return [{3: v_3_ == v_111_ + 1/2*v_1_^3, 
			5: v_5_ == 3/8*v_1_^5 + 5/2*v_1_*v_11_^2 + 5/2*v_1_^2*v_111_ + v_11111_,
			numvars+3: v__3 == v__111 + 1/2*v__1^3, 
			numvars+5: v__5 == 3/8*v__1^5 + 5/2*v__1*v__11^2 + 5/2*v__1^2*v__111 + v__11111}]
	elif switch == 'BSQ': 
		return [{3: v_3 == 0, 4: v_4 == 3*v_112 - 6*v_1*v_2, 5: v_5 == -15*v_1^3 + 135/4*v_11^2 + 45*v_1*v_111 - 15/4*v_2^2 - 9*v_11111, 6: v_6 == 0}]
	elif switch == 'GD4skew': 
		var('v2_4') #need strategy to avoid un-initialized variables
		return [{},{2: v2_22 == - 4*v1_11*v2_1 - 8*v1_1*v2_11 - v2_1111, 3: v2_3 == -1/2*v2_111 - 3*v1_1*v2_1, 4:v2_4 == 0},{}]
	else: 
		return [{},{},{}] #needs at least as many entries as components

### SPECIFY SYMPLIFYING 1-FORM ###
def c_list(switch):
	if switch == 'H1':
		if includeeven:
			return [0,0,0,4/3*v_1*v_12-2/3*v_11*v_2,10/3*v_1^2*v_11-4/9*v_11*v_111-1/9*v_1*v_1111+10/9*v_1*v_13+5/12*v_1*v_22-5/9*v_11*v_3]
		else:
			return [0,0,0,0,10/3*v_1^2*v_11-4/9*v_11*v_111-1/9*v_1*v_1111+10/9*v_1*v_13-5/9*v_11*v_3]
	if switch == 'H3double':
		if includeeven:
			return [[0, 1/4*v_1_^2 + 1/4*v__*v_2_, -1/4*v_1_*v_11_ + 3/8*v_1_*v_2_ + 1/4*v__*v_3_],[0, -1/4*v_1_^2 - 1/4*v__*v_2_, +1/4*v_1_*v_11_ - 3/8*v_1_*v_2_ - 1/4*v__*v_3_]]
		else:
			return [[0, 1/4*v_1_^2, -1/4*v_1_*v_11_ + 1/4*v__*v_3_, -1/6*v_11_^2 + 1/3*v_1_*v_3_, -5/36*v_1_^3*v_11_ + 1/9*v_11_*v_111_ + 1/36*v_1_*v_1111_ - 5/18*v_1_*v_13_ + 5/36*v_11_*v_3_ + 1/4*v__*v_5_]
,[0, -1/4*v__1^2, 1/4*v__1*v__11 - 1/4*v__*v__3, 1/6*v__11^2 - 1/3*v__1*v__3, 5/36*v__1^3*v__11 - 1/9*v__11*v__111 - 1/36*v__1*v__1111 + 5/18*v__1*v__13 - 5/36*v__11*v__3 - 1/4*v__*v__5]]


#########################################

### LIST OF EQUATIONS AND LAGRANGIANS ###

def get_equation(switch,lagrangian=False):
	quad = 0*ul(0,0)
	L = 0*ul(0,0)
	###	
	###
	if switch == 'lin':
		quad = (1/b + 1/a)*(ul(0,1) - ul(1,0)) + (1/b - 1/a)*(ul(1,1) - ul(0,0))
		###
		if lagrangian:
			L = 1/2*(ul(0,1) - ul(1,0)) * (ul(1,1) - ul(0,0)) + 1/2*(a^2 - b^2) * (quot12)^2
	###	
	###
	if switch == 'Q1':
		quad = 1/b^2*(ul(0,0)-ul(0,1))*(ul(1,0)-ul(1,1)) - 1/a^2*(ul(0,0)-ul(1,0))*(ul(0,1)-ul(1,1)) + delta^2*(a^2-b^2)
		###
		if lagrangian:
			if delta == 0:
				L = a^2*logseries(diff1/(-2*a*v_1), numvars-2) - b^2*logseries(diff2/(-2*b*v_1), numvars-2) - (a^2 - b^2)*logseries(quot12/(-2*v_1), 2*numvars-3)
			else:
				L1 = (-diff1+a^2*delta)*logseries( (-diff1+a^2*delta)/(2*a*v_1), numvars-1) - (-diff1-a^2*delta)*logseries( (-diff1-a^2*delta)/(2*a*v_1), numvars-1)
				L2 = L1.subs(a == b)
				L12 =(diff12+(a^2-b^2)*delta)*logseries( (quot12+(a+b)*delta)/(-2*v_1), 2*numvars-2) - (diff12-(a^2-b^2)*delta)*logseries( (quot12-(a+b)*delta)/(-2*v_1), 2*numvars-2)     
				L = L1 - L2 - L12
	###	
	###
	if switch == 'Q2':
		quad = 1/b^2*(ul(0,0)^2-ul(0,1)^2)*(ul(1,0)^2-ul(1,1)^2) - 1/a^2*(ul(0,0)^2-ul(1,0)^2)*(ul(0,1)^2-ul(1,1)^2) + (a^2-b^2)*(ul(0,0)^2+ul(1,0)^2+ul(0,1)^2+ul(1,1)^2) - (a^2-b^2)*(a^4-a^2*b^2+b^4)
		###
		if lagrangian:
			L1 = (sum1+a^2)*logseries((sum1+a^2)/(2*v_), numvars-1) - (sum1-a^2)*logseries((sum1-a^2)/(2*v_), numvars-1) + (-diff1+a^2)*logseries( (-diff1+a^2)/(2*a*v_1), numvars-1) - (-diff1-a^2)*logseries( (-diff1-a^2)/(2*a*v_1), numvars-1)
			L2 = L1.subs(a == b)
			L12 = (sum12+a^2-b^2)*logseries((sum12+(a^2-b^2))/(2*v_), 2*numvars-1) - (sum12-a^2+b^2)*logseries((sum12-(a^2-b^2))/(2*v_), 2*numvars-1) + (diff12+a^2-b^2)*logseries( (quot12+a+b)/(-2*v_1), 2*numvars-2) - (diff12-(a^2-b^2))*logseries( (quot12-a-b)/(-2*v_1), 2*numvars-2)     
			L = L1 - L2 - L12
	###	
	###
	if switch == 'Q3':
		if delta == 0:
			quad = -sin(a^2-b^2)*(exp(I*ul(0,0)+I*ul(1,1))+exp(I*ul(1,0)+I*ul(0,1))) + sin(a^2)*(exp(I*ul(0,0)+I*ul(1,0))+exp(I*ul(0,1)+I*ul(1,1))) - sin(b^2)*(exp(I*ul(0,0)+I*ul(0,1))+exp(I*ul(1,0)+I*ul(1,1)))
		else:
			quad = -sin(a^2-b^2)*( cos(ul(0,0))*cos(ul(1,1)) + cos(ul(1,0))*cos(ul(0,1)) ) + sin(a^2)*( cos(ul(0,0))*cos(ul(1,0)) + cos(ul(0,1))*cos(ul(1,1)) ) - sin(b^2)*( cos(ul(0,0))*cos(ul(0,1)) + cos(ul(1,0))*cos(ul(1,1)) ) + delta*sin(a^2-b^2)*sin(a^2)*sin(b^2)
		quad = quad/(a*b)^2
		###
		if lagrangian:
			L1 = a^2*diff1 - dilogseries(I*(diff1+a^2),numvars) + dilogseries(I*(diff1-a^2),numvars)
			if not(delta == 0):
				L1 += delta*( a^2*sum1 - dilogseries(I*(sum1+a^2),numvars) + dilogseries(I*(sum1-a^2),numvars) )
			L2 = L1.subs(a == b)
			L12 = (a^2-b^2)*diff12 - dilogseries(I*(diff12+a^2-b^2),2*numvars-1) + dilogseries(I*(diff12-a^2+b^2),2*numvars-1)
			if not(delta == 0):
				L12 += delta*( (a^2-b^2)*sum12 - dilogseries(I*(sum12+a^2-b^2),2*numvars-1) + dilogseries(I*(sum12-a^2+b^2),2*numvars-1) )
			L = (L1 - L2 - L12)/(2*I)
	###	
	###
	if switch == 'Q4rat':
		var('g2','g3')
		global g2
		global g3
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
	###	
	###
	if switch == 'Q4':
		g2 = 12*wp(v_)^2 - 2*wp2p(v_)
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
		###
		if lagrangian:
			L1 = a^2*logseries(wsigma(sum1)/wsigma(2*v_),numvars-2*0)
			L1 += +a*( (logseries(diff1/(-2*a*v_1)+a/(-2*v_1),numvars-1) + log(-2*v_1))*(diff1/2/a+a/2) )
			L1 += -a*( (logseries(diff1/(-2*a*v_1)-a/(-2*v_1),numvars-1) + log(-2*v_1))*(diff1/2/a-a/2) )
			L2 = L1.subs(a == b)
			L12 = ( (a^2-b^2)*logseries(wsigma(sum12)/wsigma(2*v_),2*numvars-3) - 1/6*wp(sum12)*(a^2-b^2)^3 )
			L12 += (logseries(quot12/(-2*v_1)+(a+b)/(-2*v_1),2*numvars-2) + log(-2*v_1)) * (diff12+(a^2-b^2))/2
			L12 += -(logseries(quot12/(-2*v_1)-(a+b)/(-2*v_1),2*numvars-2) + log(-2*v_1)) * (diff12-(a^2-b^2))/2
			L12 += killpowers(expand( -1 * g2/2400*( 10*diff12^4*(a^2-b^2) + 20*diff12^2*(a^2-b^2)^3 + 2*(a^2-b^2)^5 ) ))
			L12 += killpowers(expand( -1 * g3/11760*( 14*diff12^6*(a^2-b^2) ) )) #+ 280*diff12^4*(a^2-b^2)^3) ))
			L = L1 - L2 - L12	
	###	
	###
	if switch == 'H1':
		quad = (1/b + 1/a + ul(1,1) - ul(0,0))*(1/b - 1/a + ul(0,1) - ul(1,0)) - (1/b^2 - 1/a^2)
		###
		if lagrangian:
			L = 1/2*(-diff12) * (dertovar( (ul(1,1) - ul(0,0)).series(a,numvars+1).truncate().series(b,numvars).truncate() ) + 2/b + 2/a) + (1/a^2-1/b^2) * logseries(1-a*b*quot12, numvars, numvars+2)
	###
	###
	if switch == 'H3' and delta == 0:
		quad = 1/b*(ul(0,0)*ul(1,0) - ul(0,1)*ul(1,1)) - 1/a*(ul(0,0)*ul(0,1) - ul(1,0)*ul(1,1))
		L = 0
	###
	###
	if switch == 'H3tan' and delta == 0:
		quad = (1/b-1/a)*tan((ul(1,1)-ul(0,0))/4) - (1/b+1/a)*tan((ul(1,0)-ul(0,1))/4) + delta*(a/b-b/a)/(2*I)*exp(-I*(ul(0,0)+ul(1,0)+ul(0,1)+ul(1,1))/4)
		###
		if lagrangian:
			quot12ab = killpowers(expand( (a+b)*quot12 )) 
			var('x','y','z')
			L1 = 1/8*diff1^2
			L2 = L1.subs(a == b)
			integrand = arctan(x * tan(y/4)).series(y,2*numvars-2).truncate()
			L12 = expand( integrate(integrand, y).subs(x == z/y) )
			ypowers = [1]
			zpowers = [1]
			for i in [1..2*numvars-2]:
				ypowers += [killpowers(expand(ypowers[i-1]*diff12))]
				zpowers += [killpowers(expand(zpowers[i-1]*quot12ab))]
			for i in reversed([1..2*numvars-2]):
				L12 = L12.subs(z^i == zpowers[i])
				L12 = L12.subs(y^i == ypowers[i])
			L = expand(killpowers(L1 - L2 - L12))
	###	
	###
	if switch == 'H3double':
		### EXPERIMENTAL
		quad = sin((ul(0,0)+ul(1,0)+ul(0,1)+ul(1,1))/4) - 1/(a*b)*sin((ul(0,0)-ul(1,0)-ul(0,1)+ul(1,1))/4)
		###
		if lagrangian:
			var('x','y','z')
			L1 = 1/8*diff1^2 
			L2 = 1/8*sum2^2 
			integrand = arctan( -tan(z/4) + x ).series(x,2*numvars-2).truncate()
			integrand = integrand.subs(x == tan(z/4)*2*a*b/(a*b-1))
			L12 = expand( integrate(integrand, z) )
			#L12 = L12.subs( tan(z/4) == sin(z/4)/cos(z/4) )
			L12 = L12.subs( sin(z/4) == sqrt( (1-cos(z/2))/2 ) )
			#L12 = L12.subs( cos(z/4) == sqrt( (1+cos(z/2))/2 ) )
			zpowers = [1]
			for i in [1..2*numvars-2]:
				zpowers += [killpowers(expand( zpowers[i-1]*sum12 ))]
			for i in reversed([1..2*numvars-2]):
				L12 = L12.subs(z^i == zpowers[i])
			L = expand(killpowers(L1 - L2 - L12))
			L += -1/8*( ul(1,0)^2-ul(0,0)^2 + ul(1,1)^2-ul(0,1)^2 ) # d of a non-isotropic 1-form
		
		#quad = quad.subs(a == z).subs(b == a).subs(z == b)
		#L = L.subs(a == z).subs(b == a).subs(z == b)
	###	
	###
	if switch == 'BSQ':
		quad = (1/a^3-1/b^3)/(1/a-1/b+ul(1,1)-ul(2,0)) - (1/a^3-1/b^3)/(1/a-1/b+ul(0,2)-ul(1,1)) - ul(0,1)*ul(1,2) + ul(1,0)*ul(2,1) + ul(2,2)*(1/a-1/b+ul(1,2)-ul(2,1)) + ul(0,0)*(1/a-1/b+ul(0,1)-ul(1,0)) - (2/a+1/b)*(ul(1,0)+ul(1,2)) + (1/a+2/b)*(ul(0,1)+ul(2,1))
		quad = quad/a/b
		###
		if lagrangian:
			L = (1/a^3 - 1/b^3)*logseries(1 + a*b*quot12, numvars, numvars+3) - (1/a^2+1/a/b+1/b^2)*(ul(0,1)-ul(1,0)) + (1/a+1/b)*(ul(0,1)-ul(1,0))*ul(1,1) + (1/a-1/b+ul(0,1)-ul(1,0))*ul(0,0)*ul(1,1) - 1/a*ul(0,0)*ul(0,1) + 1/b*ul(0,0)*ul(1,0)
			L += 1/2*( ul(0,0)*ul(1,0)/a - ul(0,1)*ul(1,1)/a + ul(1,0)*ul(1,1)/b - ul(0,0)*ul(0,1)/b )
			L += 1/3*( ul(1,0)^3 - ul(0,1)^3 )
			L += 1/4*( (ul(1,0)^2-ul(0,0)^2)/a + (ul(1,1)^2-ul(1,0)^2)/b + (ul(0,1)^2-ul(1,1)^2)/a + (ul(0,0)^2-ul(0,1)^2)/b )
			L = dertovar(expand(L).series(a,numvars+1).truncate().series(b,numvars).truncate())
	###	
	###
	if switch == 'GD4':
		quad1 = ul(0,1,2) - ul(1,0,2) - (1/a - 1/b + ul(0,1) - ul(1,0))*ul(1,1) + 1/a*ul(0,1) - 1/b*ul(1,0)
		quad2 = ul(0,1,3) - ul(1,0,3) + (1/a - 1/b + ul(0,1) - ul(1,0))*ul(0,0) + 1/b*ul(0,1) - 1/a*ul(1,0)
		quad3 = (1/a - 1/b + ul(0,1) - ul(1,0))*ul(1,1,2) - 1/a*ul(0,1,2) + 1/b*ul(1,0,2)
		quad3 += (1/a - 1/b + ul(-1,0) - ul(0,-1))*ul(-1,-1,3) + 1/b*ul(-1,0,3) - 1/a*ul(0,-1,3)
		quad3 += (1/a^4-1/b^4)/(1/a - 1/b + ul(-1,1) - ul(0,0)) - (1/a^4-1/b^4)/(1/a - 1/b + ul(0,0) - ul(1,-1))
		quad3 += (1/a+1/b)*(ul(0,1)*ul(-1,0) - ul(1,0)*ul(0,-1))
		quad3 += -( ul(0,1,2)*ul(-1,0) - ul(1,0,2)*ul(0,-1) + ul(0,1)*ul(-1,0,3) - ul(1,0)*ul(0,-1,3) )
		quad3 += (1/a^2+1/a/b+1/b^2)*(ul(0,1)-ul(1,0)-ul(-1,0)+ul(0,-1)) - (1/a+1/b)*((1/a - 1/b + ul(0,1) - ul(1,0))*ul(1,1) - 1/a*ul(0,1) + 1/b*ul(1,0)) - (1/a+1/b)*((1/a - 1/b + ul(-1,0) - ul(0,-1))*ul(-1,-1) + 1/b*ul(-1,0) - 1/a*ul(0,-1))
		quad = [quad2,quad1,quad3/a/b]
	###	
	###
	if switch == 'GD4skew':
		quad1 = ul(0,1,2) - ul(1,0,2) - (1/a - 1/b + ul(0,1) - ul(1,0))*(ul(1,1)-ul(0,0)) + (1/a+1/b)*(ul(0,1) - ul(1,0))
		quad2 = ul(0,1,3) - ul(1,0,3) - (1/a - 1/b + ul(0,1) - ul(1,0))*(ul(1,1)+ul(0,0)) + (1/a-1/b)*(ul(0,1) + ul(1,0))
		def v(n,m): return (ul(n,m,2) + ul(n,m,3))/2
		def w(n,m): return (ul(n,m,2) - ul(n,m,3))/2
		quad3 = (1/a - 1/b + ul(0,1) - ul(1,0))*v(1,1) - 1/a*v(0,1) + 1/b*v(1,0)
		quad3 += (1/a - 1/b + ul(-1,0) - ul(0,-1))*w(-1,-1) + 1/b*w(-1,0) - 1/a*w(0,-1)
		quad3 += -( v(0,1)*ul(-1,0) - v(1,0)*ul(0,-1) + ul(0,1)*w(-1,0) - ul(1,0)*w(0,-1) ) 
		quad3 += - (1/a+1/b)*(v(0,1) - v(1,0) - w(-1,0) + w(0,-1))
		quad3 += (1/a^4-1/b^4)/(1/a - 1/b + ul(-1,1) - ul(0,0)) - (1/a^4-1/b^4)/(1/a - 1/b + ul(0,0) - ul(1,-1))
		quad3 += (1/a+1/b)*(ul(0,1)*ul(-1,0) - ul(1,0)*ul(0,-1))
		quad3 += (1/a^2+1/a/b+1/b^2)*(ul(0,1)-ul(1,0)-ul(-1,0)+ul(0,-1))
		quad = [quad1,quad3,quad2]
		###
		if lagrangian:
			var('k,l,m')
			L = -(1/a^4 - 1/b^4)*logseries(1 + a*b*quot12, numvars, numvars+3) + (1/a^3+1/a^2/b+1/a/b^2+1/b^3)*(ul(0,1)-ul(1,0)) 
			L += - (1/a^2+1/a/b+1/b^2)*(ul(0,1)-ul(1,0))*ul(1,1) + (1/a+1/b)*(ul(0,1)-ul(1,0))*v(1,1)
			L += ul(0,0)*( (1/a-1/b+ul(0,1)-ul(1,0))*v(1,1) - 1/a*v(0,1) + 1/b*v(1,0))
			L += ( (1/a+1/b)*ul(0,0)-w(0,0) )*( v(0,1) - v(1,0) - (1/a-1/b+ul(0,1)-ul(1,0))*ul(1,1) + 1/a*ul(0,1) - 1/b*ul(1,0))
			###
			L += -1/2*( ul(0,0)*ul(1,0)/a^2 - ul(0,1)*ul(1,1)/a^2 + ul(1,0)*ul(1,1)/b^2 - ul(0,0)*ul(0,1)/b^2 )
			L += -1/4*( (ul(1,0)^2-ul(0,0)^2)/a^2 + (ul(1,1)^2-ul(1,0)^2)/b^2 + (ul(0,1)^2-ul(1,1)^2)/a^2 + (ul(0,0)^2-ul(0,1)^2)/b^2 )
			###
			L += -1/6*( ul(0,0)*ul(1,0)^2/a - ul(0,1)*ul(1,1)^2/a + ul(1,0)*ul(1,1)^2/b - ul(0,0)*ul(0,1)^2/b )
			L += -1/6*( ul(0,0)^2*ul(1,0)/a - ul(0,1)^2*ul(1,1)/a + ul(1,0)^2*ul(1,1)/b - ul(0,0)^2*ul(0,1)/b )
			L +=  1/2*( ul(0,0)*ul(1,0,2)/a - ul(0,1)*ul(1,1,2)/a + ul(1,0)*ul(1,1,2)/b - ul(0,0)*ul(0,1,2)/b )
			L +=  1/2*( ul(0,0,2)*ul(1,0)/a - ul(0,1,2)*ul(1,1)/a + ul(1,0,2)*ul(1,1)/b - ul(0,0,2)*ul(0,1)/b )
			L +=  1/2*( (-ul(0,0)*ul(0,0,2))/a + (-ul(1,0)*ul(1,0,2))/b + (ul(0,1)*ul(0,1,2))/a + (ul(0,0)*ul(0,0,2))/b )
			L +=  1/2*( ul(0,0)*ul(1,0,3)/a - ul(0,1)*ul(1,1,3)/a + ul(1,0)*ul(1,1,3)/b - ul(0,0)*ul(0,1,3)/b )
			###
			L += +1/8*( ul(0,0,3)*ul(1,0,3) - ul(0,1,3)*ul(1,1,3) + ul(1,0,3)*ul(1,1,3) - ul(0,0,3)*ul(0,1,3) )
			L += -1/8*( ul(0,0,2)*ul(1,0,2) - ul(0,1,2)*ul(1,1,2) + ul(1,0,2)*ul(1,1,2) - ul(0,0,2)*ul(0,1,2) )
			L += -1/8*( ul(0,0,2)*ul(1,0,3) - ul(0,1,2)*ul(1,1,3) + ul(1,0,2)*ul(1,1,3) - ul(0,0,2)*ul(0,1,3) )
			L += -1/8*( ul(0,0,3)*ul(1,0,2) - ul(0,1,3)*ul(1,1,2) + ul(1,0,3)*ul(1,1,2) - ul(0,0,3)*ul(0,1,2) )
			L +=  1/4*( ul(0,0,2)*ul(1,0)^2 - ul(0,1,2)*ul(1,1)^2 + ul(1,0,2)*ul(1,1)^2 - ul(0,0,2)*ul(0,1)^2 )
			L +=  1/4*( ul(0,0)^2*ul(1,0,2) - ul(0,1)^2*ul(1,1,2) + ul(1,0)^2*ul(1,1,2) - ul(0,0)^2*ul(0,1,2) )
			###
			L = dertovar(expand(L).series(a,lagnumvars+1).truncate().series(b,lagnumvars).truncate())
	###	
	###
	if lagrangian:
		return [quad,L]
	else:
		return [quad]
