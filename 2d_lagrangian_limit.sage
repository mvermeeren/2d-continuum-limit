#!/usr/bin/env sage

print("Calculating Lagrangian...")
print("This could take a while...")

scaling = 1/miwaconst * diagonal_matrix([-(-1)^i*i for i in [1..lagnumvars]]) 

### CONSTRUCT SOME BUILDING BLOCKS ###

def logseries(f,order,truncorder=lagnumvars):
    g = (f-1).simplify_rational()
    powers1 = [1]
    for k in [1..order]:
        pow = powers1[k-1]*g
        pow = killpowers(pow,truncorder).simplify_rational()
        powers1 += [pow] 
    out = killpowers( expand( sum( (-1)^(i+1)/i * powers1[i] for i in [1..order]) ),truncorder)
    return out

### Ui-U and Ui+U
diff1 = dertovar( (ul(1,0)-ul(0,0)).series(a,lagnumvars+1).truncate() )
diff2 = diff1.subs(a == b)
sum1 = dertovar( (ul(1,0)+ul(0,0)).series(a,lagnumvars+1).truncate() )
sum2 = diff1.subs(a == b)

### U1-U2
diff12 = expand( dertovar( (ul(1,0) - ul(0,1)).series(a,lagnumvars+1).truncate().series(b,lagnumvars+1).truncate() ) )

### (U1-U2)/(a-b)
#adapt order to method?
quot12 = expand( dertovar( (ul(1,0) - ul(0,0)).series(a,numvars+3).truncate() ) )
var("a2")
for e in reversed([1..numvars+2]):
    quot12 = quot12.subs(a^e == sum( a2^i*b^(e-i-1) for i in [0..e-1] ) )
quot12 = expand((quot12.subs(a2 == a)))

### U1+U2
sum12 = expand( dertovar( (ul(1,0) + ul(0,1)).series(a,lagnumvars+2).truncate().series(b,lagnumvars+2).truncate() ) )

#def ulintegral(func):
#    var('dummy')
#    out = integrate(func.subs(ul(0,0) == dummy),dummy)
#    return out.subs(dummy == ul(0,0))

### LIST OF LAGRANGIANS ###
	
if switch == 'lin':
    L = 1/2*(ul(0,1) - ul(1,0)) * (ul(1,1) - ul(0,0)) + 1/2*(a^2 - b^2) * (quot12)^2
    
elif switch == 'Q1': 
    if delta == 0:
        L = a^2*logseries(diff1/(-2*a*v_1), numvars-2) - b^2*logseries(diff2/(-2*b*v_1), numvars-2) - (a^2 - b^2)*logseries(quot12/(-2*v_1), 2*numvars-3)
    else:
        L1 = (-diff1+a^2*delta)*logseries( (-diff1+a^2*delta)/(2*a*v_1), numvars-1) - (-diff1-a^2*delta)*logseries( (-diff1-a^2*delta)/(2*a*v_1), numvars-1)
        L2 = L1.subs(a == b)
        L12 =(diff12+(a^2-b^2)*delta)*logseries( (quot12+(a+b)*delta)/(-2*v_1), 2*numvars-2) - (diff12-(a^2-b^2)*delta)*logseries( (quot12-(a+b)*delta)/(-2*v_1), 2*numvars-2)     
        L = L1 - L2 - L12
        
elif switch == 'Q2': 
    L1 = (sum1+a^2)*logseries((sum1+a^2)/(2*v_), numvars-1) - (sum1-a^2)*logseries((sum1-a^2)/(2*v_), numvars-1) + (-diff1+a^2)*logseries( (-diff1+a^2)/(2*a*v_1), numvars-1) - (-diff1-a^2)*logseries( (-diff1-a^2)/(2*a*v_1), numvars-1)
    L2 = L1.subs(a == b)
    L12 = (sum12+a^2-b^2)*logseries((sum12+(a^2-b^2))/(2*v_), 2*numvars-1) - (sum12-a^2+b^2)*logseries((sum12-(a^2-b^2))/(2*v_), 2*numvars-1) + (diff12+a^2-b^2)*logseries( (quot12+a+b)/(-2*v_1), 2*numvars-2) - (diff12-(a^2-b^2))*logseries( (quot12-a-b)/(-2*v_1), 2*numvars-2)     
    L = L1 - L2 - L12

elif switch == 'Q3': 
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

	L1 = a^2*diff1 - dilogseries(I*(diff1+a^2),numvars) + dilogseries(I*(diff1-a^2),numvars)
	if not(delta == 0):
		L1 += delta*( a^2*sum1 - dilogseries(I*(sum1+a^2),numvars) + dilogseries(I*(sum1-a^2),numvars) )
	L2 = L1.subs(a == b)
	L12 = (a^2-b^2)*diff12 - dilogseries(I*(diff12+a^2-b^2),2*numvars-1) + dilogseries(I*(diff12-a^2+b^2),2*numvars-1)
	if not(delta == 0):
		L12 += delta*( (a^2-b^2)*sum12 - dilogseries(I*(sum12+a^2-b^2),2*numvars-1) + dilogseries(I*(sum12-a^2+b^2),2*numvars-1) )
	L = (L1 - L2 - L12)/(2*I)

elif switch == 'Q4w':
	#L1 = -1/2* a^2*log(wsigma(sum1/2)) - a/2*(log(-diff1/2/a+a)-1)*(-diff1/4/a+a/2) + a/2*(log(-diff1/2/a-a)-1)*(-diff1/4/a-a/2) #- 1/6*wp(sum1/2)*(a^2-b^2)^3
	#L2 = L1.subs(a == b)
	#L12 =-1/2* (a^2-b^2)*log(wsigma(sum12/2)) - (a-b)/2*(log(-quot12/2+a+b)-1)*(-quot12/4+(a+b)/2) + (a-b)/2*(log(-quot12/2-a-b)-1)*(-quot12/4-(a+b)/2) #+ 1/2400*g2*(-diff12+a^2*b^2)^5 - 1/2400*g2*(-diff12-a^2*b^2)^5
	#L = L1 - L2 - L12
	
	L1 = a^2*logseries(wsigma(sum1)/wsigma(2*v_),numvars-2*0)
	#L1 += +a*( (log(diff1/a+a)-1)*(diff1/2/a+a/2) )
	#L1 += -a*( (log(diff1/a-a)-1)*(diff1/2/a-a/2) )
	L1 += +a*( (logseries(diff1/(-2*a*v_1)+a/(-2*v_1),numvars-1) + log(-2*v_1))*(diff1/2/a+a/2) )
	L1 += -a*( (logseries(diff1/(-2*a*v_1)-a/(-2*v_1),numvars-1) + log(-2*v_1))*(diff1/2/a-a/2) )
	L2 = L1.subs(a == b)
	
	L12 = ( (a^2-b^2)*logseries(wsigma(sum12)/wsigma(2*v_),2*numvars-3) - 1/6*wp(sum12)*(a^2-b^2)^3 )
	
	#L12 += (a-b) * (log(quot12+a+b)-1)*(quot12/2+(a+b)/2) 
	#L12 += -(a-b) * (log(quot12-a-b)-1)*(quot12/2-(a+b)/2) 
	L12 += (a-b) * (logseries(quot12/(-2*v_1)+(a+b)/(-2*v_1),2*numvars-2) + log(-2*v_1))*(quot12/2+(a+b)/2) 
	L12 += -(a-b) * (logseries(quot12/(-2*v_1)-(a+b)/(-2*v_1),2*numvars-2) + log(-2*v_1))*(quot12/2-(a+b)/2) 
	#L12 += killpowers(expand((a-b) * (- g2/1200*(-diff12 + a^2-b^2)^5 - 0*g3/5880*(-diff12 + a^2-b^2)^7 )))
	#L12 += killpowers(expand(-(a-b) * (- g2/1200*(-diff12 - a^2+b^2)^5 - 0*g3/5880*(-diff12 - a^2+b^2)^7 )))
	L12 += killpowers(expand( -1 * g2/2400*( 10*diff12^4*(a^2-b^2) + 20*diff12^2*(a^2-b^2)^3 + 2*(a^2-b^2)^5 ) ))
	L12 += killpowers(expand( -1 * g3/11760*( 14*diff12^6*(a^2-b^2) ) )) #+ 280*diff12^4*(a^2-b^2)^3) ))
	L = L1 - L2 - L12
	    
elif switch == 'H1': 
    L = 1/2*(-diff12) * (dertovar( (ul(1,1) - ul(0,0)).series(a,numvars+1).truncate().series(b,numvars).truncate() ) + 2/b + 2/a) + (1/a^2-1/b^2) * logseries(1-a*b*quot12, numvars, numvars+2)#.subs([a^(numvars+2+i)==0 for i in [1..numvars]]).subs([b^(numvars+2+i)==0 for i in [1..numvars]])

elif switch == 'H3tan' or switch == 'H3sin':
	quot12ab = killpowers(expand( (a+b)*quot12 )) #.series(a,numvars+1).truncate().series(b,numvars).truncate() )
	var('x','y','z')
	L1 = 1/8*diff1^2
	L2 = L1.subs(a == b)
	#integrand = arctan(expand(x*tan(y/4))).series(y,2*numvars-2).truncate()
	#L12 = integrate(integrand, y).subs(y == diff12).subs(x == (a+b)/(a-b))
	integrand = arctan(x * tan(y/4)).series(y,2*numvars-2).truncate()
	L12 = expand( integrate(integrand, y).subs(x == z/y) )
	ypowers = [1]
	zpowers = [1]
	for i in [1..2*numvars-2]:
		print i
		ypowers += [killpowers(expand(ypowers[i-1]*diff12))]
		zpowers += [killpowers(expand(zpowers[i-1]*quot12ab))]
	for i in reversed([1..2*numvars-2]):
		print i
		L12 = L12.subs(z^i == zpowers[i])
		L12 = L12.subs(y^i == ypowers[i])
	#L12 = expand( integrate(integrand, y).subs(x == z/y) ).subs(z == (a+b)*quot12).subs(y == diff12)
	
	#@parallel
	#def expand_and_kill(func):
	#	return killpowers(expand(killpowers(func)))	
	#Lgen = expand_and_kill([t for t in summands(L12)])
	#L = L1 - L2
	#for i in Lgen:
	#	L += i[1]
	L = expand(killpowers(L1 - L2 - L12)) #expand?

elif switch == 'BSQ':
	L = (1/a^3 - 1/b^3)*logseries(1 + a*b*quot12, numvars, numvars+3) - (1/a^2+1/a/b+1/b^2)*(ul(0,1)-ul(1,0)) + (1/a+1/b)*(ul(0,1)-ul(1,0))*ul(1,1) + (1/a-1/b+ul(0,1)-ul(1,0))*ul(0,0)*ul(1,1) - 1/a*ul(0,0)*ul(0,1) + 1/b*ul(0,0)*ul(1,0)
	L += 1/2*( ul(0,0)*ul(1,0)/a - ul(0,1)*ul(1,1)/a + ul(1,0)*ul(1,1)/b - ul(0,0)*ul(0,1)/b )
	L += 1/3*( ul(1,0)^3 - ul(0,1)^3 )
	L += 1/4*( (ul(1,0)^2-ul(0,0)^2)/a + (ul(1,1)^2-ul(1,0)^2)/b + (ul(0,1)^2-ul(1,1)^2)/a + (ul(0,0)^2-ul(0,1)^2)/b )
	L = dertovar(expand(L).series(a,numvars+1).truncate().series(b,numvars).truncate())

elif switch == 'BSQ3':
	allpde = []
	L = (1/a^3 - 1/b^3)*logseries(1 + a*b*quot12, numvars, numvars+3) - (1/a^2+1/a/b+1/b^2)*(ul(0,1,1)-ul(1,0,1)) + (1/a+1/b)*(ul(0,1,1)-ul(1,0,1))*ul(1,1,1) - (ul(0,1,1)-ul(1,0,1))*ul(1,1,2) - ul(0,0,1)*(ul(0,1,2) - ul(1,0,2) - (1/a - 1/b + ul(0,1,1) - ul(1,0,1))*ul(1,1,1) + 1/a*ul(0,1,1) - 1/b*ul(1,0,1))
	L += 1/2*( ul(0,0)*ul(1,0)/a - ul(0,1)*ul(1,1)/a + ul(1,0)*ul(1,1)/b - ul(0,0)*ul(0,1)/b )
	L += 1/3*( ul(1,0)^3 - ul(0,1)^3 )
	L += -( ul(1,0)*ul(1,0,2) - ul(0,1)*ul(0,1,2) )
	L += 1/4*( (ul(1,0)^2-ul(0,0)^2)/a + (ul(1,1)^2-ul(1,0)^2)/b + (ul(0,1)^2-ul(1,1)^2)/a + (ul(0,0)^2-ul(0,1)^2)/b )
	#L += 1/4*( (ul(1,0)*ul(1,0,2)-ul(0,0)*ul(0,0,2))/a + (ul(1,1)*ul(1,1,2)-ul(1,0)*ul(1,0,2))/b + (ul(0,1)*ul(0,1,2)-ul(1,1)*ul(1,1,2))/a + (ul(0,0)*ul(0,0,2)-ul(0,1)*ul(0,1,2))/b )
	L = dertovar(expand(L).series(a,numvars+1).truncate().series(b,numvars).truncate())

elif switch == 'GD4':
	var('k,l,m')
	L = -(1/a^4 - 1/b^4)*logseries(1 + a*b*quot12, numvars, numvars+3) + (1/a^3+1/a^2/b+1/a/b^2+1/b^3)*(ul(0,1)-ul(1,0)) 
	L += - (1/a^2+1/a/b+1/b^2)*(ul(0,1)-ul(1,0))*ul(1,1) + (1/a+1/b)*(ul(0,1)-ul(1,0))*ul(1,1,2) 
	L += ((1/a+1/b)*ul(0,0)-ul(0,0,3))*(ul(0,1,2) - ul(1,0,2) - (1/a-1/b+ul(0,1)-ul(1,0))*ul(1,1) + 1/a*ul(0,1) - 1/b*ul(1,0))
	L += (1/a-1/b+ul(0,1)-ul(1,0))*ul(0,0)*ul(1,1,2) - 1/a*ul(0,0)*ul(0,1,2) + 1/b*ul(0,0)*ul(1,0,2)
	
	L += -1/2*( ul(0,0)*ul(1,0)/a^2 - ul(0,1)*ul(1,1)/a^2 + ul(1,0)*ul(1,1)/b^2 - ul(0,0)*ul(0,1)/b^2 )
	L += -1/3*( ul(0,0)*ul(1,0)^2/a - ul(0,1)*ul(1,1)^2/a + ul(1,0)*ul(1,1)^2/b - ul(0,0)*ul(0,1)^2/b )
	#L += -1/6*( ul(0,0)^2*ul(1,0)/a - ul(0,1)^2*ul(1,1)/a + ul(1,0)^2*ul(1,1)/b - ul(0,0)^2*ul(0,1)/b )
	L +=      ( ul(0,0)*ul(1,0,2)/a - ul(0,1)*ul(1,1,2)/a + ul(1,0)*ul(1,1,2)/b - ul(0,0)*ul(0,1,2)/b )
	L += -1/4*( (ul(1,0)^2-ul(0,0)^2)/a^2 + (ul(1,1)^2-ul(1,0)^2)/b^2 + (ul(0,1)^2-ul(1,1)^2)/a^2 + (ul(0,0)^2-ul(0,1)^2)/b^2 )
	#L += k/3*( ul(1,0)^3 - ul(0,1)^3)
	#L += -k/4*( (ul(1,0)*ul(1,0,2)-ul(0,0)*ul(0,0,2))/a^2 + (ul(1,1)*ul(1,1,2)-ul(1,0)*ul(1,0,2))/b^2 + (ul(0,1)*ul(0,1,2)-ul(1,1)*ul(1,1,2))/a^2 + (ul(0,0)*ul(0,0,2)-ul(0,1)*ul(0,1,2))/b^2 )
	#L += 1/4*( (ul(1,0)^2-ul(0,0)^2)/a + (ul(1,1)^2-ul(1,0)^2)/b + (ul(0,1)^2-ul(1,1)^2)/a + (ul(0,0)^2-ul(0,1)^2)/b )
	L = L*a*b
	L = dertovar(expand(L).series(a,numvars+1).truncate().series(b,numvars).truncate())
	
elif switch == 'GD4skew':
	var('k,l,m')
	L = -(1/a^4 - 1/b^4)*logseries(1 + a*b*quot12, numvars, numvars+3) + (1/a^3+1/a^2/b+1/a/b^2+1/b^3)*(ul(0,1)-ul(1,0)) 
	L += - (1/a^2+1/a/b+1/b^2)*(ul(0,1)-ul(1,0))*ul(1,1) + (1/a+1/b)*(ul(0,1)-ul(1,0))*(ul(1,1,2)+ul(1,1,3))/2
	L += ((1/a+1/b)*ul(0,0)-ul(0,0,2)/2+ul(0,0,3)/2)*(ul(0,1,2)/2 - ul(1,0,2)/2 + ul(0,1,3)/2 - ul(1,0,3)/2 - (1/a-1/b+ul(0,1)-ul(1,0))*ul(1,1) + 1/a*ul(0,1) - 1/b*ul(1,0))
	L += (1/a-1/b+ul(0,1)-ul(1,0))*ul(0,0)*(ul(1,1,2)/2+ul(1,1,3)/2) - 1/a*ul(0,0)*ul(0,1,2)/2 + 1/b*ul(0,0)*ul(1,0,2)/2 - 1/a*ul(0,0)*ul(0,1,3)/2 + 1/b*ul(0,0)*ul(1,0,3)/2
	
	L += -1/2*( ul(0,0)*ul(1,0)/a^2 - ul(0,1)*ul(1,1)/a^2 + ul(1,0)*ul(1,1)/b^2 - ul(0,0)*ul(0,1)/b^2 )
	L += -1/4*( (ul(1,0)^2-ul(0,0)^2)/a^2 + (ul(1,1)^2-ul(1,0)^2)/b^2 + (ul(0,1)^2-ul(1,1)^2)/a^2 + (ul(0,0)^2-ul(0,1)^2)/b^2 )
	
	L += -1/6*( ul(0,0)*ul(1,0)^2/a - ul(0,1)*ul(1,1)^2/a + ul(1,0)*ul(1,1)^2/b - ul(0,0)*ul(0,1)^2/b )
	L += -1/6*( ul(0,0)^2*ul(1,0)/a - ul(0,1)^2*ul(1,1)/a + ul(1,0)^2*ul(1,1)/b - ul(0,0)^2*ul(0,1)/b )
	L +=  1/2*( ul(0,0)*ul(1,0,2)/a - ul(0,1)*ul(1,1,2)/a + ul(1,0)*ul(1,1,2)/b - ul(0,0)*ul(0,1,2)/b )
	L +=  1/2*( ul(0,0,2)*ul(1,0)/a - ul(0,1,2)*ul(1,1)/a + ul(1,0,2)*ul(1,1)/b - ul(0,0,2)*ul(0,1)/b )
	L +=  1/2*( (-ul(0,0)*ul(0,0,2))/a + (-ul(1,0)*ul(1,0,2))/b + (ul(0,1)*ul(0,1,2))/a + (ul(0,0)*ul(0,0,2))/b )
	L +=  1/2*( ul(0,0)*ul(1,0,3)/a - ul(0,1)*ul(1,1,3)/a + ul(1,0)*ul(1,1,3)/b - ul(0,0)*ul(0,1,3)/b )
	
	L += +1/8*( ul(0,0,3)*ul(1,0,3) - ul(0,1,3)*ul(1,1,3) + ul(1,0,3)*ul(1,1,3) - ul(0,0,3)*ul(0,1,3) )
	L += -1/8*( ul(0,0,2)*ul(1,0,2) - ul(0,1,2)*ul(1,1,2) + ul(1,0,2)*ul(1,1,2) - ul(0,0,2)*ul(0,1,2) )
	L += -1/8*( ul(0,0,2)*ul(1,0,3) - ul(0,1,2)*ul(1,1,3) + ul(1,0,2)*ul(1,1,3) - ul(0,0,2)*ul(0,1,3) )
	L += -1/8*( ul(0,0,3)*ul(1,0,2) - ul(0,1,3)*ul(1,1,2) + ul(1,0,3)*ul(1,1,2) - ul(0,0,3)*ul(0,1,2) )
	L +=  1/4*( ul(0,0,2)*ul(1,0)^2 - ul(0,1,2)*ul(1,1)^2 + ul(1,0,2)*ul(1,1)^2 - ul(0,0,2)*ul(0,1)^2 )
	L +=  1/4*( ul(0,0)^2*ul(1,0,2) - ul(0,1)^2*ul(1,1,2) + ul(1,0)^2*ul(1,1,2) - ul(0,0)^2*ul(0,1,2) )
	
	L = dertovar(expand(L).series(a,lagnumvars+1).truncate().series(b,lagnumvars).truncate())
	
else:
	L = 0

### END OF LIST ###

### First series expansion
print walltime(w)
series = L.series(b,lagnumvars).truncate()
print walltime(w)

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
    return [0 for j in [1..i]] + [expand(dertovar( series2.coefficient(a,j) )) for j in [i+1..lagnumvars]]

cgen = expanda([i for i in [1..lagnumvars-1]])

Lcoeff = [[0 for j in [1..lagnumvars]] for i in [1..lagnumvars]]
for i in cgen:
    Lcoeff[i[0][0][0]-1] = i[1]

textadd("L disc = ",viewldisc)
latexadd(matrix(Lcoeff),viewldisc)

Lcoeff = matrix(Lcoeff) - matrix(Lcoeff).transpose()

print walltime(w)

### Euler-Maclaurin correction ###

#horizontal EM operator
@parallel
def diffEMcell(row,col,array):
	out = 0*a
	for i in [1..row]:
		if includeeven or is_odd(i):
			out += (-1)^(i+1) * miwaconst/i * vdiff(array[row-i][col], eval('t'+str(i)))
	return out

def EMdiff(array):
	w = walltime()
	out = [[0 for i in [1..lagnumvars]] for j in [1..lagnumvars]]
	cellgen = diffEMcell(flatten([[(i,j,array) for i in [0..lagnumvars-1]] for j in [0..lagnumvars-1]],max_level=1))
	for i in cellgen:
		out[i[0][0][0]][i[0][0][1]] = i[1]
	return matrix(out)

#powers of horizontal EM operator
EMdiffs = [matrix(lagnumvars,lagnumvars) for k in [0..lagnumvars]]
EMdiffs[0] = Lcoeff
for i in [1..lagnumvars]:
    EMdiffs[i] = EMdiff(EMdiffs[i-1])

#horizontal correction
lagarray = matrix(EMdiffs[0]) - 1/2*EMdiffs[1] + sum(bernoulli(2*i)/factorial(2*i)*EMdiffs[2*i] for i in [1..lagnumvars/2])

#powers of vertical EM operator
EMdiffs = [matrix(lagnumvars,lagnumvars) for k in [0..lagnumvars]]
EMdiffs[0] = lagarray
for i in [1..lagnumvars]:
    EMdiffs[i] = EMdiff(EMdiffs[i-1].transpose()).transpose()

#vertical correction
lagarray =  matrix(EMdiffs[0]) - 1/2*EMdiffs[1] + sum(bernoulli(2*i)/factorial(2*i)*EMdiffs[2*i] for i in [1..lagnumvars/2])
lagarray = scaling * lagarray * scaling

triang = utriang(lagarray)

### Ouptut
textadd("L pluri = ",viewlpluri)
latexadd(triang,viewlpluri)
