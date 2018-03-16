#!/usr/bin/env sage

### Initialize latex output array
output = []

### Open output files - data in filename
if not(switch==None):
	filename = switch
	if not(delta == 0):
		filename += '-delta'
	filename += '-' + str(numvars)
	if includeeven:
		filename += 'full'
	plaindoc = open('datadump/' + filename + '-plain','w')
	latexdoc = open('datadump/' + filename + '-tex','w')

### print to latex and plaintext files
### if all==True, print to console and pdf as well
def latexadd(eqn,all=True):
	global output
	if all:
		print latex(eqn)
		output += ['\\tiny ' + latex(eqn)]
	if not(latexdoc.closed):
		latexdoc.write(latex(eqn) + '\n\n')
	if not(plaindoc.closed):
		plaindoc.write(str(eqn) + '\n\n')
def textadd(str,all=True):
	global output
	if all:
		print str
		output += [str]
	if not(latexdoc.closed):
		latexdoc.write(str + '\n\n')
	if not(plaindoc.closed):
		plaindoc.write(str + '\n\n')

