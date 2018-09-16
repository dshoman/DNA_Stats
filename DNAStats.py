import numpy as np
from scipy.special import comb, factorial
import matplotlib.pyplot as plt
import json

import time, datetime

class DNA_Stats:
	'''This class includes functions to calculate the probability of detecting an arbitrary number of DNA variants in a sample of arbitrary size. The total number of variants is V, the sample size is N. It will also include a method to create a figure comparing the combinatoric solution with an approximation based on the Poisson distrubution.'''

	def __init__(self,V=0,N=0):
		self.V = V
		self.N = N
		if V==0:
			try: 
				self.V = int(raw_input('What is V?\n'))
			except(ValueError):
				self.V = 2
				print 'Error in Input: V set to 2'
		if N==0:
			try:
				self.N = int(raw_input('What is N?\n'))
			except(ValueError):
				self.N = 3
				print 'Error in Input: N set to 3'
		


	@staticmethod
	def break_integer(I,M,n,printing='no'):
		'''Divide an integer (I) into n integer parts, taking out the part (M) as the largest of the fractions. All other fractions should be =< M. If the method is called outright, and the results should be printed to screen, set kwarg 'printing' to 'yes'. The method will check the directory 'precalc_split' to see if the calculation has been done before, and if not add the results to a file there'''

#		t_start = time.clock()

		try:
			with open('precalc_split/V_'+str(n)+'_N_'+str(I),'r') as f_in:
				for row in f_in:
					r = json.loads(row)
					if r[0][0] == n and r[0][1] == I:
						return r[1:]
		except(IOError): pass

		a, comp, comp_prev  = [],n*[0],n*[0]
		r,m = I,M
		count, nrep = 1,0
		while count<=n:
			for ii in np.flip(range(m))+1:
#				print nrep, count,r,m,ii,comp,a, 'V_'+str(n-count+1)+'_N_'+str(r)
				if r>4 and n-count>4:
					try:
						with open('precalc_split/V_'+str(n-count+1)+'_N_'+str(r),'r') as f_in:
#							print 'hello there'
							for row in f_in:
								row_ = json.loads(row)[1:]
								for el in row_:
#									print comp[count-2], el, el[0],type(el[0])
									if el[0]<=comp[count-2]:
										a += [comp[:count-1] + el]
							if not row_[-1][0]>comp[count-2]:
								count+=-1
								nrep = 0
								m = a[-1][count-1]-1
								r += a[-1][count-1]
								break
					except(IOError):
						print 'no'
						pass
				if r-ii >= n-count:
					r = r-ii
					comp[count-1] = ii
					count+=1
					m = ii
					break
			else:
#				print 'hello'
				if count == 1:
					break
				count += -1
				r+=comp[count-1]
				m = a[-1][count-1]-1
				nrep += 1
#				print nrep, count, r, m, np.flip(range(m))+1
			if count==n+1:
				if sum(comp)==I and comp!=comp_prev:# not in a:
					a += [comp[:]]
					nrep = 0
					comp_prev = comp[:]
				count = n-1-nrep
				r+=sum(comp[count-1:])
				m = comp[count-1]-1

		with open('precalc_split/V_'+str(n)+'_N_'+str(I),'a') as f_out:
			f_out.write(str( [[n,I,len(a)]] + a )+'\n')
		if printing=='yes':
			print a
		return a


	def find_probability(self,printing='no'):
		'''Find the probability of finding V variants in sample size N (where V and N have been predefined in the class instance)'''
		if self.V>self.N: return 0
		try:
			with open('precalc_prob/V_'+str(self.V), 'r') as f_in:
				for row in f_in:
					r = json.loads(row)
					if int(r[0])==self.N:
						if printing=='yes': print 'The probability of finding all V={} in sample size N={} is: {}'.format(self.V,self.N,r[1])
						return r[1]
		except(IOError): pass

		broken_int = self.break_integer(self.N,self.N-self.V+1,self.V)
		N_poss = 0
		for bi in broken_int:
			take, occurence = np.unique(bi, return_counts=True)
			f = np.prod(factorial(occurence))
			enum, rest = 1, self.N
			for i in bi:
				enum = enum*comb(rest,i)
				rest = rest-i
			N_poss += enum/f
		res = N_poss*factorial(self.V)/float(self.V)**float(self.N)

		with open('precalc_prob/V_'+str(self.V),'a') as f_out:
			f_out.write( '[{},{}]\n'.format(self.N,res) )
		if printing=='yes':
			print 'The probability of finding all V={} in sample size N={} is: {}'.format(self.V,self.N,res)
		return res


	def mc_approx(self,V=None,N=None,N_mc=100000):
		'''Monte Carlo simulation of the selection process.'''
		if not V: V = self.V
		if not N: N = self.N
		if V>N: return 0
		draw = []
		for ii in range(N_mc):
			draw += [len(np.unique(np.random.randint(V,size=N)))]
		return float(len([ii for ii in draw if ii==V]))/N_mc

	
	def mc_alpha(self,V=None,N_start=None,N_mc=100000,alpha=.99):
		'''Find the minimum sample size for a given V, to reach the probability threshold alpha.'''
		if not V: V = self.V
		if not N_start: N_start = self.N
		p = self.mc_approx(V=V,N=N_start,N_mc=N_mc)
		if p<alpha:
			while p<alpha:
				N_start+=1
				p = self.mc_approx(V=V,N=N_start,N_mc=N_mc)
			else: print 'Required sample size for V={} and alpha={} :'.format(V,alpha,N_start)
		elif p>alpha:
			while p>alpha:
				N_start+=-1
				p = self.mc_approx(V=V,N=N_start,N_mc=N_mc)
			else: print 'Required sample size for V={} and alpha={} : {}'.format(V,alpha,N_start+1)


	def make_figure(self,V=5,N_max=10,exact='yes',poisson='no',mc='no',N_mc=100000, name=None):
		'''Creates a plot of the probability for a given number of variables, and up to a given sample size N_max. The horizontal lines are .95 and .99 confidence levels repsectively.\n
		Arguments:\n
		V: the number of variants\n
		N_max: the largest sample size to be included in the plots\n
		exact: include the exact solution in the plot (not advisable for large V or N_max)\n
		poisson: include a simple approximation, based on the assumption that finding each variant is an independent Poisson process\n
		mc: include an mc simulation of the experiment (use for large V or N_max)\n
		N_mc: the number of iterations for the MC simulation\n
		name: name for the figure, if not used, the default name will be used
		'''
		self.V = V
		xrng = np.arange(0,N_max,step=1)
		ymin,ymax = 0,1.05
		fig,ax =  plt.subplots(1,figsize=(6,4))
		ax.axis([xrng[0],xrng[-1],ymin,ymax])

		if exact=='yes':
			y_prob_comb = []
			for x in xrng:
				self.N = x
				y_prob_comb += [self.find_probability()]
			ax.plot(xrng,y_prob_comb,c='r',label='Exact Solution')
		if poisson=='yes':
			ax.plot(xrng,np.power(1-np.exp(-xrng.astype(float)/V),V),c='b',label='(1-exp(-N/{}))^{}'.format(str(V),str(V)))
		if mc=='yes':
			y_prob_mc = []
			for x in xrng:
				self.N =  x
				y_prob_mc += [self.mc_approx(N_mc=N_mc)]
			ax.plot(xrng,y_prob_mc,c='k',ls=':',label='MC')

		ax.hlines(0.99,xrng[0],xrng[-1],linestyles='-.',colors='grey')
		ax.hlines(0.95,xrng[0],xrng[-1],linestyles='--',colors='grey')
		ax.set_xlabel('N')
		ax.set_ylabel( 'Prob(all {} in N)'.format(str(V)) )
		ax.legend(loc='lower right')
		if not name: name = 'P_Finding_{}_in_{}.pdf'.format(V,N_max)
		plt.savefig('figures/'+name)



##Test runs:
#t_begin = time.clock()

#DNA_Stats().find_probability()
#DNA_Stats.break_integer(22,16,7,printing='yes')
#DNA_Stats.break_integer(23,17,7,printing='yes')
#DNA_Stats.break_integer(42,35,8,printing='yes')
#DNA_Stats(V=5,N=10).find_probability(printing='yes')
#DNA_Stats(V=9,N=60).find_probability_old(printing='yes')
#DNA_Stats(V=50,N=384).find_probability(printing='yes')
#t_now = time.clock()-t_start; print datetime.timedelta(seconds=t_now)
#DNA_Stats(V=50,N=500).find_probability(printing='yes')
#t_now = time.clock()-t_start; print datetime.timedelta(seconds=t_now)
#DNA_Stats(V=5,N=10).mc_approx()

#DNA_Stats(V=1,N=1).make_figure(V=10,N=100)
#DNA_Stats(V=1,N=1).make_figure(V=10,N_max=100,mc='yes')

#DNA_Stats(V=1,N=1).mc_alpha(V=10,N_start=69)

#t_now = time.clock()-t_begin; print 'Finished in: ', datetime.timedelta(seconds=t_now)

#for ii in range(12)[2:]:
#	for jj in range(80)[1:]:
#		DNA_Stats(V=ii,N=jj).find_probability()
#		t_now = time.clock()-t_begin; print 'Finished V={} N={} in: '.format(ii,jj), datetime.timedelta(seconds=t_now)

#t_now = time.clock()-t_begin; print 'Finished in: ', datetime.timedelta(seconds=t_now)



