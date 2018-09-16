from DNAStats import DNA_Stats
import time, datetime

t_begin = time.clock()

#Basic functionality; requires user input; print results to screen
DNA_Stats().find_probability(printing='yes')

#Create a plot of the probability of finding all 5 varieties (V) is a sample size from 5 to 40 (N_max)
#This plot will include an exact calculation and the simplified Poisson estimate
DNA_Stats(V=1,N=1).make_figure(V=5,N_max=30,exact='yes',poisson='yes',mc='no')

#Create a plot of the probability of finding all 50 varieties (V) is a sample size from 50 to 500 (N_max)
#This plot will include an MC simulation and the simplified Poisson estimate
# (NOTE: for this size V, the exact calculation is likely too time-consuming to perform)
DNA_Stats(V=1,N=1).make_figure(V=50,N_max=500,exact='no',poisson='yes',mc='yes',N_mc=10000)

#Calculate the required sample size for a 95% confidence level to detect all 5 varieties (V) in sample size 25 (N).
#This method always uses an MC simulation.
DNA_Stats(V=5,N=25).mc_alpha(alpha=.95,N_mc=10000)

t_now = time.clock()-t_begin; print 'Finished in: ', datetime.timedelta(seconds=t_now)
