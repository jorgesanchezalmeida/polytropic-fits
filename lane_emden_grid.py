"""

Computes the grid of polytropes projected in the plane of the sky, to be
used in a least-squared fit.

Sanchez Almeida, Trujilo, and Plastino, 2021, ApJ, in press 

"""
#
# - imports
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import sys
import time
#from scipy.optimize import least_squares
from scipy.integrate import simps
import pickle # some sort of save restore mechanism for variables

#
print('-----> Running:',sys.argv[0],'parameters to be introduced by hand')
#
# - ODE to be solved
#    index is the index of the politrope. Not to be confused with the Sesic index
#
def deriv(yy,ss,index): # differential equation to be solved
    phi, zeta = yy
    derphi = zeta/ss**2.
    if phi >= 0:
        derzeta = -3.*phi**index*ss**2
    else:
        derzeta = 0.
    return derphi, derzeta
#
def lane_emden(index,ss): # routine that provides the solutions
    #
    # index and ss are imput. rho is density. le_half is the 1/2 radius. nss is the pixel of the 1st zero.  
    #
    #     intitial conditions
    ss0 = ss[0]
    phi0=1.0-ss[0]**2./2.# somewere in my notes. This is the approx.
    zeta0=0.0
    yy0 = phi0, zeta0
    nss = ss.size
    #
    #     actual solution
    atol = 1e-13#1e-13
    rtol = 1e-13#1e-13 get problem if smaller.
    sol, out = odeint(deriv, yy0, ss, args=(index,), atol=atol, rtol=rtol, full_output=True) # method? Uses something called Isoda. 
    print('index=',index,out['message'])
    #
    phi = sol[:,0]
    #    max size limit
    ii=0
    while phi[ii] > 0 and ii < nss-1:
        ii=ii+1
    if ii == nss-1:
        le_size=ss[ii]
    else:
        le_size=(ss[ii]*phi[ii-1]-ss[ii-1]*phi[ii])/(phi[ii-1]-phi[ii])
    nss=ii
    phi[nss:] =0.
    rho = phi**index
    #
    #   half light radius
    yy = rho*ss**2
    syy = np.cumsum(yy)
    syy = (0.5-syy/syy.max())
    ii=0
    while syy[ii] > 0:
        ii=ii+1
    le_half=(ss[ii]*syy[ii-1]-ss[ii-1]*syy[ii])/(syy[ii-1]-syy[ii])
    #
    return phi,rho,nss,le_half,le_size
#
# - do loop in the polytropic index
#
# polytropic index for the grid
#             lnindex = np.arange(1.,5.1,2)
gridname='grid0'
lnindex = np.concatenate((np.arange(1.,1.9999,0.1),np.arange(2.,5.01,0.01),np.arange(5.1,10.1,0.1)))
gamma = 0
#
gridname='grid1'
lnindex = np.concatenate((np.arange(1.,1.9999,0.1),np.arange(2.,5.01,0.01)))
gamma = 0.
# 
gridname='grid2'
lnindex = np.concatenate((np.arange(1.,2.9999,0.1),np.arange(3.,6.01,0.01),np.arange(6.1,10.1,0.1)))
gamma=0
# 
gridname='grid2g3'
lnindex = np.concatenate((np.arange(1.,2.9999,0.1),np.arange(3.,6.01,0.01),np.arange(6.1,10.1,0.1)))
gamma=0.3
# 
gridname='grid2g1'
lnindex = np.concatenate((np.arange(1.,2.9999,0.1),np.arange(3.,6.01,0.01),np.arange(6.1,10.1,0.1)))
gamma=0.1
#
gridname='grid3'
lnindex = np.concatenate((np.arange(1.,7.,0.1),np.arange(7.,100,1)))
gamma = 0.0
#
gridname='grid3g3'
lnindex = np.concatenate((np.arange(1.,7.,0.1),np.arange(7.,100,1)))
gamma = 0.3
#
#
gridname='grid3g1'
lnindex = np.concatenate((np.arange(1.,7.,0.1),np.arange(7.,100,1)))
gamma = 0.1
#
#gridname='grid4'
#lnindex = np.concatenate((np.arange(1.,7.,0.1),np.arange(7.,100,1)))
#gamma = 0.0
#
gridname='grid5'
lnindex = np.concatenate((np.arange(1.,7.,0.1),np.arange(7.,100,1),np.arange(100,1000,10)))
gamma=0
#
gridname='grid5g3'
lnindex = np.concatenate((np.arange(1.,7.,0.1),np.arange(7.,100,1),np.arange(100,1000,10)))
gamma=0.3

gridname='grid6'
lnindex = np.concatenate((np.arange(1.,20.,0.1),np.arange(20.,100,1)))
#
gridname='grid8'
lnindex = np.arange(1,10.1,0.1)
gamma=0
gridname='grid9'
lnindex = np.arange(1,10.1,0.1)
gamma=0.3
gridname='gridtest'
lnindex = np.arange(1.,5.1,0.1)
gamma=0.0
#
print('Evaluating '+gridname)
print('Pol Index between ',lnindex.min(),' and ',lnindex.max())
#
nlnindex = len(lnindex)
#
lss = np.arange(-5,3.,0.00001) # grid for the 3D
ss = 10**lss 
nss = len(ss)
#
lrrmin=np.log10(0.01)
lrrmax=np.log10(200)
lrr = np.arange(lrrmin,lrrmax,0.02)# RR to project in the plane of sky. Abel trabsform.
rr = 10**lrr
nrr = len(rr)
#
pp = np.zeros((nrr,nlnindex),dtype='float')
rrho = np.zeros((nss,nlnindex),dtype='float')
#
start_time = time.time()
kk = -1
for index in lnindex:
    kk = kk+1
    phi,rho,nss,le_half,le_size = lane_emden(index,ss) # integrating the Lane_Emden Eq.
    #
    #print('central density=',rho.max())
    here = -1
    ii = rr*0.0
    for radius in rr:  # Abel transform
        here = here+1
        """
        # several other integration schemes that provide similar results
        #
        #select = (ss > radius).nonzero()[0] original from lane_emden5.py
        select = (np.sqrt((ss/radius)**2-1) > 0.01).nonzero()[0]
        #
        xx = ss[select]
        yy = rho[select]
        #
        integrand = yy*xx/np.sqrt(xx**2-radius**2)
        ii[here] = simps(integrand,xx,even='first') # first: Use Simpson’s rule for the first N-2 intervals with a trapezoidal rule on the last interval
        #ii[here] = np.trapz(integrand,xx) #
        
        #
        deltar = ss[select[0]]-radius
        #extra = radius*rho[select[0]-1]*(deltar/radius+np.sqrt(2*deltar/radius))
        #extra = 1+deltar/radius
        #extra = radius*rho[select[0]-1]*np.log(extra+np.sqrt(extra**2-1))
        extra = 1+deltar/radius
        extra = radius*rho[select[0]-1]*np.sqrt(extra**2-1)
        ii[here] = ii[here] + extra
        #
        """
        select = (ss >= radius).nonzero()[0] #original from lane_emden5.py
        xx = ss[select]
        yy = rho[select]**(1+gamma)
        uu = np.sqrt(xx**2-radius**2)
        ii[here] = simps(yy,uu,even='first')
        extra = (rho[select[0]-1]**(1+gamma)+rho[select[0]]**(1+gamma))/2 * np.sqrt(xx[0]**2-radius**2)
        ii[here] = ii[here] + extra
    #
    # to be saved
    pp[:,kk] = ii*2. #the factor of 2 is in the definition of Abel Transform. 
    rrho[:,kk] = rho # for testing
    print(index,le_half)
    #
end_time = time.time()
#
# -general info.
#
print('----')
print('required time [s] =',end_time - start_time)
print('time/profile [s] =',(end_time - start_time)/lnindex.size)
print('number of radial positions=',rr.size)
print('number of poly index=',lnindex.size)
#
# --------------------------------------------
# save the fits to be used elsewhere
#
# save restore mechanism
file_save = './lane_emden_grid/'+gridname+'.p'
#
data ={'rr':rr,
    'ss':ss, # to see the conditions of integration; ss.min() ss.max() etc
    'pp':pp, # this contains the actual shape in 2D. pp[rr,lnindex]
    'rrho':rrho, # too large
    'lnindex':lnindex,
    'gamma':gamma}
#
pickle.dump(data, open(file_save, "wb"))
print('grid data saved at ',file_save)
print(' ')



#----------------------------------
#
# testing the routines - rho
#
#
neq5 = np.argmin(np.abs(lnindex-5))
sapprox5 = pp[:,neq5]
approx5 = rrho[:,neq5]
analytical5 = 1/(1+ss**2)**(5./2.)
sanalytical5 = (4./3.)/(1+rr**2)**2.
residual5 = np.abs(approx5/analytical5-1.)
sresidual5 = np.abs(sapprox5/sanalytical5-1.)
numberss = ss.size
residual5[numberss-1] = residual5[numberss-2] # apporox is zero in this last pixels so problems for relative error
#
neq1 = np.argmin(np.abs(lnindex-1))
approx1 = rrho[:,neq1]
f3 = np.sqrt(3.)
analytical1 = np.sin(f3*ss)/(f3*ss)
select = (ss < np.pi/f3).nonzero()[0]
select = (approx1 > 0).nonzero()[0]
residual1 = 0.* approx1
residual1[select] = np.abs(approx1[select]/analytical1[select]-1.)
#
#
fs = 15
plt.rcParams.update({'font.size': fs})
plt.ion()
plt.close(1) # to avoid getting multiple windows
fig, ax = plt.subplots(1,1)#, sharey = False, sharex = False)
ax.set_yscale("log")
ax.set_xlabel(r"$s$, $x$")
ax.set_ylabel(r"$\psi(s)$, $f(x,m)$")
ax.set_ylim(1e-6,2.)
ax.set_xlim(1e-3,100.)
ax.set_xscale("log")
p1 = ax.plot(ss,approx1,label=r'$\psi(s)$, m='+'{:.1f}'.format(lnindex[neq1]),lw=2)
color1 = p1[0].get_color()
ax.plot(ss,residual1*1e6,label=r'Relative Error $\times 10^{6}$',ls='--',c=color1,lw=2)
p5 = ax.plot(ss,approx5,label=r'$\psi(s)$, m='+'{:.1f}'.format(lnindex[neq5]),lw=2)
color5 = p5[0].get_color()
ax.plot(ss,residual5*1e6,label=r'Relative Error $\times 10^{6}$',ls='--',c=color5,lw=2)
p5 = ax.plot(rr,sapprox5,label=r'$f(x,m)$, m='+'{:.1f}'.format(lnindex[neq5]),lw=2)
color5 = p5[0].get_color()
ax.plot(rr,sresidual5*1e6,label=r'Relative Error $\times 10^6$',ls='--',c=color5,lw=2)
ax.legend(loc=2,fontsize=fs*0.6)
outfile = 'lane_emden_grid1_test1.pdf'
print('output at '+outfile)
plt.savefig(outfile,transparent='False',bbox_inches='tight') 
plt.show()

