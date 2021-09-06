"""
procedure to fit surface density profiles using the density associated 
with polytropes

Sanchez Almeida, Trujilo, and Plastino, 2021, ApJ, in press 

"""

# - imports
#from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import sys
import time
from scipy.optimize import least_squares
from scipy.integrate import simps
from scipy.interpolate import interp2d #griddata
import pickle # some sort of save restore mechanism for variables

#
print('-----> Running:',sys.argv[0],' galaxno bbini mmini laaini rr_limit  sb_limit')
#
#
# - Routine to compute the projected politropes
#
#         reading the grid from 'lane_emden_grid1.py'
#
#gridname = 'grid0'
#gridname = 'grid1'
gridname = 'grid2'
#gridname='grid2g3'
#gridname = 'grid3'
file_save = './lane_emden_grid/'+gridname+'.p'
data = pickle.load(open(file_save,"rb"))
#
rr = data['rr'] # projected radial distance , equispaced in log
ss = data['ss'] # radial distance in 3D used to compute the grid.
pp = data['pp'] # shapes.
lnindex = data['lnindex']
print('Using '+gridname)
print('Index from ',lnindex.min(),' to ',lnindex.max())
#
lrr = np.log10(rr)
pp[np.where(pp == 0)]=1.e-50
lpp = np.log10(pp)
#
#        defines the interpolarion
polyshape = interp2d(lnindex,lrr, lpp, kind='linear', copy=False, bounds_error=False, fill_value=-1)
#
#        function to be called during LS minimization.
def lane_emden_suite(coeff,xdata,ydata,polyshape):
    bb = coeff[0]
    mm = coeff[1]
    lxx = xdata - np.log10(bb) # x axis is in
    lppinter = polyshape([mm],lxx)
    lppinter = lppinter[:,0]
    #
    laa = ydata.mean() - lppinter.mean() # it is a constant in the lograithm amplitud by least squares
    coeff[2] = laa
    #
    residual = lppinter+laa - ydata
    #print(np.shape(lppinter),np.shape(lxx),np.shape(mm))
    return residual
#
# - reading the data to be fitted
#
filedata ='./data/galaxies_edges_wnew_mass.profiles'
profile = pickle.load(open(filedata, 'rb')) 
profno = 20
profno = 450
profno = 2
profno = 123
profno = np.int(sys.argv[1]) if len(sys.argv) > 1 else 2
scale = profile['scale'][profno]/0.396 # Nacho says 'scale' is kpc/pixel. 0.396 is arcsec/pixel in SDSS images.
                                    # seems to be ok since re-scaled 'scale' agrees with the range provided the range of redshifts. 
                                    # Importantly, seeing FWHM is 1.4 in 'r', so 2 arcsec is a lot!!! (See DR7 webpage)
mstar = profile['mass'][profno]
dd = profile['step'][profno]
ldd = np.log10(dd)
arcsec = dd/scale
xdata = np.log10(dd)
#         to know what is in: profile.keys()    
mu_g_obs = profile['mu_g_obs'][profno]
rho = profile['rho'][profno]
bb = np.float(sys.argv[2]) if len(sys.argv) > 2 else 5.
mm = np.float(sys.argv[3]) if len(sys.argv) > 3 else 4.
laa = np.float(sys.argv[4]) if len(sys.argv) > 4 else 1.3
rr_limit = np.float(sys.argv[5]) if len(sys.argv) > 5 else 0.7
sb_limit = np.float(sys.argv[6]) if len(sys.argv) > 6 else 29
select = np.where(mu_g_obs < sb_limit) # ---- bright enough
xdata = xdata[select]
ydata = rho[select]
arcsec = arcsec[select]
select = np.where(arcsec > rr_limit) # ---- taking out the PSF ... very conservative. 
xdata = xdata[select]
ydata = ydata[select]
#
# - Least squares fit to PP (projected politropes)
#
coeff0 = np.array([bb,mm,laa])
bounds = np.transpose([[0.1,50.],[1.5,5.05],[-1.,np.inf]])
result = least_squares(lane_emden_suite, coeff0, args=(xdata,ydata,polyshape), method='lm')#,max_nfev = 1000)
coeff = result.x
bb = coeff[0]
laa = coeff[2]
#print('laa=',laa)
mm = coeff[1]
fit = result.fun + ydata
rfit = result.fun
lxx = np.arange(-1.5,1.4,0.01) # controls the range of full fits to be shown
fullfit =  polyshape([coeff[1]],lxx) + laa
lxx = lxx + np.log10(bb)
#
#        error estimate; see pag. 98 of my notes.
#
hess = np.matmul(np.transpose(result.jac),result.jac)
ebb = rfit.std()/np.sqrt(hess[0,0])
emm = rfit.std()/np.sqrt(hess[1,1])
#
# -least squares fit to Sersic -- Einasto profiles.
#
def einastofit(coeff,xdata,ydata):
    rho0 = coeff[0]
    rprime = coeff[1]
    alpha = coeff[2]
    return np.log10(rho0*np.exp(-(xdata/rprime)**alpha))-ydata
#
rprime = 1.
rho0 = 1.
alpha = 1.
coeff0 = [rho0,rprime,alpha]
xxdata = 10**xdata
result = least_squares(einastofit, coeff0, args=(xxdata,ydata), method='lm')#,max_nfev = 1000)
coeff_einasto = result.x
fiteinasto = 10**(result.fun + ydata)
rfit_einasto = result.fun
rho0 = coeff_einasto[0]
rprime = coeff_einasto[1]
alpha = coeff_einasto[2]
fullfit_einasto = np.log10(rho0*np.exp(-(10**lxx/rprime)**alpha))
#
hess = np.matmul(np.transpose(result.jac),result.jac)
ealpha = rfit.std()/np.sqrt(hess[2,2])
nsersic = 1/alpha
ensersic = ealpha/alpha**2
#
#-----------------------------------------------
#
#-plots
fs = 18
plt.rcParams.update({'font.size': fs})
plt.ion()
plt.close(1) # to avoid getting multiple windows
fig, ax = plt.subplots(2,1, sharex = True, sharey = False,gridspec_kw={
                           'height_ratios': [6, 1]})
fig.subplots_adjust(hspace=0)
#
ax[0].set_yscale("log")
ymin = (10**ydata).min()/5.
ymax = (10**ydata).max()*2.
ax[0].set_ylim(ymin,ymax)
ax[0].set_ylabel(r"$\Sigma(R)\, [{\rm M}_\odot {\rm pc}^{-2}]$")
ax[0].set_xscale("log")
title ='Galaxy # '+'{:.0f}'.format(profno)+'\n'+r'$\log(M_\star/{\rm M}_\odot)$ = '+'{:.2f}'.format(mstar)
p = ax[0].plot(10**xdata,10**ydata,label=title,marker='o',ls='None')
color = p[0].get_color()
ax[0].plot(10**ldd,10**rho,label='Not used to fit',marker='o',ls='None',fillstyle='none',c=color)
ax[0].plot(10**xdata,rfit*1e-3,marker='o',c='r',ls='-',alpha=0.6,label='Residuals') # mock
cartel = '$m='+'{:.2f}$'.format(mm)+'$\pm$'+'{:.2f}'.format(emm)+'\n'+'$b=$'+'{:.2f}'.format(bb)+'$\pm$'+'{:.2f}'.format(ebb)+' kpc\n'+'RMS ='+'{:.3f}'.format(rfit.std())
p = ax[0].plot(10**xdata,10**fit,ls='-',lw=3,label=cartel)
color = p[0].get_color()
ax[0].plot(10**lxx,10**fullfit,lw=2,ls='--',c=color)
#
ax[0].plot(10**lxx,10**fullfit_einasto,lw=2,ls='--',label=r'Sérsic $n='+'{:.2f}$'.format(nsersic)+'\n'+'RMS ='+'{:.3f}'.format(rfit_einasto.std()))
#  
ax[0].legend(loc=3,fontsize=fs*0.6)#,title=r'PP in ${\bf '+gridname+'}$')
#
ax[1].set_xlabel(r"$R$ [kpc]")
ax[1].set_ylabel(r"[dex]",fontsize=fs*0.7)
ax[1].tick_params(direction='in',axis='y',labelsize=fs*0.8)
ax[1].yaxis.set_major_formatter(plt.FormatStrFormatter('%4.1f'))
ax[1].yaxis.tick_right()
ax[1].plot(10**xdata,rfit,marker='o',c='r',ls='-',alpha=0.6)
ax[1].plot([-10,100],[0,0],lw=1,c='k',ls='-',alpha=0.2)
ax[1].set_ylim(-0.2,.2)
#
profno = str(profno)
outfile = './plots/lane_emden_fit_'+gridname+'_'+profno+'.pdf'
print('output at '+outfile)
plt.savefig(outfile,transparent='False',bbox_inches='tight') 
plt.show()

