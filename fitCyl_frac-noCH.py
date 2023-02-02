# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 12:47:35 2018

@author: njs2
"""

from scipy.integrate import quad
from scipy import special
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import sys#, os
np.set_printoptions(threshold='nan')
    
## Calculate the formfactor of finite cylinders 
    ## Here a is the integrating parameter, integrated from 0 to pi/2
    ##      r is the radius, L is the legnth, 
    ##      n is the power of the integrand, n=1 for calc beta correction factor, n=2 for calc pure P(q)
    
def AofQ_Cyl (a,r,L,q,n):   #   Amplitude function
    corefunc=(np.sin(q*L*np.cos(a)/2))/(q*L*np.cos(a)/2)*special.jv(1,(q*r*np.sin(a)))/(q*r*np.sin(a))
    func=corefunc**float(n)*np.sin(a)
    return func

def PofQ_Cyl (r,L,q):       # P(Q) for a straight cylinder of radius R and length L
    integral1=quad(AofQ_Cyl,0,np.pi/2, args=(r,L,q,2))[0]
    return integral1

def beta (r,L,q):           # Correction factor while multiplying to the structure factor
    integral2=quad(AofQ_Cyl,0,np.pi/2, args=(r,L,q,1))[0]
    beta_func=integral2**2./PofQ_Cyl(r,L,q)
    return beta_func
    
def SofQ_fractal (q,e_pore,e_cluster,d): # Structure factor due to fractal formation 
                                         # in hard spheres of radius e_pore and cluster size e_clsuter
    frac_func = 1 + d*special.gamma(d-1)*np.sin((d-1)*np.arctan(q*e_cluster))/[(q*e_pore)**d*(1+(q*e_cluster)**(-2))**(d-1)/2]
    return frac_func

def IofQ_Cyl_corrhole_frac (q,*params): # For fractal forming straight cylinders
    I0,r,L,rb,v,e_pore,e_cluster,d=params
    ff = np.vectorize (PofQ_Cyl)
    sf = np.vectorize (SofQ_fractal)
    betaf= np.vectorize (beta)
    #print (ff,sf,betaf)
    final_func=I0*np.multiply(ff(r,L,q), (1.+betaf(r,L,q)*(sf(q,e_pore,e_cluster,d)-1.)))
    return final_func

def Iqfitting (q,data,data_err,guessval,guesslow,guesshigh): # Fot fitting IofQ_Cyl_frac to data with initial guess values
    gfit_params=list(guessval)
    glow=list(guesslow)
    ghigh=list(guesshigh)
    print (gfit_params)
    Ivec=np.vectorize(IofQ_Cyl_corrhole_frac)
    opt_params,pcovs_params=curve_fit(f=Ivec,xdata=q,ydata=data,p0=gfit_params,sigma=data_err,absolute_sigma=True,
                                      bounds=([glow,ghigh]), x_scale=gscale0, verbose=2)
    unc=np.sqrt(np.diag(pcovs_params)) # standard deviation in parameters
    return opt_params,unc

def chisqred (q,y,yerr): # Calculating chisqdred for the fit
    ratio= ((IofQ_Cyl_corrhole_frac(q,*fit_params1)-y)/yerr)**2.
    chisqrd=np.sum(ratio)
    red_chisqrd=chisqrd/(len(y)-len(fit_params1))
    return red_chisqrd

def DoPlot (q,y,yerr): # For plotting the data and calculating fit
    fig=plt.figure()
    fig.suptitle(datafile, fontsize=14,weight="bold")
    fig.add_subplot(111)
    plt.xscale('log')
    plt.yscale('log')
    plt.errorbar(q,y,yerr, marker='o', c='r',label='data')
    Iqfit=IofQ_Cyl_corrhole_frac(q,*fit_params1)
   #Iqfit=IofQ_Cyl_frac(q,*fit_params1)
    plt.plot(q, Iqfit, '--', c='b', label='fit')
    plt.show()
    return q,Iqfit

        # Reading and fitting and plotting the data
datafile=sys.argv[1]
Q, IQ, IQ_ERR = np.loadtxt(datafile, dtype=float, skiprows=3, usecols=(0,1,2), unpack=True)
guessval0=np.array([2.5,9.2,400,150,10000.,2.50])
guesslow0=np.array([0.,8.,40,40,900,1])
guesshigh0=np.array([100,15,2000,np.inf, np.inf, 5.00 ])
gscale0=[0.00001,1,1,1, 1, 0.01]
fit_params0, unc_params0=Iqfitting(q=Q,data=IQ, data_err=IQ_ERR, 
                                 guessval=guessval0, guesslow=guesslow0,guesshigh=guesshigh0)
guesslow1=[]
guesshigh1=[]
new_params=np.multiply(fit_params0,0.95)
for v in range(len(unc_params0)):
    lowval=new_params[v]-unc_params0[v]
    guesslow1.append(lowval)
    highval=fit_params0[v]+unc_params0[v]
    guesshigh1.append(highval)

fit_params1, unc_params1=Iqfitting(q=Q,data=IQ, data_err=IQ_ERR, 
                                 guessval=new_params, guesslow=guesslow1,guesshigh=guesshigh1)

#fit_params1=[3.632351795133535, 9.291749344976191, 597.6635095350161, 39.94262285274863, 0.7646523514717469, 254.6087368034978, 264319679.43412867, 2.522696194594406]
#unc_params1=[0.21896201455628636, 0.048089906757187745, 36.30298072506388, 2.5352126011459197, 0.06084742997808868, 11.090816128293195, 3.1675851532098686e-11, 0.07280736723684778]
a,b,c,d,e,f=fit_params1
fitQ,fitIQ=DoPlot (Q,IQ,IQ_ERR)
chisqr=chisqred(Q,IQ,IQ_ERR)

        # Writing to file
filename="Fit_noCH_"+datafile
f=open(filename,'w')
f.write("Printing parameters I0 (/cm), R (A), L(A), e_pore (A), e_cluster (A), d (-) and their standardv deviations" )
f.write("First round of fits for paramters\n")
#for k in range(len(fit_params0)):
#    f.write("%f +/- %f\n"% (fit_params0[k],unc_params0[k]))
f.write("\nFinal fits\n")
for l in range(len(fit_params1)):
    f.write("%f +/- %f\n"% (fit_params1[l],unc_params1[l]))
f.write("Reduced chi-squared = %f\n"%chisqr)
f.write("\nq\tI(Q)\terror\tfit\n\n")
for i in range(len(Q)):
    f.write ("%f\t%f\t%f\t%f\n"% (fitQ[i], IQ[i], IQ_ERR[i], fitIQ[i]))
f.close()
    




    



    
    


