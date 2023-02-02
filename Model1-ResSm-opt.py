# -*- coding: utf-8 -*-
"""
Created on Sun Dec 16 21:01:27 2018

@author: njs2
"""

from scipy.integrate import quad 
from scipy import special
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import sys, os
from math import pow
np.set_printoptions(threshold='nan')
    
## For correlation hole these constants are evaluated in SI units
    
e0=8.85*pow(10,-12) #permittivity of space in SI units
ew=80.1 #dielectric constant of water
kbt=1.38*pow(10,-23)*(273.15+20) #Botlzmann constant x 20C RT, SI units  
ios=2*pow(10,-7)*pow((1.6022*pow(10,-19)),2.)*6.023*pow(10,26) #ionic strength of pH7 water, concentration is in #ions/m^3, hence avagadro's number
Ld=pow((8*np.pi*6.023*pow(10,23)*pow(1.6022,2)*pow(10,-38)/(ew*e0*kbt)),-0.5) #debye length using Mezzenga's expression, in meters
dsldsq=pow((6.393-2.5),2.)*pow(10,-12) #(squared SLD in ang^-4)

## Calculate the formfactor of finite cylinders 
    ## Here a is the integrating parameter, integrated from 0 to pi/2
    ##      r is the radius, La is the apparent length, 
    ##      n is the power of the integrand: n=1 for calc beta correction factor, n=2 for calc pure P(q)    
    
def AofQ_Cyl (a,r,La,q,n):   #   Amplitude function
    corefunc=(np.sin(q*La*np.cos(a)/2))/(q*La*np.cos(a)/2)*2*special.jv(1,(q*r*np.sin(a)))/(q*r*np.sin(a))
    func=(corefunc**float(n))*np.sin(a)
    return func

def PofQ_Cyl (r,La,q):       # P(Q) for a straight cylinder of radius R and length L
    integral1=quad(AofQ_Cyl,0,np.pi/2, args=(r,La,q,2))[0]
    return integral1

def betanum (r,La,q):           # Correction factor while multiplying to the structure factor
    integral2=quad(AofQ_Cyl,0,np.pi/2, args=(r,La,q,1))[0]
    beta_func=(integral2**2.)
    return beta_func

def si_func(t):             # The Si(x) function for calculating scattering from infinitely thin rod
    corefunc1=(1/t)*np.sin(t)
    return corefunc1

def PofQ_long_thinCyl (La,q):  # The P(Q) for rods L>>r dominated by scattering from the length alone
    par= np.multiply(q,La)
    si_par=quad(si_func, 0, par)[0]
    PofQ=2*si_par/(par) - 4*(np.sin(par/2.))**2./(par**2.)
    return PofQ
        
def corr_hole_func (mu,v,La,q):      # The Schneider's expression to account for correltn hole, charge density mu, Length La of rods 
    qadj=q*pow(10,10)                # convert q from Ang-1 to meters-1
    corr_func=(v*pow(mu*La,2.))/(pow(qadj,2.)+pow(Ld,-2.))/(kbt*ew*e0) #dimesnionless
    #print(q,corr_func)
    return corr_func

def PofQ_corr_hole (r,La,mu,v,q):  #The P(Q) from interacting rigid rods using a corrhole and P(Q) for long cylinders
                                    # Based on polymer reference interaction site model (PRISM) theory (Schweizer & Curro, 1994
                                    # v is number/volume of the rods in #/m^3
    pCyl=np.vectorize(PofQ_Cyl)
    pthinCyl=np.vectorize(PofQ_long_thinCyl)
    corr= np.vectorize(corr_hole_func)
    PofQ= v*np.divide(pCyl(r,La,q),(1+np.multiply(corr(mu,v,La,q),pthinCyl(La,q))))
    return PofQ

def SofQ_fractal (q,La,e_cluster,df):    # Structure factor due to fractal formation 
                                         # in rods with apparent length L and fractal cluster size e_cluster
    frac_func = df*special.gamma(df-1)*np.sin((df-1)*np.arctan(q*e_cluster))/(pow((q*La),df)*pow((1+pow((q*e_cluster),-2.)),(df-1.)/2.))
    return frac_func

def IofQ_Cyl_corrhole_frac (q,*params): # For fractal forming straight cylinders
    r,v,La,mu,df=params #volume frac, radius, length and cluster size
    df=fracdim
    e_cluster=eclust
    #v=conc*v1w*800./La #number of rods/m^3
    #print(params)
    ff = np.vectorize (PofQ_corr_hole)
    sf = np.vectorize (SofQ_fractal)
    betaf= np.vectorize (betanum)
    I0=(dsldsq*((np.pi*pow(r,2)*La)**2.))*pow(10,-22)   #v is the number fraction in 1/m^3--> required to be fit for good convergence
                                                        #dsldq is the difference in sld squared, dimensions is Ang^-4
                                                        #The r and L are in angstroms from the fitting
                                                        #10^-10 is to convert the dimensions to cm^-1    
    pos=np.where(Q==q)[0][0]
    if rq0qm[pos][0]==q:    
        qj=rq0qm[pos][1]
        smearqj=rq0qm[pos][2]
        if len(qj)==len(smearqj):
            pq=ff(r,La,mu,v,qj)
            beta=betaf(r,La,qj)/pq
            I_mod=np.multiply(pq, (1.+np.multiply(beta,sf(qj,La,e_cluster,df))))
            I_mod_res=np.multiply(I_mod,smearqj)
            Iqs=np.sum(I_mod_res)
           #print(Iqs)
    else:
        print("Something is wrong with the index!!")
    final_func=I0*Iqs
    return final_func

def Iqfitting (q,data,data_err,guessval,guesslow,guesshigh): # Fot fitting IofQ_Cyl_frac to data with initial guess values
    gfit_params=list(guessval)
    glow=list(guesslow)
    ghigh=list(guesshigh)
    print (gfit_params)
    Ivec=np.vectorize(IofQ_Cyl_corrhole_frac)
    opt_params,pcovs_params=curve_fit(f=Ivec,xdata=q,ydata=data,p0=gfit_params,sigma=data_err,absolute_sigma=True,
                                      bounds=([glow,ghigh]), verbose=2, xtol=1e-8,gtol=1e-8,ftol=1e-8, x_scale=gscale0)
    unc=np.sqrt(np.diag(pcovs_params))# standard deviation in parameters
    return opt_params,unc

def chisqred (q,y,yerr,Iq,pars):#Calculating chisqdred for the fit
    ratio= ((Iq-y)/yerr)**2.
    chisqrd=np.sum(ratio)
    red_chisqrd=chisqrd/(len(y)-len(pars))
    return red_chisqrd

def DoPlot (q,y,yerr,pars): # For plotting the data and calculating fit
    fig=plt.figure()
    fig.suptitle(datafile, fontsize=14,weight="bold")
    fig.add_subplot(111)
    plt.xscale('log')
    plt.yscale('log')
    plt.errorbar(q,y,yerr, marker='o', c='r',label='data')
    Ivec=np.vectorize(IofQ_Cyl_corrhole_frac)
    Iqfit=Ivec(q,*pars)
    plt.plot(q, Iqfit, '--', c='b', label='fit')
    plt.show()
    return q,Iqfit

        # Reading and fitting and plotting the data

print ("path is " + os.getcwd() + ", Correct?\n\n")
datafile=sys.argv[1]
conc=float(sys.argv[2])
v1w=1*10/pow(10,-3)/(3600*80)*6.023*pow(10,23) #number density of 80 nm rods made of peptides of mw 3600g/mole
                            #v is in #/m^3, in the order of 2*10^22

Q, IQ, IQ_ERR, SIGQ, QMEAN, FS = np.loadtxt(datafile, dtype=float, skiprows=0, usecols=(0,1,2,3,4,5), unpack=True)

def gauss_smear (q0,qmean,sigq,fs):      #For a given qmean, sigq and fs corresponding to a value q0 in the experimental list of q values,
    #print ( q0, qmean, sigq, fs )      #this calculates the smear contribution R(q,qmean) for all experimentally measured q value.
    dq=sigq/5.
    qj=np.arange(qmean-3.*sigq,qmean+3.*sigq,dq)
    smear=fs/((2*np.pi*(sigq**2.))**0.5)*np.exp(-((qj-qmean)**2.)/(2*(sigq**2.)))*dq
    for i in range(len(qj)):
        #print(qj[i],smear[i])
      # print('\n')
        if qj[i]<0:
            qj[i]=0
            smear=fs/(2*np.pi*sigq**2.)**0.5*np.exp(-(-qmean)**2./(2*sigq**2.))*dq
            print ("Negative value encountered! Changed to 0")
    #print ('\n')
    return q0,qj,smear    

##Create the array that lists in col 0=Q, col 1= list of qj for Q, col 2 = list of smear values for qj

rq0qm=[]
for i in range(len(Q)):
    r_val=gauss_smear(Q[i],QMEAN[i],SIGQ[i],FS[i])
    rq0qm.append(r_val)
    
## Fit a line to background at high Q
    
bkgQ=np.log10(Q[-25:,])
bkg_IQ=np.log10(IQ[-25:,])
bkg_val,bkg_cov=np.polyfit(bkgQ,bkg_IQ,1,cov=True)
bkg_slope=bkg_val[0]
bkg_int=bkg_val[1]
print(bkg_slope,bkg_int)

## Subtract background from IQ

IQ_n=IQ-np.power(10.,(bkg_slope*np.log(Q)+bkg_int))
IQ_ERR_n=IQ_ERR

## Fit fractal dimension to low Q

fracx=np.log10(Q[1:10])
fracy=np.log10(IQ[1:10])
val,cova=np.polyfit(fracx,fracy,1, cov=True)
fracdim=-1*val[0]
fracdim_err=np.sqrt(cova[0][0])
low_frac=fracdim-fracdim_err*3
high_frac=fracdim+fracdim_err*3
eclust=pow(10.,4)
print(fracdim,fracdim_err,eclust)

## Define initial values and guesses

rad=10.0 #Ang
Lapp=85 #Ang always greater than 40 A (length of bundle)
nD=conc*v1w*800./Lapp #number of rods/m^3
cd= 1.6*pow(10,-20) #charge density in ~ 0.75 e/nm in e/A

guessval0=np.array([rad,nD,Lapp,cd,fracdim])
guesslow0=np.array([8,0.5*nD,70,1.6*pow(10,-21),low_frac])
guesshigh0=np.array([12,2*nD,140,1.6*pow(10,-19),high_frac])
gscale0=[0.1,2*pow(10,23),1,pow(10,-23),0.1]
fit_params0, unc_params0=Iqfitting(q=Q,data=IQ_n, data_err=IQ_ERR_n, 
                                 guessval=guessval0, guesslow=guesslow0,guesshigh=guesshigh0)


b,c,d,e,f=fit_params0 #rad,nd,length,chargeDen,fractDim                                                                 
fitQ0,fitIQ0=DoPlot(Q,IQ_n,IQ_ERR_n,fit_params0)
chisqr0=chisqred(Q,IQ_n,IQ_ERR_n,fitIQ0,fit_params0)

new_params=np.multiply(fit_params0,1)
fit_params1, unc_params1=Iqfitting(q=Q,data=IQ_n, data_err=IQ_ERR_n, 
                                 guessval=new_params, guesslow=fit_params0*0.9,guesshigh=fit_params0*1.1)

b,c,d,e,f=fit_params1 #b=rad,c=v,d=length,e=chargeDen,f=fractDim

fitQ,fitIQ=DoPlot (Q,IQ_n,IQ_ERR_n,fit_params1)
chisqr=chisqred(fitQ,IQ_n,IQ_ERR_n,fitIQ,fit_params1)
I0fitted=(dsldsq*((np.pi*pow(b,2)*d)**2.))*pow(10,-22)*c
final_conc=c*d/800/v1w

        # Writing to file
filename="ResSmr-Fit-"+sys.argv[1]
f=open(filename,'w')
f.write("\n Given params: fracdim = %.3e +/- %.3e \t"%(fracdim,fracdim_err))
f.write("\n \t ecluster = %.3e \t"%eclust)
#f.write("\n \t Bkg value= %.3e +/- %.3e, \t"%(bkg,bkg_std))
f.write("\n Given guess values: %.2e,%.2e,%.2e,%2e,%2e"%(nD,rad,Lapp,cd,fracdim))
f.write("\n Fitted params radius r(A),number density v(#/m3), apparent length Lapp(A), charge desnity (e/A), fractalDimens () \n")
f.write("\n I round results:\n")
for ll in range(len(fit_params0)):
    f.write("%.3e +/- %.3e \n"% (fit_params0[ll],unc_params0[ll]))
f.write("\n II round results (final):\n")
for l in range(len(fit_params1)):
    f.write("%.3e +/- %.3e \n"% (fit_params1[l],unc_params1[l]))
f.write("\n I0fitted = %.4e \n"%I0fitted)
f.write("\n fitted conc = %.4e \n"%final_conc)
f.write("\n Reduced chi-squared =\n1. %f \n2. %f\n"%(chisqr0,chisqr))
f.write("\n q\t I(Q)\t error\t fit\t nor_IQ\t nor_fit\n\n")
for i in range(len(Q)):
    f.write ("%f \t%f \t%f \t%f \t%0.3e \t%0.3e \n"% (Q[i], IQ_n[i], IQ_ERR_n[i], fitIQ[i], (IQ_n[i])/I0fitted, (fitIQ[i])/I0fitted))
f.close()