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
import sys, os
from math import pow
np.set_printoptions(threshold='nan')
    
## Calculate the formfactor of finite cylinders 
    ## Here a is the integrating parameter, integrated from 0 to pi/2
    ##      r is the radius, La is the apparent length, 
    ##      n is the power of the integrand: n=1 for calc beta correction factor, n=2 for calc pure P(q)
    
e0=8.85*pow(10,-12) #permittivity of space in SI units
ew=80.1 #dielectric constant of water
kbt=1.38*pow(10,-23)*(273.15+20) #Botlzmann constant x 20C RT, SI units  
ios=2*pow(10,-7)*pow((1.6022*pow(10,-19)),2.)*6.023*pow(10,26) #ionic strength of pH7 water, concentration is in #ions/m^3, hence avagadro's number
Ld=pow((8*np.pi*6.023*pow(10,23)*pow(1.6022,2)*pow(10,-38)/(ew*e0*kbt)),-0.5) #debye length using Mezzenga's expression, in meters
dsldsq=pow((6.393-2.5),2.)*pow(10,-12) #(squared SLD in ang^-4)
    
def AofQ_Cyl (a,r,La,q,n):   #   Amplitude function
    corefunc=(np.sin(q*La*np.cos(a)/2))/(q*La*np.cos(a)/2)*2*special.jv(1,(q*r*np.sin(a)))/(q*r*np.sin(a))
    func=(corefunc**float(n))*np.sin(a)
    return func

def PofQ_Cyl (r,La,q):       # P(Q) for a straight cylinder of radius R and length L
    integral1=quad(AofQ_Cyl,0,np.pi/2, args=(r,La,q,2))[0]
    return integral1

def beta (r,La,q):           # Correction factor while multiplying to the structure factor
    integral2=quad(AofQ_Cyl,0,np.pi/2, args=(r,La,q,1))[0]
    beta_func=(integral2**2.)/PofQ_Cyl(r,La,q)
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

def PofQ_PRISM (r,La,mu,v,vf,s,q):  #The P(Q) from interacting rigid rods using a corrhole and P(Q) for long cylinders
                                    # Based on polymer reference interaction site model (PRISM) theory (Schweizer & Curro, 1994
                                    # v is number/volume of the rods in #/m^3
    pCyl=np.vectorize(PofQ_Cyl)
    pthinCyl=np.vectorize(PofQ_long_thinCyl)
    corr= np.vectorize(corr_hole_func)
    sbundle=(1.+vf*special.j0(q*s*2*r))/2.
    PofQ= v*np.divide(pCyl(r,La,q),(1+np.multiply(corr(mu,v,La,q),pthinCyl(La,q))))*sbundle
    return PofQ 
    
def SofQ_fractal (q,La,e_cluster,df):    # Structure factor due to fractal formation 
                                         # in rods with apparent length L and fractal cluster size e_cluster
    frac_func = df*special.gamma(df-1)*np.sin((df-1)*np.arctan(q*e_cluster))/(pow((q*La),df)*pow((1+pow((q*e_cluster),-2.)),(df-1.)/2.))
    return frac_func

def IofQ_Cyl_corrhole_frac (q,*params): # For fractal forming straight cylinders
    s,vf,v,r,La,mu,e_cluster=params
    scale=1.
    df=fracdim
    #print(params)
    ff = np.vectorize (PofQ_PRISM)
    sf = np.vectorize (SofQ_fractal)
    betaf= np.vectorize (beta)
    I0=scale*(dsldsq*((np.pi*pow(r,2.)*La)**2.))*pow(10,-22)   #scale is a fudge factor,
                                                        #v is the number fraction in 1/m^3
                                                        #dsldq is the difference in sld squared, dimensions is Ang^-4
                                                        #The r and L are in angstroms from the fitting
                                                        #10^-10 is to convert the dimsions to cm^-1
    final_func=I0*np.multiply(ff(r,La,mu,v,vf,s,q), (1.+np.multiply(betaf(r,La,q),sf(q,La,e_cluster,df)))) + bkg
    return final_func

def Iqfitting (q,data,data_err,guessval,guesslow,guesshigh): # Fot fitting IofQ_Cyl_frac to data with initial guess values
    gfit_params=list(guessval)
    glow=list(guesslow)
    ghigh=list(guesshigh)
    print (gfit_params)
    Ivec=np.vectorize(IofQ_Cyl_corrhole_frac)
    opt_params,pcovs_params=curve_fit(f=Ivec,xdata=q,ydata=data,p0=gfit_params,sigma=data_err,absolute_sigma=True,
                                      bounds=([glow,ghigh]),x_scale=gscale0, verbose=2)
    unc=np.sqrt(np.diag(pcovs_params)) # standard deviation in parameters
    return opt_params,unc

def chisqred (q,y,yerr,pars): # Calculating chisqdred for the fit
    ratio= ((IofQ_Cyl_corrhole_frac(q,*pars)-y)/yerr)**2.
    chisqrd=np.sum(ratio)
    red_chisqrd=chisqrd/(len(y)-len(pars))
    return red_chisqrd

def DoPlot (q,y,yerr,pars): # For plotting the data and calculating fit
    fig=plt.figure()
    fig.suptitle(datafile, fontsize=14,weight="bold")
    fig.add_subplot(111)
    plt.xscale('log')
    plt.yscale('log')
    plt.errorbar(q,y-bkg,yerr, marker='o', c='r',label='data')
    Iqfit=IofQ_Cyl_corrhole_frac(q,*pars)
    plt.plot(q, Iqfit-bkg, '--', c='b', label='fit')
    plt.show()
    return q,Iqfit


        # Reading and fitting and plotting the data
print ("path is " + os.getcwd() + ", Correct?")
#datafolder=sys.argv[1] #Name of the folder 
#datafile=datafolder +"/ReducedFiles/"+ sys.argv[2] #Name of the file within ReducedFiles folder of the datafolder.
#Lfit=float(sys.argv[3])
datafile=sys.argv[1]
conc=float(sys.argv[2])
v1w=1*10/pow(10,-3)/(3600*80)*6.023*pow(10,23) #number density of 80 nm rods made of peptides of mw 3600g/mole
                            #v is in #/m^3, in the order of 2*10^22

Q, IQ, IQ_ERR, Q_ERR = np.loadtxt(datafile, dtype=float, skiprows=0, usecols=(0,1,2,3), unpack=True)

#bkgy=np.log(IQ[-25:,])
#bkgfit=np.polyfit(bkgx,bkgy,1)
bkg=np.average(IQ[-25:,])
print(bkg)

fracx=np.log(Q[1:12])
fracy=np.log(IQ[1:12])
fracslope,fracintercept=np.polyfit(fracx,fracy,1)
fracdim=-1*fracslope
print(fracdim,fracslope,fracintercept)
dimerfrac=0.5
dis=1.
scale=1 #dimesnionless
rad=10. #Ang
Lapp=80. #Ang
nD=conc*v1w*800./Lapp #number of rods/m^3
eclust=2*np.pi/Q[0]*100 #Ang
#fracdim=2. #dimensionless
cd= 9*pow(10,-20)/10 #charge density in ~5e/A

guessval0=np.array([dis,dimerfrac,nD,rad,Lapp,cd,eclust])
guesslow0=np.array([1.,0,0.1*nD,7,50,pow(10,-22),1000])
guesshigh0=np.array([1000,1,10*nD,11,Lapp*5,pow(10,1)*cd,pow(10,5)*eclust])
gscale0=[0.01,0.01,pow(10,21),0.1,0.1,pow(10,-22),100.]
fit_params0, unc_params0=Iqfitting(q=Q,data=IQ, data_err=IQ_ERR, 
                                 guessval=guessval0, guesslow=guesslow0,guesshigh=guesshigh0)


a,b,c,d,e,f,g=fit_params0 #scale,nD,rad,Lapp,cd,eclust,fracdim                                                                  
fitQ,fitIQ=DoPlot (Q,IQ,IQ_ERR,fit_params0)
chisqr=chisqred(Q,IQ,IQ_ERR,fit_params0)

guesslow1=[]
guesshigh1=[]
new_params=np.multiply(fit_params0,1.0)

fit_params1, unc_params1=Iqfitting(q=Q,data=IQ, data_err=IQ_ERR, 
                                 guessval=new_params, guesslow=guesslow0*0.95,guesshigh=guesshigh0*1.05)

#fit_params1=[3.632351795133535, 9.291749344976191, 597.6635095350161, 39.94262285274863, 0.7646523514717469, 254.6087368034978, 264319679.43412867, 2.522696194594406]
#unc_params1=[0.21896201455628636, 0.048089906757187745, 36.30298072506388, 2.5352126011459197, 0.06084742997808868, 11.090816128293195, 3.1675851532098686e-11, 0.07280736723684778]
a,b,c,d,e,f,g=fit_params1 #scale,nD,rad,Lapp,cd,eclust,fracdim
fitQ,fitIQ=DoPlot (Q,IQ,IQ_ERR,fit_params1)
chisqr=chisqred(Q,IQ,IQ_ERR,fit_params1)
I0fitted=(dsldsq*((np.pi*pow(c,2)*d)**2.))*pow(10,-22)*b
        # Writing to file
filename="Fit_bundle_"+sys.argv[1]
f=open(filename,'w')
f.write("\n Given params: fracdim = %.2e, \t"%fracdim)
f.write("\n Bkg value= %.3e, \t"%bkg)
f.write("\n Given guess values: %.2e,%.2e,%.2e,%.2e,%.2e,%2e"%(dis,nD,rad,Lapp,cd,eclust))
f.write("\n Fitted params interparticle distance,number density v(#/m3), radius r(A),\napparent length Lapp(A),charge density cd(C/A), cluster size eclust(A)\n")
for l in range(len(fit_params1)):
    f.write("%.3e +/- %.3e \n"% (fit_params1[l],unc_params1[l]))
f.write("\n I0fitted = %.5e \n"%I0fitted)
f.write("\n Reduced chi-squared = %f\n"%chisqr)
f.write("\n q\t I(Q)\t error\t fit\t nor_IQ\t nor_fit\n\n")
for i in range(len(Q)):
    f.write ("%f\t%f\t%f\t%f\t%f\t%f\n"% (fitQ[i], IQ[i]-bkg, IQ_ERR[i], fitIQ[i]-bkg, (IQ[i]-bkg)/I0fitted, (fitIQ[i]-bkg)/I0fitted))
f.close()
    




    



    
    


