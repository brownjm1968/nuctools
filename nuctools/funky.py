import pandas as pd
import numpy as np
import scipy.integrate as intgr
import time
from . import physics_tools as pt

# run all functions?
RUNALL = False

__all__ = ["c_edge","threshFinder",
           "decayFinder","six","four","two","sint",'ten','ten_exp',
           "gauss","noisy_gauss","egrp_avg","giveIt","areal_density"]

mnc2 = 939565413.3  # NIST neutron mass [eV]
c = 299.792458      # NIST speed of light [m/us]

###### Functions ##############################################################

def c_edge(Eg):
    """
    Compute the Compton edge energy associatedwith a given 
    gamma energy
    
    Parameters
    ----------
    Eg : float
        Energy of gamma that compton scatters
    
    Returns
    -------
    Ec : float
        Energy of Compton edge

    Notes
    -----
    .. math:: E_c(E_\\gamma) = \\frac{2.0 E_\\gamma^2}{511.0+2.0E_\\gamma}
    """
    Ec = 2.0*Eg**2/(511.0+2.0*Eg)
    return Ec

#------------------------------------------------------------------------------
# input: pulse vector (floats), threshold where pulse
#        increases above 20% for example
# output: index where threshold is FIRST crossed
def threshFinder(pulse,threshold):
    """
    For a positive "pulse", define a threshold and return the 
    index immediately after that threshold is passed.

    Parameters
    ----------
    pulse: array_like
        1-d set of values that make up the ordinates of a 
        positive pulse
    threshold : float
        A value > the minimum pulse value, and < the max pulse
        value
    
    Returns
    -------
    index : int
        The index of the pulse value immediately after the 
        pulse value is greater than the threshold

    """
    index = 0
    for i in range(len(pulse)):
        if pulse[i] > threshold:
            index = i
            break
    return index
#------------------------------------------------------------------------------
# input: pulse vector (floats), threshold where pulse
#        decreases below 20% for example
# output: index where threshold is SECOND crossed
def decayFinder(pulse,threshold):
    """
    For a positive "pulse", define a threshold and return the 
    index immediately after that threshold is passed for a 
    second time. 

    Parameters
    ----------
    pulse: array_like
        1-d set of values that make up the ordinates of a 
        positive pulse
    threshold : float
        A value > the minimum pulse value, and < the max pulse
        value
    
    Returns
    -------
    index : int
        The index of the pulse value immediately after the 
        pulse value is < than the threshold, if the 
        threshold has already been passed once.
    """
    index = 0
    flag = False
    peakIndex = max(range(len(pulse)), key=pulse.__getitem__)
    for i in range(len(pulse)):
        if (flag == True) and (pulse[i] < threshold) and (i > peakIndex):
            index = i
            break
        if pulse[i] > threshold:  # flag the first crossing
            flag = True
    return index
#------------------------------------------------------------------------------
# string formatter
def six(x):                     # (6 digits . 6 digits) float
    return "{:6.6f}".format(x)
def four(x):                     # (6 digits . 4 digits) float
    return "{:6.4f}".format(x)
def two(x):                     # (6 digits . 2 digits) float
    return "{:6.2f}".format(x)
def sint(x):                    # int
    return "{:d}".format(x)
def ten(x):
    return "{:9.5f}".format(x) # (10 digits total, 4 digits after dec.)
def ten_exp(x):
    return "{:9.3e}".format(x) # (10 digits total, 4 digits after dec.)
#------------------------------------------------------------------------------

def gauss(x,mu,sig):
    """
    Generate a Gaussian distribution

    Gaussian is generated with user input average and 
    standard deviation

    Parameters
    ----------
    x : numpy array 
        Vector of absiccae to plot the gaussian along
    mu : float
        The average value of the gauss, where the peak
        will occur.
    sig : float
        The standard deviation from the average, varying
        the width of the "bell" curve.
    
    Returns
    -------
    gauss : numpy array
        The ordinates of same length x, which give the
        Gaussian or normal distribution for the specified
        parameters.
    
    Examples
    --------
    >>> x = np.linspace(0,100,100)
    >>> mu,sig = 50,10
    >>> gauss1 = gauss(x,mu,sig)
    >>> plt.figure(3)
    >>> plt.plot(x,gauss1,color="k",ls=" ",marker=".")
    >>> plt.show()

    """
    # gauss curve
    gauss = 1/(2*np.pi*sig**2)**0.5*np.exp(-(x-mu)**2/(2*sig**2))
    
    return gauss

#------------------------------------------------------------------------------

def noisy_gauss(x,mu,sig):
    """
    Generate a Gaussian distribution with noise.

    Gaussian is generated with user input average and 
    standard deviation, noise is added to all data pts.
    with the numpy random package.

    Parameters
    ----------
    x : numpy array 
        Vector of absiccae to plot the gaussian along
    mu : float
        The average value of the gauss, where the peak
        will occur.
    sig : float
        The standard deviation from the average, varying
        the width of the "bell" curve.
    
    Returns
    -------
    gauss : numpy array
        The ordinates of same length x, which give the
        Gaussian or normal distribution for the specified
        parameters.
    
    Examples
    --------
    >>> x = np.linspace(0,100,100)
    >>> mu,sig = 50,10
    >>> gauss = noisy_gauss(x,mu,sig)
    >>> plt.figure(3)
    >>> plt.plot(x,gauss,color="k",ls=" ",marker=".")
    >>> plt.show()

    """
    # gauss curve
    ngauss = gauss(x,mu,sig)+np.random.rand(len(x))*0.01
    
    return ngauss

#------------------------------------------------------------------------------
def egrp_avg(E,Eg,data,err_data):
    """
    Group quantities and average them along specified
    energy grid Eg.

    Grouping is done with integration by scipy quad()
    function and division by group size.

    Parameters
    ----------

    E : array_like
        Energy grid data is currently on
    Eg : array_like
        Grouped energy grid
    data : array_like
        Quantity to be group averaged
    err_data : array_like 
        The error on the data array. Must be same
        length as data array.

    Returns
    -------

    gdata : pandas DataFrame
        The grouped data is attribute 'group_avg', and the
         **estimated** absolute error on the grouped data
         is attribute 'err_group_avg'.

    Notes
    -----
    The estimated error on the grouped data is simply the error
    on a simple average on the group. The actual grouping is done
    by the scipy.integrate.quad function.

    .. math:: \\frac{\\Delta \\langle x\\rangle}{\\langle x\\rangle}= \\frac{\\sqrt{\\sum^N (\\Delta x )^2}}{\\sum^Nx}  

    """
    # cast inputs as numpy arrays
    E,Eg,data = np.array(E),np.array(Eg),np.array(data)
    # if E is not in ascending order throw error
    if any(np.diff(E) < 0) or any(np.diff(Eg) < 0):
        raise ValueError("E,Eg must be in ascending order.")
    # create grid with 20 times the energy points (linear)
    interpE = np.linspace(0,max(E),len(E)*20)
    # create function to integrate over
    integ = lambda interpE: np.interp(interpE,E,data)
    # make energy group average array 
    Eg_avg = np.zeros(len(Eg))
    rel_err = np.zeros(len(Eg))
    # ---------------------------------------
    # --- error estimation ------------------
    # ---------------------------------------
    df = pd.DataFrame({'E':E,'data':data,'err_data':err_data})
    # for each energy group, integrate and normalize
    for i in range(len(Eg)):
        Eg_avg[i] = intgr.quad(integ,Eg[i-1],Eg[i])[0]/(Eg[i]-Eg[i-1])
        # --- error ----------
        squared_error = df.err_data[(df.E>Eg[i-1]) & (df.E<Eg[i])]**2
        simple_grp = np.sum(df.data[(df.E>Eg[i-1]) & (df.E<Eg[i])])
        rel_err[i] = np.sqrt(squared_error.sum())/simple_grp


    return pd.DataFrame({'group_avg':Eg_avg,'err_group_avg':rel_err*Eg_avg})

#------------------------------------------------------------------------------
def frohner_cor(sig1,sig2,n1,n2):
    """
    Takes cross-sections [barns] and atom densities [atoms/barn] for 
    two thicknesses of the same sample, and returns extrapolated cross
    section according to Frohner.

    Parameters
    ----------
    sig1 : array_like
        Cross section of the thinner of the two samples.
    sig2 : array_like
        Cross section of the thicker of the two samples.
    n1 : float
        Atom density of the thinner sample
    n2 : float 
        Atom density of the thicker sample

    Returns
    -------
    sig0 : array_like
        The extrapolated cross section from sig1 and sig2
    """

    return (n2*sig1-n1*sig2)/(n2-n1)

#------------------------------------------------------------------------------
def frohner_cor_3rd_order(sig1,sig2,sig3,n1,n2,n3):
    """
    Takes cross-sections [barns] and atom densities [atoms/barn] for 
    three thicknesses of the same sample, and returns extrapolated 
    cross section according to Frohner.

    Parameters
    ----------
    sig1 : array_like
        Cross section of the thinnest of the three samples.
    sig2 : array_like
        Cross section of the mid-thickness of the three samples.
    sig3 : array_like
        Cross section of the thickest of the three samples.
    n1 : float
        Atom density of the thinnest sample
    n2 : float 
        Atom density of the mid-thickness sample
    n3 : float
        Atom density of the thickest sample

    Returns
    -------
    sig0 : array_like
        The extrapolated cross section from sig1, sig2, and sig3
    """

    # two terms in the numerator
    numer1 = (n1*sig2-n2*sig1)*(n3**2-n1**2-(n1-n3)/(n1-n2)*(n2**2-n1**2))
    numer2 = (n1*n2**2-n1**2*n2)*(sig3-sig2-(n1-n3)/(n1-n2)*(sig2-sig1))
    denom = (n1-n2)*(n3**2-n1**2) - (n1-n3)*(n2**2-n1**2)
    
    return (numer1-numer2)/denom

def areal_density(diam,unc_diam,mass,unc_mass,molar_mass):
    """
    Calculate areal density [at/b] of a sample based on sample
    diameter, mass, and molar mass. Assumes a disc-shaped sample.

    Does not consider isotopic abundances or multiple elements

    Parameters
    ----------
    diam : float
        The diameter of the sample [cm]
    unc_diam : float
        The uncertianty in the sample diameter [cm]
    mass : float
        The mass of the sample [g]
    unc_mass : float
        The uncertainty in the mass measurement of the sample [g]
    molar_mass : float
        The molar mass for the element/isotope

    Returns
    -------
    atoms_per_barn : float
        The sample thickness/areal-density [at/b]
    dapb : float
        The uncertainty in areal-density [at/b]
    """
    cm2pb = 1e-24 # [cm^2/b]

    area = (diam/2)**2*np.pi
    dadd = diam/2*np.pi       # partial deriv. w.r.t. diam
    darea = np.sqrt( ( unc_diam*dadd )**2 )

    avo = pt.Na # avogadros number

    atoms_per_barn = avo/molar_mass*mass/area*cm2pb # [at/mol][mol/g][g][1/cm^2][cm^2/b] = at/b
    dapbdarea = -avo/molar_mass*mass*cm2pb/area**2
    dapbdmass = avo/molar_mass/area*cm2pb
    dapb = np.sqrt( ( darea*dapbdarea )**2 + ( unc_mass*dapbdmass )**2 )

    return atoms_per_barn,dapb 


#############################################################

#############################################################

def giveIt():
    return 900
