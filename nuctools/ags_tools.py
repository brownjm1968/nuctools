import numpy as np
import pandas as pd

__all__ = ['ags','e1p0','e2p0','e3p0']

class ags:
    """
    Python class to contain AGS-like functions and information

    Attributes
    ----------
    tof : array-like
        A 1-d array that is a Pandas Series of time-of-flight values,
        typically given in [ns]
    counts : array-like
        A 1-d array that is a Pandas Series of binned counts
    cps : array-like
        A 1-d array that is a Pandas Series of count rate given in 
        counts per second
    energy : array-like
        A 1-d array that is a Pandas Series of energy values transformed 
        from the time-of-flight values

    """
    def __init__(self):
        self.tof = None
        self.counts = None
        self.cps = None
        self.energy = None
        
    def read_grouped_counts(self,spectrum,comp_pt,comp_fct,binsize):
        """
        Read in the grouped counts file from AGL and populate the 
        counts, tof, and cps attributes of the ags class
        
        Parameters
        ----------
        filename : str
            The full file path to the AGL grouped counts file
        comp_pt : array-like
            Compression points given in bin numbers (integers). Must
            be of length == len(comp_fct)+1
        comp_fct : array-like
            Compression factors for each group. These are specified in 
            a separate text file from AGL. Typically in powers of 2 
            (2^N).
        binsize : float
            The width of the base bin in [ns]

        Returns
        -------
        nothing : None
            Populates the attributes of the class: tof, counts, cps

        Examples
        --------

        Notes
        -----

        """
        cp = np.array(comp_pt)
        cf = 2**np.array(comp_fct)
        #hist = pd.read_csv(filename,delim_whitespace=True,names=['bin','counts'])
        hist = pd.DataFrame({
            'bin'    : spectrum[0],
            'counts' : spectrum[1],
            })
        hist['cps'],hist['tof'] = 0,0
        
        for i in range(len(cp)-1):
            hist.loc[cp[i]:cp[i+1]-1,'cps'] = hist.counts[cp[i]:cp[i+1]]/(cf[i])
            if i==0:
                hist.loc[cp[i]:cp[i+1]-1,'tof'] = 0+np.arange(cp[i+1]-cp[i])*binsize*cf[i]
            else:
                hist.loc[cp[i]:cp[i+1]-1,'tof'] = hist.tof[cp[i]-1]+np.arange(cp[i+1]-cp[i])*binsize*cf[i]
                
        self.tof = hist.tof
        self.cps = hist.cps
        self.counts = hist.counts

def e1p0(tof,p1,p2,p3):
    """
    Background function for TOF spectra

    Parameters
    ----------
    tof : array-like
        The time-of-flight spectrum
    p1 : float
        constant background
    p2 : float
        multiplier on 1st exponential
    p3 : float
        multiplier on time-of-flight in 1st exponent
    p4 : float
        constant added to 1st exponent

    Returns
    -------
    e1p0 : array-like
        The function in the length of t (see notes)

    Notes
    -----
    .. math:: f(t) = p1 + p2e^{p3t+p4}
    """
    return p1 + p2*np.exp(p3*tof)

def e2p0(tof,p1,p2,p3,p5,p6):
    """
    Background function for TOF spectra

    Background function for TOF spectra

    Parameters
    ----------
    tof : array-like
        The time-of-flight spectrum
    p1 : float
        constant background
    p2 : float
        multiplier on 1st exponential
    p3 : float
        multiplier on time-of-flight in 1st exponent
    p4 : float
        constant added to 1st exponent
    p5-p7 : float
        (see equation in notes)

    Returns
    -------
    e2p0 : array-like
        The function in the length of t (see notes)

    Notes
    -----
    .. math:: f(t) = p1 + p2e^{p3t+p4} + p5e^{p6t+p7}
    """
    return p1 + p2*np.exp(p3*tof) + p5*np.exp(p6*tof)

def e3p0(tof,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10):
    """
    Background function for TOF spectra

    Parameters
    ----------
    tof : array-like
        The time-of-flight spectrum
    p1 : float
        constant background
    p2 : float
        multiplier on 1st exponential
    p3 : float
        multiplier on time-of-flight in 1st exponent
    p4 : float
        constant added to 1st exponent
    p5-p10 : float
        (see equation in notes)

    Returns
    -------
    e3p0 : array-like
        The function in the length of t (see notes)

    Notes
    -----
    .. math:: f(t) = p1 + p2e^{p3t+p4} + p5e^{p6t+p7} + p8e^{p9t+p10}
    """
    return p1 + p2*np.exp(p3*tof+p4) + p5*np.exp(p6*tof+p7) + p8*np.exp(p9*tof+p10)



