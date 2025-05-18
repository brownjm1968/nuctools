import numpy as np
import pandas as pd
from . import tof_tools as tt
from . import physics_tools as pt

__all__ = ['ags','e1p0','e2p0','e3p0','rebin_tof_shift','rebin_tof_shift_slow']

class ags:
    """
    Python class to contain AGS-like functions and information. Please see the AGS Manual 
    (B. Becker, C. Bastian, J. Heyse, S. Kopecky, P. Schillebeeckx, AGS - Analysis of Geel
    Spectra, European Commission, Joint Research Centre, September 2014) for more detail: 
    https://www.oecd-nea.org/jcms/pl_19568/ags-analysis-of-geel-spectra-users-manual

    Attributes
    ----------
    tof : array-like
        A 1-d array that is a Pandas Series of time-of-flight values,
        typically given in [us]
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
        self.dtof = None
        self.counts = None
        self.cps = None
        self.energy = None
        self.obs = None
        self.unc_obs = None
        
    def read_grouped_counts(self,spectrum,comp_pt,comp_fct,binsize,triggers):
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
            The width of the base bin in [us]
        triggers : integers
            the number of times the linac fired

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
        #hist = pd.read_csv(filename,sep=r'\s+',names=['bin','counts'])
        hist = pd.DataFrame({
            'counts'    : spectrum,
            'dtof' : np.zeros(len(spectrum))
            })
        hist['cps'],hist['dcps'],hist['tof'] = 0,0,0
        
        for i in range(len(cp)-1):
            hist.loc[cp[i]:cp[i+1]-1,'cps'] = hist.counts[cp[i]:cp[i+1]]/(cf[i]*triggers*binsize*1e-6)
            hist.loc[cp[i]:cp[i+1]-1,'dcps'] = np.sqrt(hist.counts[cp[i]:cp[i+1]])/(cf[i]*triggers*binsize*1e-6)
            hist.loc[cp[i]:cp[i+1]-1,'dtof'] = cf[i]*binsize
            if i==0:
                hist.loc[cp[i]:cp[i+1]-1,'tof'] = 0+np.arange(cp[i+1]-cp[i])*binsize*cf[i]
            else:
                hist.loc[cp[i]:cp[i+1]-1,'tof'] = hist.tof[cp[i]-1]+np.arange(cp[i+1]-cp[i])*binsize*cf[i]
                
        self.tof = hist.tof
        self.dtof = hist.dtof
        self.cps = hist.cps
        self.dcps = hist.dcps
        self.counts = hist.counts

    def calc_energy(self,FP,t0):
        """
        Calculate energy from the existing TOF spectrum

        Parameters
        ----------
        FP : float
            The flight path length the neutron traveled [m]
        t0 : float
            The time-zero used to correct the TOF spectrum [us]
        """
        self.energy = tt.tofe(self.tof-t0,FP)

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


def dtco(dtof,counts,dead_time,trigs):

    lc = len(counts)
    SUM = np.zeros(lc)
    time_window = 0.0
    for i in range(lc):
        if tof[i] < dead_time:
            SUM[i] = 0.0
        else: 
            i0 = None
            time_sum = 0
            tind = i-1
            while time_sum < dead_time:
                time_sum += dtof[tind]
                tind-1
            SUM = 1.
    corr_counts = trigs * ( -np.log( 1-(counts/trigs)/(1-SUM) ) )


def rebin_tof_shift_slow(spectrum,newtof,fp1,fp2,gf1,gf2):
    """
    Rebin counts from `spectrum` into the `new_tof` grid
    based on differences in flight paths 1 and 2 and gflash
    locations 1 and 2. Based on "old" AGS shift function. 

    Parameters
    ----------
    spectrum : Pandas DataFrame
        Must have columns of: "tof","counts","dcounts"
    newtof : array
        1 dim array of tof, ascending order, units [us]
    fp1 : float
        Original position, units [m]
    fp2 : float
        New position, units [m]
    gf1 : float
        Original gamma flash location, units [us]
    gf2 : float
        New gamma flash location, units [us]

    Returns
    -------
    new_spect : Pandas DataFrame
        DataFrame with 3 columns: tof, counts, error counts in new binning
    """

    c = pt.c # [m/us]

    orig_binwidth = spectrum.tof.diff().shift(-1)
    ls = len(spectrum)
    newspec = pd.DataFrame()
    newspec['counts'] = np.zeros(ls)
    newspec['dcounts'] = np.zeros(ls)
    newspec['tof'] = newtof
    newspec['binwidth'] = newspec.tof.diff().shift(-1)

    new_chan = 0
    lim_chan = 0
    orig_chan = 0
    last_bin = ls
    while orig_chan < last_bin:

        oll = spectrum.tof[orig_chan]        # original tof lower limit
        oul = oll + orig_binwidth[orig_chan] # original tof upper limit
        nll = newspec.tof[lim_chan]
        nul = nll + newspec.binwidth[lim_chan]
        oc = spectrum.loc[orig_chan,'counts']
        doc = spectrum.loc[orig_chan,'dcounts']

        # transform orig limits and width to new FP
        delta_ll = (oll-gf1+fp1/c)*fp2/fp1 - (oll-gf2+fp2/c)
        delta_ul = (oul-gf1+fp1/c)*fp2/fp1 - (oul-gf2+fp2/c)
        orig_width = oul+delta_ul - (oll+delta_ll) 

        new_chan = lim_chan
        while lim_chan < last_bin:
            nll = newspec.loc[lim_chan,'tof']
            nul = nll + newspec.loc[lim_chan,'binwidth']
            lim_chan += 1
            if oll+delta_ll >= nul:
                new_chan += 1
            else:
                break

        if (new_chan==0) and (oll+delta_ll < nll) and (oul+delta_ul > nll):
            weight = (oul+delta_ul - nll)/orig_width
            newspec.loc[new_chan,'counts']  += oc*weight
            newspec.loc[new_chan,'dcounts'] += doc*weight

        # -------- old bin straddles left side of new bin
        if (oll+delta_ll > nll):
            # -------- all of old bin within new bin 
            if(oul+delta_ul <= nul):
                weight = 1
                newspec.loc[new_chan,'counts']  += oc*weight
                newspec.loc[new_chan,'dcounts'] += doc*weight
            # -------- old bin straddles right side of new bin
            else:
                weight = (nul-(oll+delta_ll))/orig_width
                newspec.loc[new_chan,'counts']  += oc*weight
                newspec.loc[new_chan,'dcounts'] += doc*weight
                # -------- keep incrementing new channels until old upper limit above new upper limit
                while lim_chan < last_bin:
                    nll = newspec.tof[lim_chan]
                    nul = nll + newspec.binwidth[lim_chan]
                    lim_chan += 1
                    if (new_chan >= last_bin) or (oul+delta_ul <= nul):
                        break
                    new_chan += 1
                    new_width = nul-nll
                    weight = new_width/orig_width
                    newspec.loc[new_chan,'counts']  += oc*weight
                    newspec.loc[new_chan,'dcounts'] += doc*weight
                # -------- after all that we landed straddling lower limit of the new bin
                if (new_chan < last_bin-1) and ((oul+delta_ul) > nll):
                    new_chan += 1
                    weight = (oul+delta_ul-nll)/orig_width
                    newspec.loc[new_chan,'counts']  += oc*weight
                    newspec.loc[new_chan,'dcounts'] += doc*weight
                    
        orig_chan += 1
        lim_chan = new_chan - 1
        if lim_chan < 0: 
            lim_chan = 0

    return newspec


def rebin_tof_shift(spectrum, newtof, fp1, fp2, gf1, gf2):
    """
    Rebin counts from `spectrum` into the `new_tof` grid
    based on differences in flight paths 1 and 2 and gflash
    locations 1 and 2. Based on "old" AGS shift function. 

    Parameters
    ----------
    spectrum : Pandas DataFrame
        Must have columns of: "tof","counts","dcounts"
    newtof : array
        1 dim array of tof, ascending order, units [us]
    fp1 : float
        Original position, units [m]
    fp2 : float
        New position, units [m]
    gf1 : float
        Original gamma flash location, units [us]
    gf2 : float
        New gamma flash location, units [us]

    Returns
    -------
    new_spect : Pandas DataFrame
        DataFrame with 3 columns: tof, counts, error counts in new binning
    """

    c = pt.c # [m/us]

    # Original TOF edges
    orig_tof = spectrum['tof'].values
    orig_counts = spectrum['counts'].values
    orig_dcounts = spectrum['dcounts'].values
    orig_binwidth = np.diff(orig_tof, append=orig_tof[-1])  # handle last bin gracefully
    oul = orig_tof + orig_binwidth  # upper limits

    # New TOF edges
    newtof = np.asarray(newtof)
    new_binwidth = np.diff(newtof, append=newtof[-1])
    new_upper = newtof + new_binwidth

    # Initialize new spectrum arrays
    counts_out = np.zeros_like(newtof)
    dcounts_out = np.zeros_like(newtof)

    # Shifted TOF bounds
    delta_ll = (orig_tof - gf1 + fp1/c) * fp2/fp1 - (orig_tof - gf2 + fp2/c)
    delta_ul = (oul - gf1 + fp1/c) * fp2/fp1 - (oul - gf2 + fp2/c)
    shifted_ll = orig_tof + delta_ll
    shifted_ul = oul + delta_ul
    orig_width = shifted_ul - shifted_ll

    # Loop through each original bin and distribute into new bins
    for i in range(len(orig_tof)):
        sl = shifted_ll[i]
        su = shifted_ul[i]
        w = orig_width[i]
        if w <= 0:
            continue  # skip bad bins

        overlaps = (newtof < su) & (new_upper > sl)
        for j in np.where(overlaps)[0]:
            overlap_start = max(sl, newtof[j])
            overlap_end = min(su, new_upper[j])
            frac = (overlap_end - overlap_start) / w
            counts_out[j] += frac * orig_counts[i]
            dcounts_out[j] += frac * orig_dcounts[i]

    return pd.DataFrame({
        'tof': newtof,
        'counts': counts_out,
        'dcounts': dcounts_out,
    })



