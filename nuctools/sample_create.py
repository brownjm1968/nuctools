import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
from . import tof_tools
import os,sys

__all__ = ["sample",'get_montrig','count_rate','dtc','inc']

class sample:
    """
    Store all relevant information for a sample being 
    studied at the RPI Linac. 

    Stores and calculates data about a sample studied at
    the RPI Linac and follows the file structure and naming 
    conventions that are found in codes originating from 
    Brian McDermott, such as "40mprocess" and "h5query".
    The output from these codes is used as input for this
    class.

    Parameters
    ----------
    samp_name : str
        Name of sample of interest, e.g. Ta-2mm, MUST be 
        the same as it appears in file structure of initial
        raw binary files as it is typically propagated 
        through to the HDF5 file.
    data_type : str
        Type of data which was collected, e.g. yield,
        or transmission
    tof_arr : numpy array
        Time of flight data, 1st column is ascending
        time bins [us], 2nd column is bin width, 3rd
        column is the summed counts
    hdf_address : str, optional
        File address of the hdf file output from 40mprocess
        which contains: monitors, triggers, and pulse events;
        must end with slash
    phys_dims : array_like, optional
        Physical dimensions [cm] of sample being studied, must be 
        in format [thickness,length,width] or [thickness,radius]
        where thickness is the length of sample parallel to the
        Linac beam
    tau : float, optional
        The dead time of the system [us] used to collect the counts
        data, default value is 0 
    paralyzeable : bool, optional
        If system used to collect counts data was paralyzeable 
        a different dead time model must be used, default is 
        non-paralyzeable dead-time model
    plot_montrig : bool, optional
        Plot the monitor ratio through the montrig function
    factor : 1-d numpy array, optional
        Compression factors for the tof groups, 1st compression factor
        corresponds to time of flights below c_pt[0], 2nd corresponds to tofs above
        c_pt[0] and below c_pt[1], and so on until end of range is found.
    c_pts : 1-d numpy array, optional
        Compression time of flights that seperate the bins that are to be grouped.
        List points from lowest to highest.
    avg : bool, optional
        Return an average of group values instead of
        the sum in the group. group error is 1/N*(sum(error^2))^(1/2)
    binedge : bool, optional
        Return TOF column for grouped data as the first
        tof in group instead of the mean tof of the
        group.
    verbose : bool, optional
        Choose whether to print out approximate comp pt energy locations
        and indices. Default is False.
    FP : float, optional
        Flight path needed to convert the time of flight values of bins to
        energies for the compression point selection.
    verbose : bool, optional
        Option to print out extra information

    Notes
    -----
    Will eventually require break points and grouping 
    metadata.

    Allows for yield calculations, must later include 
    transmission.
    """
    def __init__(self,samp_name,data_type,tof_arr,hdf_address=None,
                 phys_dims=None,tau=0.0,paralyzeable=False,
                 plot_montrig=False,factor=None,c_pts=None,
                 avg=False,binedge=False,FP=None,verbose=False):
        
        if verbose == True:
            print("--------\n {} data for {}\n--------".format(data_type,
                  samp_name))

        self.samp_name = samp_name
        self.phys_dims = phys_dims
        self.data_type = data_type

        #------------------------------------------------------------
        #------------ Monitors and triggers -------------------------
        #------------------------------------------------------------
        default_montrig = False

        if hdf_address == None:

            default_montrig = True
            print("----------\n----------")
            print("Monitors and triggers not given, cps calculated with")
            print("default 1e8 monitor and trigger counts.")
            print("----------\n----------")

            self.montrig = np.array([1e8,1e8])

        else:
            self.montrig = get_montrig(hdf_address,samp_name,plot=plot_montrig,
                                       verbose=verbose)

        #------------------------------------------------------------
        #------------ Count rate ------------------------------------
        #------------------------------------------------------------
        if default_montrig:
            self.cps = count_rate(tof_arr,self.montrig[1],tau,paralyzeable=paralyzeable,
                 factor=factor,c_pts=c_pts,avg=avg,binedge=binedge,verbose=verbose,FP=FP)
        else: 
            self.cps = count_rate(tof_arr,self.montrig[1][1],tau,paralyzeable=paralyzeable,
                 factor=factor,c_pts=c_pts,avg=avg,binedge=binedge,verbose=verbose,FP=FP)
        #------------------------------------------------------------


#------------------------------------------------------------------------------

def get_montrig(hdf_address,samp_name,plot=False,verbose=False):
    """
    Collect monitors and triggers used to reduce
    capture data to proper dead time corrected
    monitor normalized count rates

    Parameters
    ----------
    hdf_address : str
        File address of the HDF5 file, must end with "/"
    samp_name : str
        Name of the sample being analyzed. Should
        match the sample name in the HDF5 filename
    plot : bool, optional
        Choose whether to plot the monitors and triggers
        as a function of cycle number
    verbose : bool, optional
        Choose whether to print out trigger and monitor 
        information.

    Returns
    -------
    montrig : numpy array
        Vector of monitor counts and triggers

    Notes
    -----
    The montrig array is printed out when function
    is called, user can see order of monitors and
    triggers.
    """

    with h5.File(hdf_address+samp_name+"_master.h5","r") as hdf:

        mon_counts = []
        det_counts = []
        sum_mon = np.zeros(8)
        sum_trig = np.zeros(2)

        for cycle in hdf[samp_name]:
            mon_vector = hdf[samp_name+"/"+cycle+"/MON"]
            trig_vector = hdf[samp_name+"/"+cycle+"/TRIG"][0]
            datag_array = hdf[samp_name+"/"+cycle+"/DATA/GOOD"]
            datat_array = hdf[samp_name+"/"+cycle+"/DATA/TRUNC"]
            sum_detcts = len(datag_array)+len(datat_array)
            mon_counts.append(mon_vector)
            det_counts.append(sum_detcts)

            # number of monitors stored in /MON != to 8
            if (len(sum_mon) != len(mon_vector)) or (len(sum_trig) != len(trig_vector)):
                raise ValueError("Unexpected number of mon's in /MON")

            sum_mon += mon_vector
            sum_trig += trig_vector

        num_cyc = len(hdf[samp_name])
        mon_counts = np.transpose(mon_counts)
        # monitor ratios 
        mr0,mr1,mr2 = np.zeros(num_cyc),np.zeros(num_cyc),np.zeros(num_cyc)
        for i in range(num_cyc):
            mr0[i] = det_counts[i]/mon_counts[0][i]
            mr1[i] = det_counts[i]/mon_counts[1][i]
            mr2[i] = det_counts[i]/mon_counts[2][i]
        # std of monitor cycles / mean of mon. cycles
        stdmean = [0,0,0]
        stdmean[0] = np.std(mr0)/np.mean(mr0)
        stdmean[1] = np.std(mr1)/np.mean(mr1)
        stdmean[2] = np.std(mr2)/np.mean(mr2)
        
        if plot:
            x = range(num_cyc)
            plt.plot(x,mr0/np.mean(mr0),label="Sum/Mon 0")
            plt.plot(x,mr1/np.mean(mr1),label="Sum/Mon 1")
            plt.plot(x,mr2/np.mean(mr2),label="Sum/Mon 2")
            plt.legend()
            plt.show()
        
        montrig = np.concatenate((sum_mon,sum_trig))
        if verbose == True:
            print("         Counts          Std/mean")
            for i in range(len(stdmean)):
                print("Mon. {} : {}  ,  {}".format(i,sum_mon[i],stdmean[i]))
            if len(sum_trig>1):
                print("Triggers : {}  {}".format(montrig[6],montrig[7]))

    return montrig


def count_rate(tof_arr,trig,tau=None,paralyzeable=False,factor=None,c_pts=None,
               avg=False,binedge=False,verbose=False,FP=None):
    """
    Converts counts to dead-time corrected count
    rates

    Parameters
    ----------
    tof_arr : 2-d numpy array
        Time of flight data, 1st column is ascending
        time bins [us], 2nd column is the summed counts,
        and the 3rd column (if desired) is error on the 
        counts. 3rd column only necessary if avg==True.
        !! time of flight data **MUST** be equal bin sizes !!
    trig : float
        Triggers recorded for the sample
    tau : float
        The dead time of the system used to collect 
        counts data [us]
    paralyzeable : bool,optional
        If the system used to collect the counts data
        was paralyzeable, a different dead time model
        must be used
    factor : 1-d numpy array, optional
        Compression factors for the tof groups, 1st compression factor
        corresponds to time of flights below c_pt[0], 2nd corresponds to tofs above
        c_pt[0] and below c_pt[1], and so on until end of range is found.
    c_pts : 1-d numpy array
        Compression time of flights that seperate the bins that are to be grouped.
        List points from lowest to highest.
    avg : bool, optional
        Return an average of group values instead of
        the sum in the group. group error is 1/N*(sum(error^2))^(1/2)
    binedge : bool, optional
        Return TOF column for grouped data as the first
        tof in group instead of the mean tof of the
        group.
    verbose : bool, optional
        Choose whether to print out approximate comp pt energy locations
        and indices. Default is False.
    FP : float
        Flight path needed to convert the time of flight values of bins to
        energies for the compression point selection.

    Returns
    -------
    dc_arr : 2-d numpy array
        TOF data, 1st column ascending time bins [us],
        2nd column is count rate [cps], 3rd column is 
        the error on that count rate

    Notes
    -----
    Units matter!!

    Error for count rates are sqrt of the counts vector user 
    provides.

    Dead time correction error given in the dead time correction
    function `dtc()`.

    .. math:: CPS = \\frac{C}{b_w Trig}

    .. math:: \\Delta C = \\sqrt{C}

    .. math:: \\Delta CPS = \\sqrt{\\left(\\frac{\\Delta C}{b_w Trig}\\right)^2}

    .. math:: CPS_{DC} = (CPS) dtcf 

    .. math:: \\Delta CPS_{DC}=\\sqrt{\\left(\\frac{\\Delta CPS}{CPS} \\right)^2 + \\left(\\frac{\\Delta dtcf}{dtcf} \\right)^2}

    """
    use_single_group = False
    use_comp_group = False

    # -------------------------------------------------------------------------
    # ---- run some checks ----------------------------------------------------
    # -------------------------------------------------------------------------
    # did the user give one grouping input but not the other?
    mismatch1 = ((factor is None) and (c_pts is not None))
    mismatch2 = ((factor is not None) and (c_pts is None))
    if mismatch1:
        raise ValueError("c_pts given, but factor is not.")
    # if mismatch2, we may need to use single group if factor is not an array
    if mismatch2:
        try:
            len(factor)
            raise ValueError()
        except TypeError:
            print('Using single group')
            use_single_group = True
        except:
            print('Factor is an array and c_pts is not given. Aborting.')
            raise ValueError("Check input for argument: factor.")
    if factor is not None and c_pts is not None:
        use_comp_group = True
    if use_comp_group:
        # we need the flight path for comp_group
        if FP is None:
            raise ValueError("Specify the flight path for comp_group()")

    # cast array into np array
    tof_arr = np.array(tof_arr)
    # if input is pandas, may need to transpose
    if tof_arr.shape[0] > tof_arr.shape[1]:
        tof_arr = np.transpose(tof_arr)

    # -------------------------------------------------------------------------
    # ---- Calculate and apply dead time correction ---------------------------
    # -------------------------------------------------------------------------
    if tau is None:
        tau = 0.0     # will be "corrected" by 1.0, with error 0.0

    tau = tau*1e-6    # [us --> s]

    if paralyzeable == False:
        bin_width = tof_arr[0][1]-tof_arr[0][0]
        dtcf = dtc(tof_arr[1],tau,trig,bin_width)
        dc_counts = tof_arr[1]*dtcf[0]
        # sum of the squares of relative errors in counts and dead-time corr.
        # time the corrected counts
        if avg:
            ecounts = tof_arr[2]
        else:
            ecounts = np.sqrt(tof_arr[1])

        rel_err_dtcf = (dtcf[1]/dtcf[0])
        rel_err_counts = (ecounts/tof_arr[1])
        edc_counts = np.sqrt(rel_err_dtcf**2 + rel_err_counts**2) * dc_counts
        # --------------------------------------------
        # check if error on correction is significant
        # --------------------------------------------
        ratio = rel_err_dtcf/(rel_err_dtcf+rel_err_counts)
        if (ratio > 0.02).any():
            print("Max % of dead time error: ",max(ratio))
            print("!! Dead time correction factor accounts for > 2% error in the counts !!")
            print("!! This is not accounted for in grouping if avg is not used !!")
            print("To fix, you have to change the hard-coded grouping function.")

    else:
        raise NotImplementedError("paralyzeable dead time not implemented.")

    # -------------------------------------------------------------------------
    # ---- Group the counts if user has specified -----------------------------
    # -------------------------------------------------------------------------
    if use_single_group:
        # include the input error on the counts
        gtof_array = tof_tools.single_group([tof_arr[0],dc_counts,ecounts],
            factor,avg=avg,binedge=binedge)
    if use_comp_group:
        # include the input error on the counts
        gtof_array = tof_tools.comp_group([tof_arr[0],dc_counts,ecounts],
            factor,c_pts,FP,avg=avg,binedge=binedge,verbose=verbose)

    # -------------------------------------------------------------------------
    # ---- Calculate the count rate -------------------------------------------
    # -------------------------------------------------------------------------
    if use_comp_group is False and use_single_group is False:
        cps = dc_counts/bin_width/trig/1e-6
        # error in the cps 
        err_cps = 1/bin_width/trig/1e-6 * edc_counts
        # stage a tof array for the return statement
        time_of_flight = tof_arr[0]
    else:
        cps = gtof_array[2]/gtof_array[1]/trig/1e-6
        # error in cps
        err_cps = 1/gtof_array[1]/trig/1e-6 * gtof_array[3]
        # stage a tof array for the return statement
        time_of_flight = gtof_array[0]

    cps_arr = np.array([time_of_flight,cps,err_cps])


    return cps_arr



def dtc(counts,dead_time,trigs,bin_width):
    """
    Calculate the TOF dead-time correction factor according
    to the non-paralyzable model

    Parameters
    ----------
    counts : array_like
        A 1-d array of counts as a function of tof (tof 
        in ascending order, equal bin spacing throughout)
    dead_time : float
        The dead time [us] associated with the data collection 
        system.
    trigs : int
        The number of triggers = the number of times channel 
        was collecting data. For RPI the number of LINAC 
        pulses.
    bin_width : float
        The width of the bins in time [us]. Should be constant
        for the array of counts.

    Returns
    -------
    dtcf : 2-d numpy array
        The dead time correction factor for the counts array
        based on the non-paralyzeable dead time model is in 
        column 0 of array dtcf, the error associated with it
        is in column 1. 
        len(dtcf[0]) = len(dtcf[1]) = len(counts).

    Notes
    -----
    dtc() returns 0's for the first bins that span the dead 
    time. 1's would keep the signal closer to reality, but 
    0's make it obvious to the user that those channels do 
    not have a correction applied, and data in this time of
    flight range is rarely used since it usually occurs before
    the gamma flash.

    .. math:: dtcf_{\\tau>b_w} = \\frac{1}{1-\\frac{SUM}{trig}} = \\frac{1}{1-\\frac{\\sum_iw_iC_i}{trig}}

    where the weights are all 1 except the first and last bin. 
    First weight is 0.5, last weight is the fraction of bins
    the dead time minus 0.5 covers (e.g. 3.1 minus 1/2 = 2.6) 
    minus the whole bins dead time minus 0.5 covers (2), so 
    the last weight is 0.6.

    Now the error.

    .. math:: \\Delta SUM = \\sqrt{\\sum_i(w_i\\sqrt{C_i})^2}

    .. math:: \\Delta dtcf_{\\tau>b_w} = \\sqrt{\\left(\\frac{\\partial dtcf_{\\tau>b_w}}{\\partial SUM}\\Delta SUM\\right)^2} = \\sqrt{\\left(\\frac{\\Delta SUM}{trig\\left(1-\\frac{SUM}{trig}\\right)^2}\\right)^2} 

    .. math:: = \\frac{\\Delta SUM}{trig\\left(1-\\frac{SUM}{trig}\\right)^2} = \\frac{\\Delta SUM}{trig}dtcf_{\\tau>b_w}^2

    If the dead time is less than the bin width:

    .. math:: dtcf_{\\tau<b_w} = \\frac{1}{1-\\frac{C\\tau}{b_wtrig}}

    .. math:: \\Delta dtcf_{\\tau<b_w} = \\sqrt{\\left(\\frac{\\partial dtcf_{\\tau<b_w}}{\\partial C}\\Delta C\\right)^2} = \\frac{\\tau\\sqrt{C}}{b_wtrig\\left(1-\\frac{C\\tau}{b_wtrig}\\right)^2}
    
    Algorithm is a Pythonized version of the dead-time-correction algorithm by Y. Danon, "Design and Construction of the RPI Enhanced Thermal Neutron Target and Thermal Cross Section Measurements of Rare Earth Isotopes.", Doctoral Thesis, (1993).
    """

    # cast as numpy = major speedup
    counts = np.array(counts)
    lc = len(counts)
    dtcf = np.zeros(lc)
    ddtcf = np.zeros(lc)
    if bin_width>=dead_time:
        x = counts*dead_time/bin_width/trigs
        dtcf = 1/(1-x)
        ddtcf = dead_time*np.sqrt(counts)/bin_width/trigs/(1-x)**2
    else:
        # how many bins does the dead time span - 0.5 a bin
        an = dead_time/bin_width-0.5
        # how many whole bins does the dead time span
        n = int(an)
        # bins 0-n will be left 0 
        print("Channels 0 -",n,"set to 0")
        # fraction of earliest bin that dead time spans
        f1 = an-n
        # add half of current bin
        SUM = counts[n+1:lc] * 0.5
        # (w_i(C_i)^0.5)^2 = C_i*0.5*0.5
        sSUM = SUM*0.5
        # add preceding whole bins within span of dead time
        for j in range(1,n+1):
            SUM += counts[n+1-j:lc-j]
            sSUM += counts[n+1-j:lc-j]
        # add the fraction of the last bin dead time spans
        SUM += f1*counts[0:lc-n-1]
        sSUM += f1**2*counts[0:lc-n-1]
        dSUM = np.sqrt(sSUM)

        dtcf[n+1:lc] = 1/(1-SUM/trigs)
        ddtcf[n+1:lc] = dSUM/trigs*dtcf[n+1:lc]**2
    return np.array([dtcf,ddtcf])


def inc(x):
    return x + 1





