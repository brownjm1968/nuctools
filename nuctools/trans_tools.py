import numpy as np
import pandas as pd
from . import tof_tools as tt 
from . import math_tools as mt
from . import sam_tools as st

__all__ = ['Trans','calc_cov','write_covar_file','read_covar_file']

# TODO: create sample/open objects

class Trans:
    """
    This is a class designed to hold all the data and methods necessary to calculate
    neutron transmission from a time-of-flight measurement at RPI.
    
    Parameters
    ----------
    
    open_data_name : string
        Full path name for the file containing the OPEN data
    open_mon_name : string
        Full path name for the file containing the OPEN monitors
    sample_data_name : string
        Full path name for the file containing the SAMPLE data
    sample_mon_name : string
        Full path name for the file containing the SAMPLE data
    datafile_names : list, optional
        List of strings for the column names for the data. OPEN and
        SAMPLE DATA should typically have the same structure
    tau : float, optional
        Tau is the dead time in micro seconds
    correct_deadtime : bool, optional
        Whether the deadtime of the system needs to be corrected for 
        in the counts data
    mon_index : int, optional
        Which monitor to use. Historically it is 0-10, with the Linac 
        triggers being 0,1
    bkgfunction : string, optional
        Options are: \'power_law\', \'double_exp\', \'double_power_law\',
        and \'exp_log\'.
        There are only four background shapes, the coefficients are 
        controlled with o_bkgcoeff, and s_bkgcoeff parameters
    o_bkgcoeff : list, optional
        Coefficients for the either of the two background fucntions 
        programmed for the OPEN data. power_law needs 3 coeff's, 
        double_exp needs 5 coeffecients, double_power_law requires
        5 coefficients, exp_log requires 4 coeff's
    s_bkgcoeff : list, optional
        Coefficients for the either of the two background fucntions 
        programmed for the SAMPLE data. power_law needs 3 coeff's, 
        double_exp needs 5 coeffecients, double_power_law requires
        5 coefficients, exp_log requires 4 coeff's
    factor : 1-d numpy array, optional
        Compression factors for the tof groups, 1st compression factor
        corresponds to time of flights below c_pt[0], 2nd corresponds to tofs above
        c_pt[0] and below c_pt[1], and so on until end of range is found.
    c_pts : 1-d numpy array, optional
        Compression time of flights that seperate the bins that are to be grouped.
        List points from lowest to highest (should re-write to sort and remove user
        capability to fail.)
    FP : float, optional
        Flight path needed to convert the time of flight values of bins to
        energies for the compression point selection.
    samp_bkg_err : float, optional
        The fractional error for every point on the background. Default is
        0.2 (20%). This is a temporary way to conservatively estimate error 
        on the transmission.
    open_bkg_err : float, optional
        The fractional error for every point on the background. Default is
        0.2 (20%). This is a temporary way to conservatively estimate error 
        on the transmission.
    optstat_group : bool, optional
        Whether or not you want to use an optimized grouper 
    D : float, optional
    res_per_bin : int, optional
    E_cpts : array-like
    
    Attributes
    ----------
    
    open_data : DataFrame
        DataFrame of time-of-flight (tof), bin widths (bin_width), counts 
        (counts), count rate (cps) and error on the count rate (dcps) from
        method calc_count_rate, and background (bkg) for the OPEN. All 
        time values should be in units of microseconds.
    sample_data : DataFrame
        DataFrame of time-of-flight (tof), bin widths (bin_width), counts 
        (counts), count rate (cps) and error on the count rate (dcps) from
        method calc_count_rate, and background (bkg) for the SAMPLE. All 
        time values should be in units of microseconds.
        All time values should be in units of microseconds. Also contains
        the count rate and error in the count rate as given by the method
        calc_count_rate
    open_mon : DataFrame
        Monitor detector counts for OPEN. 1-d list of integers
    sample_mon : DataFrame
        Monitor detector counts for SAMPLE. 1-d list of integers
    open_trig : int
        The triggers used to calculate the count rate for the OPEN. 
        Typically this should be the number of Linac triggers, given at
        mon_index=0.
    sample_trig : int
        The triggers used to calculate the count rate for the SAMPLE. 
        Typically this should be the number of Linac triggers, given at
        mon_index=0.
    
    Methods
    -------
    
    calc_count_rate
        Calculates the count rate for the sample given, such as OPEN or
        SAMPLE
    calc_open_rate
        Sets the 'cps' and 'dcps' columns in DataFrame open_data with the
        calculation from calc_count_rate
    calc_sample_rate
        Sets the 'cps' and 'dcps' columns in DataFrame sample_data with the
        calculation from calc_count_rate
    power_law 
        Uses coefficients given to calculate the power law function at the 
        time-of-flight values given
    double_exp 
        Uses coefficients given to calculate the double_exp function at the 
        time-of-flight values given
    calc_trans
        Sets the attribute 
    calc_sigma
        Calculates the cross section from transmission and atoms/barn thickness
        
    Notes
    -----
    Grouping for the transmission goes as follows:

    .. math:: C_g = \\sum_i^NC_i, \\qquad b_w = (N)(b_{w,raw})

    .. math:: \\dot C_{g} = \\frac{C_g}{b_w trig}

    .. math:: T_g = \\frac{\\dot C_{g,s}}{\\dot C_{g,O}} = \\frac{\\frac{\\sum_i^NC_s}{trig_sb_w}}{\\frac{\\sum_i^NC_O}{trig_Ob_w}}

    Background has four functions for now:
    
    `power_law`

    .. math:: a(tof)^{-b}+c

    `double_exp`

    .. math:: ae^{-b(tof)}+ce^{-d(tof)}+e

    `double_power_law`

    `exp_log`

    .. math:: f(t) = e^{ a + (b)t + c/ln(t) }

    Transmission from background subtracted count rates:

    .. math:: T_g = \\frac{\\dot C_{g,s}-\\dot B_{s}}{\\dot C_{g,O}-\\dot B_O}

    The error on transmission is propagated from the error on the counts (square root) and uses an 
    input constant fraction error on the background.

    .. math:: \\Delta T = \\sqrt{\\left(\\frac{\\partial T}{\\partial x}\\Delta x\\right)^2 + ...} 

    .. math:: = \\sqrt{\\left(\\frac{\\Delta \\dot C_s}{\\dot C_O - \\dot B_O}\\right)^2 + \\left(\\frac{\\Delta \\dot B_s}{\\dot C_O - \\dot B_O}\\right)^2 +\\left(\\frac{(\\dot C_s - \\dot B_s)\\Delta \\dot C_O}{\\dot C_O - \\dot B_O}\\right)^2 +\\left(\\frac{(\\dot C_s - \\dot B_s)\\Delta \\dot B_O}{\\dot C_O - \\dot B_O}\\right)^2}

    If pandas throws an "unsupported type" error, check your input file that is being
    read. It's possible you gave a file to be grouped, but you didn't tell the program
    to group.
    
    """
    def __init__(self,open_data_name,open_mon_name,sample_data_name,sample_mon_name,
                 datafile_names=None,tau=None,correct_deadtime=False,mon_index=0,
                 bkgfunction='power_law',o_bkgcoeff=None,s_bkgcoeff=None,factor=None,
                 c_pts=None,FP=None,samp_bkg_err=0.2,open_bkg_err=0.2,optstat_group=False,
                 D=None,res_per_bin=None,E_cpts=None,verbosity=False,binedge=False):
    
    
        self.mon_index = mon_index
        self.correct_deadtime = correct_deadtime
        self.bkgfunction = bkgfunction
        self.samp_bkg_err = samp_bkg_err
        self.open_bkg_err = open_bkg_err

        # ---------------------------------------------------------------------------------
        # ----- decide how to choose some attributes and defaults -------------------------
        # ---------------------------------------------------------------------------------
        if o_bkgcoeff is not None:
            self.o_bkgcoeff = o_bkgcoeff
        if s_bkgcoeff is not None:
            self.s_bkgcoeff = s_bkgcoeff
        
        if bkgfunction == 'power_law':
            if o_bkgcoeff is None:
                self.o_bkgcoeff = [8e5,1.4,102]
            if s_bkgcoeff is None:
                self.s_bkgcoeff = [8e5,1.45,83]
        elif bkgfunction == 'double_exp':
            if o_bkgcoeff is None:
                self.o_bkgcoeff = [4e5,2e-1,4e3,8e-3,120]
            if s_bkgcoeff is None:
                self.s_bkgcoeff = [4e5,2e-1,4e3,8e-3,120]
        elif bkgfunction == 'double_power_law':
            if o_bkgcoeff is None or s_bkgcoeff is None:
                raise ValueError("User must specify double_power_law bkg coefficients.")
        elif bkgfunction == 'exp_power_mix':
            if o_bkgcoeff is None or s_bkgcoeff is None:
                raise ValueError("User must specify exp_power_mix bkg coefficients.")
        elif bkgfunction == 'exp_log':
            if o_bkgcoeff is None or s_bkgcoeff is None:
                raise ValueError("User must specify exp_log bkg coefficients.")

        # the names for the Ogrp function from rpixdr (tof,bin_width,counts)
        if datafile_names is None:
            self.datafile_names = ['tof','bin_width','counts']
        else:
            self.datafile_names = datafile_names

        if E_cpts is None:
            E_cpts = [1e2,1e3,1e4]
        
        # ---------------------------------------------------------------------------------
        # ----- read in the counts, triggers, and monitors --------------------------------
        # ---------------------------------------------------------------------------------
        
        grouping_is_used = False
        # read in the counts data
        if factor is None and c_pts is None and optstat_group==False:
            # -----------------------------------------------------------------
            # ---- read in optgroup files that have already been grouped ------
            # -----------------------------------------------------------------
            self.open_data   = pd.read_csv(open_data_name  ,sep=r'\s+',names=self.datafile_names)
            self.sample_data = pd.read_csv(sample_data_name,sep=r'\s+',names=self.datafile_names)
        elif optstat_group == True:
            # -----------------------------------------------------------------
            # ---- group by "optimum statistics" ------------------------------
            # -----------------------------------------------------------------
            if D==None or res_per_bin==None or FP==None:
                print("Fix input parameters.")
                raise ValueError("Missing: level spacing (D) or res. desired per bin (res_per_bin) or FP.")
            # ---- reset the column names -------------------------------------------------
            self.datafile_names = ['ch','counts','dcounts']
            # ---- read in the fixed bin width data ---------------------------------------
            ug_open_data = pd.read_csv(open_data_name,sep=r'\s+',names=self.datafile_names,
                                         skiprows=14)
            ug_samp_data = pd.read_csv(sample_data_name,sep=r'\s+',names=self.datafile_names,
                                         skiprows=14)
            # ---- read in the bin width --------------------------------------------------
            open_bw = np.genfromtxt(open_data_name  ,unpack=True,skip_header=12,max_rows=1)
            samp_bw = np.genfromtxt(sample_data_name,unpack=True,skip_header=12,max_rows=1)
            # ---- create the open and sample data ----------------------------------------
            self.open_data   = pd.DataFrame(tt.optstat_group([ug_samp_data.ch*float(open_bw[0]),ug_open_data.counts],
                               FP,D,res_per_bin,E_cpts=E_cpts,verbose=verbosity,binedge=binedge).T,
                               columns=['tof','bin_width','counts','dcounts','tof_low'])
            self.sample_data = pd.DataFrame(tt.optstat_group([ug_samp_data.ch*float(samp_bw[0]),ug_samp_data.counts],
                               FP,D,res_per_bin,E_cpts=E_cpts,verbose=verbosity,binedge=binedge).T,
                               columns=['tof','bin_width','counts','dcounts','tof_low'])
        elif factor is not None and c_pts is not None and FP is not None:
            # -----------------------------------------------------------------
            # ---- Group with compression point grouping ----------------------
            # -----------------------------------------------------------------
            grouping_is_used = True
            # ---- reset the column names -------------------------------------------------
            self.datafile_names = ['ch','counts','dcounts']
            # ---- read in the fixed bin width data ---------------------------------------
            ug_open_data   = pd.read_csv(open_data_name,sep=r'\s+',names=self.datafile_names,
                                         skiprows=14)
            ug_samp_data = pd.read_csv(sample_data_name,sep=r'\s+',names=self.datafile_names,
                                         skiprows=14)
            # ---- read in the bin width --------------------------------------------------
            open_bw = np.genfromtxt(open_data_name  ,unpack=True,skip_header=12,max_rows=1)
            samp_bw = np.genfromtxt(sample_data_name,unpack=True,skip_header=12,max_rows=1)
            # ---- create the open and sample data ----------------------------------------
            self.open_data   = pd.DataFrame(tt.comp_group([ug_samp_data.ch*float(open_bw[0]),ug_open_data.counts],
                               factor,c_pts,FP,binedge=binedge).T,columns=['tof','bin_width','counts','dcounts'])
            self.sample_data = pd.DataFrame(tt.comp_group([ug_samp_data.ch*float(samp_bw[0]),ug_samp_data.counts],
                               factor,c_pts,FP,binedge=binedge).T,columns=['tof','bin_width','counts','dcounts'])
        else:
            err_str = "Data could not be read. If grouping is desired define factor, c_pts, and FP."
            print('If grouping, NORSUM file, constant bin width must be used.')
            raise ValueError(err_str)
        
        # read in the monitors and triggers
        self.open_mon   = pd.read_csv(open_mon_name  ,names=['mon'])
        self.sample_mon = pd.read_csv(sample_mon_name,names=['mon'])
        
        # ---------------------------------------------------------------------------------
        # ----- test data that was read just read -----------------------------------------
        # ---------------------------------------------------------------------------------
        # check if data is a 2-d array
        if (len(self.open_data.shape) > 2) or (len(self.sample_data.shape) > 2):
            raise ValueError("Data should be 2-d: columns and rows")
        # check if data has 3 columns (this is true for Ogrp data)
        if (self.open_data.shape[1] != 3) or (self.sample_data.shape[1] != 3):
            if factor == None and optstat_group == False:
                raise ValueError("Data should only have 3 columns")
        # check if user gave Ogrp files but also gave grouping input
        if grouping_is_used and 'Ogr' in open_data_name:
            raise ValueError('Name mismatch for grouping')
        # check if monitor files are 11 monitors
        # -- consider revising this requirement.
        if (len(self.open_mon) > 11) or (len(self.sample_mon) > 11):
            raise ValueError("There should be 11 monitors in the mon file")
        # check if there are extra columns in the monitor files
        if (self.open_mon.shape[1] > 1) or (self.sample_mon.shape[1] > 1):
            raise ValueError("The monitor file should be one column of 10 data")
            
        # if triggers is not specified by user, set it to first element of monitor array
        self.sample_trig = self.sample_mon['mon'][self.mon_index] 
        self.open_trig = self.open_mon['mon'][self.mon_index]
        
        # ---------------------------------------------------------------------------------
        # ----- convert the counts to count rate for each ---------------------------------
        # ---------------------------------------------------------------------------------
        self.calc_open_rate()
        self.calc_sample_rate()
        
        # ---------------------------------------------------------------------------------
        # ----- calculate bkg function along the tof grid ---------------------------------
        # ---------------------------------------------------------------------------------
        if bkgfunction == 'double_exp':
            self.open_data['bkg'] = self.double_exp(self.open_data.tof,self.o_bkgcoeff[0],
                                                    self.o_bkgcoeff[1],self.o_bkgcoeff[2],
                                                    self.o_bkgcoeff[3],self.o_bkgcoeff[4])
            self.sample_data['bkg'] = self.double_exp(self.sample_data.tof,self.s_bkgcoeff[0],
                                                        self.s_bkgcoeff[1],self.s_bkgcoeff[2],
                                                        self.s_bkgcoeff[3],self.s_bkgcoeff[4])
        elif bkgfunction == 'power_law':
            self.open_data['bkg'] = self.power_law(self.open_data.tof,self.o_bkgcoeff[0],
                                                   self.o_bkgcoeff[1],self.o_bkgcoeff[2])
            self.sample_data['bkg'] = self.power_law(self.sample_data.tof,self.s_bkgcoeff[0],
                                                       self.s_bkgcoeff[1],self.s_bkgcoeff[2])
        elif bkgfunction == 'double_power_law':
            self.open_data['bkg'] = self.double_power_law(self.open_data.tof,FP,
                                                    self.o_bkgcoeff[0],
                                                    self.o_bkgcoeff[1],self.o_bkgcoeff[2],
                                                    self.o_bkgcoeff[3],self.o_bkgcoeff[4])
            self.sample_data['bkg'] = self.double_power_law(self.sample_data.tof,FP,
                                                    self.s_bkgcoeff[0],
                                                    self.s_bkgcoeff[1],self.s_bkgcoeff[2],
                                                    self.s_bkgcoeff[3],self.s_bkgcoeff[4])
        elif bkgfunction == 'exp_power_mix':
            self.open_data['bkg'] = self.exp_power_mix(self.open_data.tof,self.o_bkgcoeff[0],
                                                    self.o_bkgcoeff[1],self.o_bkgcoeff[2],
                                                    self.o_bkgcoeff[3],self.o_bkgcoeff[4])
            self.sample_data['bkg'] = self.exp_power_mix(self.sample_data.tof,self.s_bkgcoeff[0],
                                                        self.s_bkgcoeff[1],self.s_bkgcoeff[2],
                                                        self.s_bkgcoeff[3],self.s_bkgcoeff[4])
        elif bkgfunction == 'exp_log':
            self.open_data['bkg'] = self.exp_log(self.open_data.tof,self.o_bkgcoeff[0],
                                                 self.o_bkgcoeff[1],self.o_bkgcoeff[2],
                                                 self.o_bkgcoeff[3])
            self.sample_data['bkg'] = self.exp_log(self.sample_data.tof,self.s_bkgcoeff[0],
                                                   self.s_bkgcoeff[1],self.s_bkgcoeff[2],
                                                   self.s_bkgcoeff[3])
        else:
            raise ValueError('Options for bkgfunction are: \'double_exp\',\'power_law\',',
                '\'double_power_law\',\'exp_power_mix\'')
        
        # ---------------------------------------------------------------------------------
        # -----calculate transmission -----------------------------------------------------
        # ---------------------------------------------------------------------------------
        self.calc_trans()

        # ---------------------------------------------------------------------------------
        # ----- finished with initialization ----------------------------------------------
        # ---------------------------------------------------------------------------------
        
    def calc_count_rate(self,COUNTS,BIN_WIDTH,TRIGGERS):
        """"
        Function to convert counts to count rate.
        
        Parameters
        ----------
        
        counts : array-like
            1-d vector of counts that were collected for each bin
        bin_width : array-like
            1-d vector of widths in micro-seconds for each bin associated
            with the counts in the counts vector
        triggers : float
        The number of times that the bin was "open" or collecting 
            counts data
        correct_deadtime : bool, optional
            If true, the counts data will be corrected for the dead time
            of the detection system.
        tau : float, optional
            Tau is the dead time of the data collection system in micro-
            seconds. If correct_deadtime is true, tau is required to correct
            for the dead time of the system
        
        Returns
        -------
        
        cps : array-like
            1-d vector of count rate in counts per second
        dcps : array-like
            1-d vector of the error in the count rate

        Notes
        -----
        .. math:: CPS = \\frac{C}{b_w\\cdot trig}
        
        """
        
        if self.correct_deadtime == True:
            if tau == None:
                raise ValueError("If dead time correction is wanted, specify tau.")
            # TODO : write the dead time correction
            else:
                raise ValueError("Doh! Dead time has not been implemented yet.")
            # correct counts for dead time
            
            # calculate count rate
        
            # calculate error on the count rate
            
        else:
            # calculate count rate
            cps = COUNTS/BIN_WIDTH/TRIGGERS/1e-06
            
            # calculate error on the count rate
            dcps = np.sqrt(COUNTS)/BIN_WIDTH/TRIGGERS/1e-06
        
            return cps,dcps
                
    def calc_open_rate(self):
        """
        Call calc_count_rate for the open_data
        """
        self.open_data['cps'],self.open_data['dcps'] = self.calc_count_rate(self.open_data.counts,
                                                self.open_data.bin_width,self.open_trig)
    def calc_sample_rate(self):
        """
        Call calc_count_rate for the sample_data
        """
        self.sample_data['cps'],self.sample_data['dcps'] = self.calc_count_rate(self.sample_data['counts'],
                                                self.sample_data['bin_width'],self.sample_trig)
        
        
        
    # ---------------------------------------------------------------------------------
    # ----- determine a bkg function for open and sample ------------------------------
    # ---------------------------------------------------------------------------------
    # --> for now I'll have a couple functions to adjust on the fly
        
    def double_exp(self,t,a,b,c,d,e):
        """
        .. math:: f(t) = (a)e^{-t(b)} + (c)e^{-t(d)} + e
        """
        return a*np.exp(-t*b) + c*np.exp(-t*d) + e
        
    def power_law(self,t,a,b,c):
        """
        .. math:: f(t) = (a)t^{-(b)} + c
        """
        return a*t**(-b) + c

    def double_power_law(self,t,FP,a1,b1,a2,b2,c):
        """
        .. math:: f(t) = (a1)t^{(b1)} + (a2)\\left(\\frac{72.3FP}{t}\\right)^{2(b2)} + c
        """
        return a1*t**b1 + a2*tt.tofe(t,FP)**b2 + c

    def exp_power_mix(self,t,a,b,c,d,e):
        """
        .. math:: f(t) = (a)e^{-t(b)} + (c)t^{-(d)} + e
        """
        return a*np.exp(-t*b) + c*t**(-d) + e

    def exp_log(self,t,a,b,c,d):
        """
        .. math:: f(t) = e^{ a + (b)t + c/ln(t) } + d
        """
        return np.exp(a + b*t + c/np.log(t)) + d
        
    # ---------------------------------------------------------------------------------
    # ----- calculate the transmission ------------------------------------------------
    # ---------------------------------------------------------------------------------
    def calc_trans(self):
        """
        .. math:: T = \\frac{\\dot C_s - \\dot B_s}{\\dot C_o - \\dot B_o}
        """
        self.trans = ((self.sample_data.cps - self.sample_data.bkg)/
                             (self.open_data.cps - self.open_data.bkg))

        self.dtrans = np.sqrt( 
            (self.sample_data.dcps/(self.open_data.cps-self.open_data.bkg))**2 +

            (self.samp_bkg_err*self.sample_data.bkg/(self.open_data.cps-self.open_data.bkg))**2 +
            
            ((self.sample_data.bkg-self.sample_data.cps)/(self.open_data.cps-self.open_data.bkg)**2 * 
                self.open_data.dcps)**2 +

            ((self.sample_data.cps-self.sample_data.bkg)/(self.open_data.cps-self.open_data.bkg)**2 * 
                self.open_data.bkg*self.open_bkg_err)**2
                              )

    def calc_sigma(self,N):
        """
        .. math:: \\sigma_t = \\frac{-1}{N}log(T) \\pm \\frac{\\Delta T}{NT}
        """
        # cross section for this samples transmission
        self.sig = -1/N*np.log(self.trans)
        # error on the cross section for the sample transmission
        self.dsig = self.dtrans/N/self.trans
        

def calc_cov(mon_list,dmon_list,room_bkgs,droom_bkgs,norm_list,dnorm_list,tof,sample_cps,
             sample_dcps,open_cps,open_dcps,bkg_function,bkg_pars,bkg_pars_cov,split=False,
             trans_idc_output=False):
    """
    Calculate covariance for transmission T = ( a1*Cs - a2*ks*B - B0s )/( a2*Co - a4*ko*B - B0o )

    To avoid code complexity, for now I'll force the user to specify everything up front instead
    of passing a Trans class. Available background functions: "power_law", "exp"

    Parameters
    ----------
    mon_list : list
        A list of the monitor values in order: [a1,a2,a3,a4]
    dmon_list : list
        A list of the uncertainty on monitor values in matching order
    room_bkgs : list
        A list of the room/ambient backgrounds a.k.a. B-zeros in order: [B0s,B0o]
    droom_bkgs : list
        A list of the uncertainty on room/ambient backgrounds in matching order
    norm_list : list
        A list of the normalizations for background in order: [ks,ko]
    dnorm_list : list
        A list of the uncertainty on normalizations for backgrounds in matching order
    tof : array-like
        A 1-dim array of time-of-flight for sample and open count rates. This must
        be "true" TOF, i.e. time-zero subtracted in order to get the background 
        correct.
    sample_cps : array-like
        A 1-dim array of sample count rate in order of time-of-flight (Cs)
    sample_dcps : array-like
        A 1-dim array of unc. on sample count rate in order of time-of-flight (dCs)
    open_cps : array-like
        A 1-dim array of open count rate in order of time-of-flight (Co)
    open_dcps : array-like
        A 1-dim array of unc. on open count rate in order of time-of-flight (Co)
    bkg_function : str
        The name of the background function used, matching the names in class Trans 
    bkg_pars : list
        The parameters needed for the background function (B-zero subtracted!)
    bkg_pars_cov : array-like
        A 2-dim array of the covariance for the fitted parameters.
    split : bool, optional
        Whether to return two parts of covariance instead of the sum, the statistical
        and the systematic
    trans_idc_output : bool, optional
        Extra output to be used with function `write_sammy_idc_file()`. In addition to
        default output `cov`, 3 extra variables are returned: `stat`, `syst_full`, and
        `sys_der_list`. 
    
    Returns
    -------
    cov : array-like
        A 2-dim numpy array with full covariance
    stat,syst : tuple
        A tuple of 2 2-dim arrays, the statistical and systematic covariance

    """
    available_bkg_funcs = np.array(['power_law','exp'])
    if not (available_bkg_funcs == bkg_function).any():
        raise ValueError("background function not available")
    # ----------------------------
    # Rename to simplify
    # ----------------------------
    ks,ko,dks,dko = norm_list[0],norm_list[1],dnorm_list[0],dnorm_list[1]
    a1,a2,a3,a4 = mon_list[0],mon_list[1],mon_list[2],mon_list[3]
    da1,da2,da3,da4 = dmon_list[0],dmon_list[1],dmon_list[2],dmon_list[3]
    b0s,b0o,db0s,db0o = room_bkgs[0],room_bkgs[1],droom_bkgs[0],droom_bkgs[1]
    cs,co,dcs,dco = sample_cps,open_cps,sample_dcps,open_dcps
    N,D,bfit = 1,1,0
    # ----------------------------
    # Bkg specific variables 
    # ----------------------------
    if bkg_function == 'power_law':
        A,B = bkg_pars[0],bkg_pars[1]
        bfit = A*tof**(-B)
        N = (a1*cs-a2*ks*bfit-b0s)
        D = (a3*co-a4*ko*bfit-b0o)
        dda = (-D*ks + N*ko)*tof**-B/D**2 
        ddb = (-D*ks + N*ko)*-1*bfit*np.log(tof)/D**2
        fit_ders = [dda,ddb]
    if bkg_function == 'exp':
        A,B = bkg_pars[0],bkg_pars[1]
        bfit = A*np.exp(-tof*B)
        N = (a1*cs-a2*ks*bfit-b0s)
        D = (a3*co-a4*ko*bfit-b0o)
        dda  = -1*(ks*D+ko*N)/(A*D**2)
        ddb  = (ks*D+ko*N)*bfit*tof/D**2
        fit_ders = [dda,ddb]
    # ----------------------------
    # Bkg agnostic variables
    # ----------------------------
    # statistical
    ddcs = a1/D
    ddco = -a3*N/D**2
    # systematic
    ddks = -a2*bfit/D
    ddko = a4*N*bfit/D**2
    ddb0s = -1/D
    ddb0o = 1*N/D**2
    dda1 = cs/D
    dda2 = -ks*bfit/D
    dda3 = -co*N/D**2
    dda4 = ko*bfit*N/D**2
    # statistical err and derivs *order matters!*
    stat_err_list = [dcs,dco]
    stat_der_list = [ddcs,ddco]
    # systematic err and derivs *order matters!*
    sys_err_list = [dks,dko,db0s,db0o,da1,da2,da3,da4] # doesn't inc. A,B they come with covar mat
    sys_der_list = np.concatenate([fit_ders,[ddks,ddko,ddb0s,ddb0o,dda1,dda2,dda3,dda4]])
    
    stat = mt.stat_cov(stat_err_list,stat_der_list)
    syst = mt.sys_cov(sys_err_list,sys_der_list,cov_mat=bkg_pars_cov) # Note bkg. cov. matrix here
    
    cov = stat+syst

    # ----------------------------
    # Return options
    # ----------------------------
    if trans_idc_output:
        # Form the full syst par cov 
        K = len(bkg_pars_cov)+len(sys_err_list)
        syst_full = np.zeros((K,K))
        syst_full[0:len(bkg_pars_cov),0:len(bkg_pars_cov)] = bkg_pars_cov
        for i in range(len(sys_err_list)):
            index = i+len(bkg_pars_cov)
            syst_full[index,index] = sys_err_list[i]**2

        return cov,stat,syst_full,sys_der_list
    if split:
        return stat,syst
    else:
        return cov

def write_covar_file(filename,energy,t,dt,sys_err_list,stat_err_list,sys_der_list,stat_der_list,
                 sys_err_str_list,ptwise_str_list,sys_cov=None,sys_cov_str_list=None,
                 high_precision=False):
    """
    Write a file with the information necessary to calculate the point by point covariance
    of a function.

    Parameters
    ----------
    filename : str
        The full file path to the file that is being written to
    energy : array
        The 1-d vector of energy. This should be the same length as all
        the other pointwise quantities given.
    t : array-like
        The observable, in this case transmission. Should be the same length as
        energy
    dt : array-like
        The uncertainty on the observable (for convenience). Should be the same
        length as energy, and should be the diagonal of the covariance matrix.
    sys_err_list : array-like
        The list of systematic errors for the function. These are constant 
        along the length of the function so are only single floats.
    stat_err_list : array-like
        A 2-d array where (len(stat_err_list) = no. of statistical variables),
        and (len(stat_err_list[0]) = len(energy)). This is a list of vectors 
        of the pointwise statistical error for each statistical variable.
    sys_der_list : array-like
        A list very similar to stat_err_list, but the number of vectors 
        corresponds to the number of systematic variables in the function.
        Each vector in this list is the derivative of the function with 
        respect to that variable. **The order of the derivative lists
        (systematic and statistical) should match the order of the error
        lists, in terms of which variables are listed**. 
    stat_der_list : array-like
        A list very similar to sys_der_list, but contains derivatives of 
        the function with respect to the statistical variables. It bears
        repeating that  **The order of the derivative lists
        (systematic and statistical) should match the order of the error
        lists, in terms of which variables are listed**.
    sys_err_str_list : list
        A list of strings in the order of the given systematic errors. This
        is used to label the information in the file.
    ptwise_str_list : list
        A list of strings in the order of [stat_err_list],[sys_der_list],
        [stat_der_list] that will be used to label the columns of data 
        listed in the file.
    sys_cov : array-like, optional
        If there are correlated systematic input parameters (e.g. parameters
        in the fitted function for background in transmission) the user can 
        provide the covariance matrix here. The covariance matrix is put at 
        the first rows and columns of the input systematic covariance matrix
        so **care must be taken to ensure the order of the derivatives matches
        the new order of errors given the inserted variables in sys_cov**. As
        an example if i have [dcs] as sys_err_list and [[da**2,dadb],[dbda,db**2]]
        as the sys_cov, my sys_der_list will look like [dda,ddb,ddcs].
    sys_cov_str_list : list, optional
        If the sys_cov is given, label the columns in the file with this list. 
        Optional unless the sys_cov variable is given.
    high_precision : bool, optional
        If the precision written to the file is too low (default 1e-6), write
        the floats with higher precision (1e-12).

    Returns
    -------
    none : ascii file
        Returns a somewhat annotated file with all the information necessary to
        reproduce the point-by-point covariance of the function. This output file
        is readable by read_covar_file().


    """
    def sixe(x):
        return '{:10.6e}'.format(x)
    if high_precision:
        def sixe(x):
            return '{:10.12e}'.format(x)

    with open(filename,'w+') as f:
        f.write('######\n# Systematic errors\n')
        f.write('# ')
        for i,string in enumerate(sys_err_str_list):
            f.write(string)
            if i != len(sys_err_str_list)-1:
                f.write(',')
        f.write('\n######\n')
        for i,err in enumerate(sys_err_list):
            f.write(sixe(err))
            if i != len(sys_err_list)-1:
                f.write(',')
        f.write('\n')
        if sys_cov is not None:
            f.write('######\n# Correlated syst. errors\n')
            f.write('# ')
            for i,string in enumerate(sys_cov_str_list):
                f.write(string)
                if i != len(sys_cov_str_list)-1:
                    f.write(',')
            f.write('\n######\n')
            for i,row in enumerate(sys_cov):
                for j,elem in enumerate(row):
                    f.write(sixe(elem))
                    if j != len(row)-1:
                        f.write(',')
                f.write('\n')
        f.write('######\n# E,t,dt,[stat. error],[sys derivatives],[stat derivatives]\n')
        f.write('# ')
        for i,string in enumerate(ptwise_str_list):
            f.write(string)
            if i != len(ptwise_str_list)-1:
                f.write(',')
        f.write('\n######\n')
        for i in range(len(sys_der_list[0])):
            f.write(sixe(energy[i])+','+sixe(t[i])+','+sixe(dt[i])+',')
            for j in range(len(stat_err_list)):
                f.write(sixe(stat_err_list[j][i])+',')
            for j in range(len(sys_der_list)):
                f.write(sixe(sys_der_list[j][i])+',')
            for j in range(len(stat_der_list)):
                f.write(sixe(stat_der_list[j][i]))
                if j != len(stat_der_list)-1:
                    f.write(',')
            f.write('\n')

def read_covar_file(filename,stat_err_cols,sys_der_cols,stat_der_cols,corr_sys_input=0):
    """
    Read a covariance file written by the write_covar_file() function in nuctools

    Parameters
    ----------
    filename : str
        The full file path to the file being read
    stat_err_cols : int
        The number of columns of statistical error in the file. Should be the same as 
        the number of statistical variables in the fucntion.
    sys_der_cols : int
        The number of columns of systematic derivatives in the file. Should be the same
        as the number of systematic variables in the fucntion.
    stat_der_cols : int
        The number of columns of statistical derivatives in the file. Should be the 
        same as the number of statistical variables in the fucntion.
    corr_sys_input : int, optional
        The nubmer of correlated variables given in a covariance matrix for systematic
        variables. Default is 0, meaning there is no covariance matrix in the file.

    Returns
    -------
    cy : numpy array
        The point-by-point covariance of the function

    """

    # read in the systematic error list
    sys_err_list = np.genfromtxt(filename,unpack=True,skip_header=4,delimiter=',',max_rows=1)
    # is there a matrix for correlated systematic error inputs?
    if corr_sys_input > 0:
        sys_cov = np.genfromtxt(filename,skip_header=9,delimiter=',',max_rows=corr_sys_input)
    # read in the point-wise data for statistical error and systematic/statistical derivatives
    ptwise = np.genfromtxt(filename,skip_header=9+corr_sys_input+4,delimiter=',',unpack=True)
    
    stat_err_list = ptwise[3:3+stat_err_cols,:]
    sys_der_list = ptwise[3+stat_err_cols:3+stat_err_cols+sys_der_cols,:]
    stat_der_list = ptwise[3+stat_err_cols+sys_der_cols:3+stat_err_cols+sys_der_cols+stat_der_cols,:]
    
    stat = mt.stat_cov(stat_err_list,stat_der_list)
    if corr_sys_input > 0:
        syst = mt.sys_cov(sys_err_list,sys_der_list,cov_mat=sys_cov)
    else:
        syst = mt.sys_cov(sys_err_list,sys_der_list)

    E = ptwise[0]
    observable = ptwise[1]
    unc_observable = ptwise[2]
    
    cy = stat+syst
    
    return E,observable,unc_observable,cy




