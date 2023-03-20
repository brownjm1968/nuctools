import numpy as np
import glob
import time
import warnings
import os



__all__ = ['describe_norsum_header','Rpy_xdr']

def describe_norsum_header(norsum_file,t0):
    x = 1


class Rpy_xdr:
    
    """
    Read all the data files in a specified folder with a specified extension 
    to sum counts over all cycles and run monitor statistics.

    Python program designed to mimic the GUI based analysis program at RPI called
    RPI XDR.
    
    TODO: consider removing read_files() from initialization and making methods
    to read and write files to avoid summing every time.

    TODO: write monitor file

    TODO: include dead time correction with dtc from sample_create

    Parameters
    ----------
    folder : str
        The folder where the raw data files live
    ext : str
        The file extension type (e.g. '.mp')
    mon_line_num : int
        The line number where the monitors are specified. *This is very 
        TOF clock dependent*, but for the files generated at the RPI 
        linac ('*.mp' and '*.889') monitors are listed 5 monitor counts
        per line for 2 lines. The mon_line_num parameter should be set 
        to that first line.
    num_samp : int
        How many samples are contained in the files in the specified 
        folder
    samp_line_num : int
        The line number on which the sample number is specified. This 
        program assumes that the sample number is the third element 
        of an array split by the space character (e.g. 'cmline3=Sample #  1')
    raw_bw : float
        The width of the counting bins in time [us]
    min_tof : float
        The minimum cutoff on the TOF spectrum [us] used calculate the total 
        detector counts in the monitor ratios. Used to cut out the gamma 
        flash


    Attributes
    ----------
    monitors : list
        A list of numpy arrays. Each element is a numpy array which corresponds
        to the 10 monitor counts for each cycle for a given sample.
    det_sum : list
        A list of numpy arrays. Each element is a sum (over cycles) of TOF 
        dependent detector counts for a given sample.
    mon_std : list
        A list of numpy arrays. Each element corresponds to different samples,
        each containing the fractional standard deviation of the monitor to 
        detector sum ratio for all 10 monitors.

    Methods
    -------
    read_files
        Read the files for the specified extension and folder. This function
        reads the files and stores the monitors and det_sum for each sample 
        and cycle.
    mon_frac_std
        Calculates the fractional standard deviation of the monitor to
        detector sum ratios for each sample and stores this information
        in mon_std
    write_tof_files
        Writes the contents of each element (and therefore each sample) of 
        det_sum array to a separate file along with the associated time
        of flight.

    Notes
    -----
    Fractional standard deviations are computed by:

    .. math:: \\frac{\\sigma\\left(\\frac{C_{Mon}}{C_{Det}}\\right)}{\\langle\\frac{C_{Mon}}{C_{Det}}\\rangle}

    where :math:`\\sigma()` is the standard deviation as computed by the numpy package.

    """
    
    def __init__(self,folder,ext='mp',mon_line_num=43,num_samp=None,samp_line_num=40,
                 raw_bw=None,head_len=None,min_tof=0.0,verbose=False,mon_header=None):
        
        self.folder = folder
        self.ext = ext
        self.mon_line_num = mon_line_num
        self.samp_line_num = samp_line_num
        self.head_len = head_len
        self.mon_header = mon_header

        if num_samp is None:
            self.num_samp = 4
            warnings.warn("Setting number of samples to 4")
        else:
            self.num_samp = num_samp
            if verbose:
                print("Number of samples: {}".format(self.num_samp))
        
        if raw_bw is None:
            self.raw_bw = 0.0064
            warnings.warn("Setting bin width to default: 0.0064 us")
        else:
            self.raw_bw = raw_bw

        if head_len is None:
            self.head_len = 59
            warnings.warn("Setting header length to 59")

        if mon_header is None:
            self.mon_header = "m1 m2 m3 m4 m5 m6 m7 m8 m9 m10 det"
            warnings.warn("Setting monitor header to default.")

        self.min_tof = min_tof

        # ---------------------------
        # Failures
        # ---------------------------
        if os.path.isdir(folder) is False:
            raise ValueError("Specified folder does is not a directory.")
        
        # ---------------------------
        # read the files with specified extension
        # ---------------------------
        self.monitors = [None]*self.num_samp  # <-- Each index could be diff. len()
        self.det_sum  = [None]*self.num_samp  # <-- Each index could be diff. len()
        for samp_num in range(self.num_samp):
            self.read_files(samp_num,verbose)
        
        # ---------------------------
        # calculate the percent std of monitor/det_counts ratio
        # ---------------------------
        self.mon_std = [None]*self.num_samp
        for samp_num in range(self.num_samp):
            self.mon_frac_std(samp_num)
    
    def read_files(self,samp_num,verbose):
        
        """
        Grab all the files with the specified extension and
        sample number. Read the monitor counts and detector
        TOF counts and store them into memory
        
        """
        
        first_file = True
        
        for i,file in enumerate(np.sort(glob.glob('{}*.{}'.format(self.folder,self.ext)))):

            if verbose: 
                print("File {}, {}".format(i,file))
            
            # ---------------------------
            # continue if wrong sample
            # ---------------------------
            sample = np.genfromtxt(file,max_rows=1,skip_header=self.samp_line_num-1,comments='*')[2]
            if sample != samp_num:
                continue
            
            # ---------------------------
            # grab the 10 monitors listed on the two rows, flatten to one list
            # ---------------------------
            temp = np.genfromtxt(file,max_rows=2,skip_header=self.mon_line_num-1)[:,2:].flatten()
            
            # ---------------------------
            # read the detector tof counts
            # ---------------------------
            det = np.genfromtxt(file,skip_header=self.head_len)
            
            # ---------------------------
            # gate the tof of total detector counts for monitor normalization
            # ---------------------------
            tof_array = np.arange(len(det))*self.raw_bw
            det_counts = np.sum(det[tof_array>=self.min_tof])
            temp = np.hstack((temp,det_counts))
            
            # ---------------------------
            # place into memory
            # ---------------------------
            if first_file:
                # set the monitors list
                self.monitors[samp_num] = temp
                first_file = False
                # set the detector sum
                self.det_sum[samp_num] = det
            else:
                self.monitors[samp_num] = np.vstack((self.monitors[samp_num],temp))
                self.det_sum[samp_num] += det
        
        # ---------------------------
        # transpose to make first index of that samples monitors the number of cycles
        # ---------------------------
        self.monitors[samp_num] = self.monitors[samp_num].T
        
    def mon_frac_std(self,samp_num):
        """
        Calculate the percent stand. deviation for all monitors of a 
        specific sample: samp_num
        
        """
        # ---------------------------
        # set length to len of monitor array
        # ---------------------------
        self.mon_std[samp_num] = np.empty(len(self.monitors[samp_num]))
        self.mon_std[samp_num][:] = np.nan
        
        # ---------------------------
        # standard deviation of mon/det ratio over mean(mon/det ratio)
        # ---------------------------
        for i,mon in enumerate(self.monitors[samp_num]):
            self.mon_std[samp_num][i] = (np.std(mon/self.monitors[samp_num][10]) / 
                               np.mean(mon/self.monitors[samp_num][10]))
            
    def write_tof_files(self,output_folder,name_string):
        
        """
        Takes the whatever is in the detector sum attribute: det_sum and
        writes it to file.
        
        Parameters
        ----------
        output_folder : str
            The folder to write the data files to
        name_string : str
            The string that should be placed at the beginning of each samples
            data file (experiment title, isotope of interest, etc....)
        
        """
        
        # ---------------------------
        # function overwrites files
        # ---------------------------
        warnings.warn("This function will overwrite old files.")
        
        for samp_num in range(self.num_samp):
            
            # ---------------------------
            # Add tof vector and transpose array to go in columns
            # ---------------------------
            tof_array = np.arange(len(self.det_sum[samp_num]))*self.raw_bw
            save_arr = np.vstack((tof_array,self.det_sum[samp_num])).T

            # ---------------------------
            # Get the monitors
            # ---------------------------
            monstring = ""
            for mon in self.monitors[samp_num].sum(axis=1):
                monstring += "{} ".format(mon)

            # ---------------------------
            # write the monitors to file
            # ---------------------------
            np.savetxt('data/{}_mon_{}.dat'.format(name_string,samp_num),
                       self.monitors[samp_num].T,fmt="%d",
                       header=self.mon_header,comments='')
            
            # ---------------------------
            # Save with 6 pt. precision
            # ---------------------------
            np.savetxt(output_folder+"{}_samp_{}.dat".format(name_string,samp_num),
                       save_arr,fmt='%.6f',header=monstring)
                
                