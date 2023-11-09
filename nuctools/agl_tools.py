import time
import json
import numpy as np
import pandas as pd

__all__ = ['read_and_add','sum_tof','read_rejected','read_ecals','read_gain_adj']

def read_and_add(filename,hist_df,ecal,gain,llduld,wfcoeff,numadc=4,
                 unweighted=False,verbose=False,mca_df=None,aglgroup=False):
    """
    TODO: incorporate gain adjustment
    
    Read the "list-mode" binary data files and add the events
    to the histograms in the Pandas DataFrame `hist_df`

    Parameters
    ----------
    filename : str
        The full file path to the binary file
    hist_df : DataFrame
        Pandas DataFrame object with attributes (or columns) of `tof`, `adc1`,
        `adc2`, `adc3`, `adc4`. This DataFrame is filled with histograms
        for each channel designated by adc numbers in this function. **This 
        histogram DataFrame must be set up properly** for this function to 
        work. That means the `hist_df.tof` vector must have appropriate bin 
        widths.
    ecal : array-like
        A 2-dimensional array of N rows and 3 columns. Each column is for 
        parameter `A0,A1,A2`. The number of rows is determined by the number of
        ADCs in the data file. Accessing `A1` for ADC1, e.g., would be: 
        `ecal[0][1]`
    gain : array-like
        A 1-d array of gain values, one for every ADC. These are calibrated to match 
        value for beginning run set to source calibration.
    llduld : array-like
        A 2-dimensional array of N rows and 2 columns. Column 0 represents the 
        lower-level discriminator (LLD) and column 1 the upper-level discriminator
        (ULD). The number of rows is determined by the number of ADCs in the file
    wfcoeff : array-like
        A 1-d array for weighting function coefficients defined by the GELINA 
        weighting function. There are 8 parameters, :math:`a_i` from :math:`i=-3,4`
    numadc : int, optional
        The number of ADCs listed **in the binary file. This will change the way**
        **the binary file is read.** Often there will be either 2 or 4 ADCs, 4 is
        the default.
    unweighted : bool, optional
        Whether to add unweighted histograms to the DataFrame 
    verbose : bool, optional
        Whether to increase print output
    mca_df : DataFrame
        Pandas DataFrame object with attributes (or columns) of `ed`, `adc1`,
        `adc2`, `adc3`, `adc4`. This DataFrame is filled with histograms
        for each channel designated by adc numbers in this function. **This 
        histogram DataFrame must be set up properly** for this function to 
        work. `ed` is for "energy deposition"
    agl_group : bool
        Whether to group using bin-edges defined AGL-style (`True`), as opposed to 
        setting TOF compression points for bin edge definition RPI-style (`False`)

    Examples
    --------
    >>> import pandas as pd
    >>> import nuctools as nuc
    >>> user_max_tof = 2.5e3 # [us]
    >>> bin_width = 0.001    # [us]
    >>> numbins = int(user_max_tof/bin_width)+1
    >>> bin_width_arr = np.repeat(bin_width,numbins)
    >>> ecal = [[ 144.85,6.4355,1],
                [ 145.66,6.6452,1],
                [ 179.75,6.3147,1],
                [ 56.069,6.3898,1]] # [keV]
    >>> gain = [1,1,1,1]
    >>> llduld = [[16,1100],
                  [13,1160],
                  [12,1150],
                  [30,1150]] # [keV]
    >>> wfcoeff = [68.7,-513.4,1367.7,-1614.6,993.7,-159.5,16.5,-0.568] # [1/MeV]
    >>> filename = "example.lst"
    >>> data = pd.DataFrame({
            'tof'  : np.cumsum(bin_width_arr),
            'adc1' : np.zeros(numbins),
            'adc2' : np.zeros(numbins),
            'adc3' : np.zeros(numbins),
            'adc4' : np.zeros(numbins)
        })
    >>> nuc.read_and_add(filename,data,ecal,gain,llduld,wfcoeff)

    Notes
    -----

    Formula for the energy calibration is given by:

    .. math:: E_d = A_0 + A_1 \\cdot PA^{A_2}

    where :math:`PA` is pulse area and :math:`E_d` is in units of keV. The 
    weighting function is given by [1]:

    .. math:: W(E_d) = \\sum_{i=-3}^{4} a_i \\cdot E_d^i.

    The :math:`PA` is multiplied by the gain, as calibrated to a calibration source, 
    and logged during runs.

    [1] Borella et. al, The use of C6D6 detectors for neutron induced capture 
    cross-section measurements in the resonance region, Nuclear Instruments and
    Methods in Physics Research A 577 (2007) 626–640
    """

    #------------------------
    # Define binary format
    #------------------------
    # dt = [('tcmsb',np.dtype(">u2")),('tclsb',np.dtype(">u2")),('adc1',np.dtype(">u2")),
    #       ('adc2',np.dtype(">u2")),('adc3',np.dtype(">u2")),('adc4',np.dtype(">u2"))]
    dt = [('tcmsb',np.dtype(">u2")),('tclsb',np.dtype(">u2"))]
    for i in range(1,numadc+1):
        dt.append(('adc{}'.format(i),np.dtype(">u2")))
    #------------------------
    # Read single file
    #------------------------
    d = np.fromfile(filename,dtype=dt)
    
    dat = pd.DataFrame({
        'tc': np.left_shift((d['tcmsb'] & 0x01FF).astype(np.uint32),16) + d['tclsb'],
    })
    for i in range(1,numadc+1):
        dat['adc{}'.format(i)] = d['adc{}'.format(i)].astype(float)
    
    #------------------------
    # Bad TOFs
    #------------------------
    initial_num_events = dat.shape[0]
    dat = dat[dat.tc!=0]
    if(verbose):
        print("Rejected for TOF=0:       {:.3f}%; Number counts lost: {}".format(
            (initial_num_events-len(dat))/initial_num_events*100,initial_num_events-len(dat)))

    datdict = {}
    #------------------------
    # for every adc run these functions
    #------------------------
    for i in range(1,numadc+1):
        if( verbose ):
            if( (dat['adc{}'.format(i)]!=0).any() ):
                print("Found events for ADC{}".format(i))
        #------------------------
        # reject coincidence
        #------------------------
        temp = dat[(dat.loc[:,~dat.columns.isin(['tc'])]>0).sum(axis=1)==1]
        #------------------------
        # split into adc's
        #------------------------
        temp = temp[['tc','adc{}'.format(i)]][(temp['adc{}'.format(i)]!=0)].reset_index(drop=True).copy()
        #------------------------
        # save mca before lld and uld cuts (but bin-structure can still cut!) and apply gain
        #------------------------
        if gain[i-1] == 1.0:
            datdict["mca{}".format(i)] = temp['adc{}'.format(i)]
        else:
            # add uniform rand [0,1) if gain != 1 (prevent non-physical zeros)
            rnd = np.random.rand(len(temp['adc{}'.format(i)]))
            temp['adc{}'.format(i)] = gain[i-1]*(temp['adc{}'.format(i)] + rnd)-1
            datdict["mca{}".format(i)] = temp['adc{}'.format(i)]
        #------------------------
        # apply lower and upper level discriminators
        #------------------------
        temp = temp[(temp['adc{}'.format(i)] > llduld[i-1][0]) & (temp['adc{}'.format(i)] < llduld[i-1][1])]
        #------------------------
        # Calibrate pulse area (PA) to energy deposited [keV]
        #------------------------
        temp['Ed'] = (ecal[i-1][0] + ecal[i-1][1]*(temp['adc{}'.format(i)])**ecal[i-1][2])
        #------------------------
        # apply weighting function to get the weighted counts (coeff are [1/MeV])
        #------------------------
        temp['Cw'] = np.zeros(len(temp['Ed']))
        for j in range(-3,4+1):
            temp.Cw += wfcoeff[j+3]*(temp['Ed']/1000)**(j)
        datdict["d{}".format(i)] = temp
        
    #------------------------
    # add to histograms
    #------------------------
    max_tof_ns = hist_df.tof.max() # should be in [ns]
    numbins = len(hist_df)
    bin_def = numbins
    binrange = (0,max_tof_ns)
    if( aglgroup ): 
        last_tof = 2*hist_df.tof[numbins-1]-hist_df.tof[numbins-2]
        bin_def = np.insert(np.array(hist_df.tof),numbins,last_tof)
        bin_range = (0,last_tof)

    #------------------------
    # add weighted
    #------------------------
    for i in range(1,numadc+1):
        hist_df['adc{}'.format(i)] += np.histogram(datdict["d{}".format(i)].tc,weights=datdict["d{}".format(i)].Cw,bins=bin_def)[0]

    #------------------------
    # add unweighted if requested
    #------------------------
    if( unweighted ):
        for i in range(1,numadc+1):
            hist_df['uwadc{}'.format(i)] += np.histogram(datdict["d{}".format(i)].tc,bins=bin_def)[0]

    #------------------------
    # add mca if requested
    #------------------------
    if( mca_df is not None ):
        numbins = len(mca_df)
        for i in range(1,numadc+1):
            mca_df['adc{}'.format(i)] += np.histogram(datdict["mca{}".format(i)],bins=numbins,range=(0,numbins))[0]

    
def sum_tof(agl_inp_file,file_list):
    """
    Sum the TOF histograms from all events in the file list into a Pandas
    DataFrame and return the DataFrame

    Parameters
    ----------
    agl_inp_file: str
        The file name that describes necessary input for reading list-mode 
        files like AGL
    file_list : list
        A list of strings with full filepaths. Suggest to glob with Python: 
        `glob.glob("my/directory/*.lst")`
    
    Returns
    -------
    data : DataFrame
        A Pandas DataFrame with column names: "tof", "adc1", "adc2",
        "adc3", "adc4"
    data,mca : tuple
        A tuple of Pandas DataFrames with columns as above for `data` and columns
        of "ed", "adc1", "adc2", "adc3", "adc4" for `mca`

    Examples
    --------
    >>> import glob
    >>> import nuctools as nuc
    >>> filelist   = glob.glob('La/la_fp4a_t03_r03_*.lst')
    >>> 
    >>> data = nuc.sum_tof("test-files/test-agl-inp.json",filelist)

    """

    # -----------------------------------------------------
    # Read dictionary into memory and move to variables
    # -----------------------------------------------------
    with open(agl_inp_file,'r') as f:
        agldict = json.load(f)
    numadc = agldict['numadc']
    llduld = agldict['llduld']
    mcabins = agldict['mcabins']
    unweighted = agldict['unweighted']
    verbose = agldict['verbose']
    aglgroup = agldict['useAGLgrouping']
    ecal = agldict['ecal']
    gain_adj = agldict['gain_adj']
    wfcoeff = agldict['wfcoeff']
    bin_width = agldict['bin_width']
    user_max_tof = agldict['max_tof']
    cfct = agldict['cfct']
    zones = agldict['zones']
    rundict = agldict['rundict']


    print("Total number of files: ",len(file_list))
    if(verbose):
        print("Number of ADCs: ",numadc)
        print("Number bins for MCA: ",mcabins)

    # ------------------------------
    # create the TOF histo
    # - internally we use only [ns] to allow integer math
    # - TODO: move all input to YAML file
    # ------------------------------
    bin_width *= 1000 # go to [ns]
    if( aglgroup ):
        numbins = np.sum(zones*1024)
        toflist = np.array([0])
        for i,mult in enumerate(cfct):
            zone_array = np.cumsum(np.repeat(2**mult*bin_width,zones[i+1]*1024))
            if i==0:
                toflist = np.insert(zone_array,0,0.0)
            else:
                toflist = np.concatenate([toflist,zone_array+toflist[len(toflist)-1]])
        data = pd.DataFrame({
            'tof'  : toflist,
            'adc1' : np.zeros(len(toflist)),
            'adc2' : np.zeros(len(toflist)),
        })
        data.drop(data.tail(1).index,inplace=True)
    else:
        numbins = int(user_max_tof/bin_width)+1
        bin_width_arr = np.repeat(bin_width,numbins)
        bin_width_arr[0] = 0 # begin first tof bin at zero
        data = pd.DataFrame({
            'tof'  : np.cumsum(bin_width_arr),
        })
    for i in range(1,numadc+1):
        data['adc{}'.format(i)] = np.zeros(numbins)
        if( unweighted ):
            data['uwadc{}'.format(i)] = np.zeros(numbins)
    if( verbose ):
        print("Number bins for TOF: ",numbins)
        print(data.head())
    # ------------------------------
    # create the MCA histo
    # ------------------------------
    mca = pd.DataFrame({
        'ed' : np.linspace(0,np.max(llduld),mcabins)
        })
    for i in range(1,numadc+1):
        mca['adc{}'.format(i)] = np.zeros(mcabins)
    if( mcabins==0 ):
        mca = None
    # ------------------------------
    # Sum over the files in a folder
    # ------------------------------
    rundict_keys = np.sort(list(rundict.keys()))
    start = time.time()
    i,j,k=0,0,0
    for filename in np.sort(file_list):
        bad_run_found = False
        for key in rundict_keys:
            if key in filename:
                for bad_run in rundict[key]:
                    if( "{:0>4d}".format(bad_run) in filename and key in filename ):
                        j+=1
                        bad_run_found = True
        if(bad_run_found):
            continue
        # set gain values to proper run
        for k,key in enumerate(rundict_keys):
            if key in filename:
                gain = gain_adj[k]
        if(verbose):
            print("-------------------------\nFile: {}\n-------------------------".format(filename))
            print("Gain: ",gain)
        if(i==10):
            time_per_file = (time.time()-start)/10
            print("Time/file at 10 files: {:.2f} secs".format(time_per_file))
        if(i==len(file_list)//2):
            print("50% done.")

        read_and_add(filename,data,ecal,gain,llduld,wfcoeff,numadc,unweighted,verbose,mca,aglgroup)
        i+=1
    end = time.time()
    print("    Cycles thrown out: {}/{}".format(j,len(file_list)))
    print("  Total time required: {:.2f} secs".format(end-start))

    if( mcabins>0 ):
        return data,mca
    else:
        return data

def read_rejected(agl_outfile,expid):
    """
    Read existing AGL output file for rejected run IDs 
    and the rejected run, place into string for pyagl

    Parameters
    ----------
    agl_outfile : str
        Full path to AGL report file, typically ending with "*rep.txt"
    expid : str
        unique string for run ids, for example "zr91_fp14a" when the 
        run id strings are "zr91_fp14a_c01_r01", etc.

    Returns
    -------
    json_string : str
        A string in json format that can be copy and pasted into the
        pyagl input file to match the same behavior as AGL

    Examples
    --------
    >>> import nuctools
    >>> agl_file = "ornl_fp14_WF_Fe50mm_2022_BCoNaS_an_rep.txt"
    >>> agl_file_path = "/some/dir/"
    >>> expid = "fe50_fp14a"
    >>> rejected_string = nuc.read_rejected(agl_file_path+agl_file,expid)

    """
    with open(agl_outfile,'r') as f:
        lines = f.readlines()
    rejected_dict = {}
    for line in lines:
        if "Runtime" in line:
            break
        if expid in line:
            runid = line.split()[0]
            rejected_dict[runid] = []
        if "rejected" in line:
            reject_list = line.split()
            reject_list = np.array(reject_list[1:len(reject_list)]).astype(int).tolist()
            rejected_dict[runid] = reject_list
    return json.dumps(rejected_dict,sort_keys=True)

def read_ecals(agl_outfile,numadc,numruns):
    """
    Read the energy calibrations for every adc and every run

    Parameters
    ----------
    agl_outfile : str
        Full path to AGL report file, typically ending with "*rep.txt"
    numadc : int
        The number of ADCs that were used in the experiment
    numruns : int
        The number of runs that have energy calibrations listed. Keep 
        in mind AGL produces 1 extra listed set of parameters 
    """
    with open(agl_outfile,'r') as f:
        lines = f.readlines()
    adcsfound = 0
    runsfound = 0
    numparam = 3
    ecal = np.zeros((numruns,numadc,numparam)) # no error catch for def too small
    # load the ecals for every adc for every run, skip last listed (repeat)
    for line in lines:
        # break if extra run found
        if "1EN_CAL_A0" in line:
            runsfound+=1
            if runsfound > numruns:
                break
        if "EN_CAL_A0" in line:
            adc = int(line[0])
            if adc>adcsfound:
                adcsfound = adc
            a0 = float(line.split()[1])
            ecal[runsfound-1,adc-1,0] = a0
        if "EN_CAL_A1" in line:
            adc = int(line[0])
            if adc>adcsfound:
                adcsfound = adc
            a1 = float(line.split()[1])
            ecal[runsfound-1,adc-1,1] = a1
        if "EN_CAL_A2" in line:
            adc = int(line[0])
            if adc>adcsfound:
                adcsfound = adc
            a2 = float(line.split()[1])
            ecal[runsfound-1,adc-1,2] = a2
    # test if ADCs is unexpected
    if adcsfound != numadc:
        print("ADCs found   = ",adcsfound)
        print("ADCs defined = ",numadc)
        raise ValueError("Number of ADCs defined is not equal to ADCs found.")
    print(ecal)

def read_gain_adj(agl_outfile,numadc,numruns):
    """
    Read the pulse area (MCA) gain adjustments for every adc and every run

    Parameters
    ----------
    agl_outfile : str
        Full path to AGL report file, typically ending with "*rep.txt"
    numadc : int
        The number of ADCs that were used in the experiment
    numruns : int
        The number of runs that have energy calibrations listed. Keep 
        in mind AGL produces 1 extra listed set of parameters 
    """
    with open(agl_outfile,'r') as f:
        lines = f.readlines()
    gain_adj = np.zeros((numruns,numadc))
    runsfound = 0
    for line in lines:
        if "1AMP_GAIN_ADJUSTMENT" in line:
            runsfound += 1
        if runsfound > numruns:
            break
        if "AMP_GAIN_ADJUSTMENT" in line:
            adc = int(line[0])-1
            gain = line.split()[1]
            gain_adj[runsfound-1,adc] = gain
    np.set_printoptions(precision=5)
    print(gain_adj)

