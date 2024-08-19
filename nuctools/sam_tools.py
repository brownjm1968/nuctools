import pandas as pd
import numpy as np
import h5py as h5
from . import tof_tools as tt
from . import math_tools as mt
from . import funky
import time

__all__ = ['read_h5rpcm','read_h5xscm','read_pds','write_sammy_idc_file',
           'print_chann','create_sam_inp','fmt_twenty','fmt_ten','fmt_par',
           'switch_par_flags','describe_norsum_header','LPTtoCOV','mini_cov']


def read_h5rpcm(filename):
    """
    New SAMMY covariance file is in HDF5 in the form of an upper 
    triangular matrix. Return covariance, parameters, etc.

    Parameters
    ----------
    filename : str
        The path and filename for the SAMMY covariance file in HDF5 format

    Returns
    -------
    ucov,upar,covind,ispup : tuple
        A tuple of numpy arrrays with the reduced parameters covariance matrix, 
        reduced parameters, indices for parameters into the covariance matrix, 
        and a list of 0,1 booleans on whether the parameter was PUP'd
    """
    with h5.File(filename,"r+") as f:
        flatcov = np.array(f['covar'])
        upar = np.array(f['param'])
        covind = np.array(f['covind'])
        try: 
            ispup = np.array(f['ispup'])
        except:
            ispup = None

    lp = len(upar)
    ucov = np.zeros((lp,lp))
    k = 0
    for i in range(lp):
        for j in range(i,lp):
            ucov[i,j] = flatcov[k]
            ucov[j,i] = flatcov[k]
            k+=1

    if ~np.all(ucov == ucov.T):
        raise ValueError("Covariance matrix is not symmetric!")

    return ucov,upar,covind,ispup


def read_h5xscm(filename):
    """
    New SAMMY XS covariance file is in HDF5 in the form of a lower
    triangular matrix (different than RPCM above!). Return cross section
    covariance

    Parameters
    ----------
    filename : str
        The name of the full file path to the H5 covariance file
    Returns
    -------
    xscov : array-like
        The 2-d symmetric array of covariance on the theoretical observable
        in SAMMY.
    """

    with h5.File(filename,"r+") as f:
        try: 
            flatxscov = np.array(f['xs_cov'])
        except:
            flatxscov = None
    
    lxc = len(flatxscov)
    lp = int(np.floor(np.sqrt(1+8*lxc)-1)/2)
    xscov = np.zeros((lp,lp))
    k = 0
    for i in range(lp):
        for j in range(0,i+1):
            xscov[i,j] = flatxscov[k]
            xscov[j,i] = flatxscov[k]
            k+=1

    if ~np.all(xscov == xscov.T):
        raise ValueError("Covariance matrix is not symmetric!")

    return xscov

def read_pds(pds_file):
    """
    Read the PDS partial derivative file and load it into a pandas dataframe. 

    Parameters
    ----------
    pds_file : str
        The path to the PDS file

    Returns
    -------
    df : DataFrame
        A pandas data frame, with columns: 
        "expdata","expunc","theory","dtdp1","dtdp2","dtdp3","dtdp{}",...
        all the way to the number of parameters used to create the 
        cross section
    """
    with open(pds_file,"r") as f:
        lines = f.readlines()
    evenlinesplit = []
    deriv_table = []
    pandas_cols = ["expdata","expunc","theory","dtdp1","dtdp2","dtdp3"]
    for i,line in enumerate(lines):
        if i==0:
            # first line
            numpars = int(line)
            print("Number of pars = ",numpars)
            for parnumber in range(4,numpars+1):
                pandas_cols.append("dtdp{}".format(parnumber))
            print(pandas_cols)
            continue
        if i==1:
            # second line
            reduced_pars = np.array(line.split()).astype(float)
            print("Reduced pars = ",reduced_pars)
            continue
        if (i-2)%2 == 0:
            # odd lines (as shown in SAMMY manual, but even when 0-indexed)
            evenlinesplit = line.split()
        if (i-2)%2 == 1:
            # even lines (as shown in SAMMY manual, but odd when 0-indexed)
            oddlinesplit = line.split()
            evenlinesplit.extend(oddlinesplit)
            newrow = np.array(evenlinesplit).astype(float)
            deriv_table.append(newrow)
    df = pd.DataFrame(deriv_table,columns=pandas_cols)
    return df

def write_sammy_idc_file(filename_abbrev,energy,obs,obserr,stat,syst,sys_der_list,high_precision=False):
    """
    Write a file with the information requested by the SAMMY code for calculating covariance.
    This function will also write the "data file" needed by SAMMY

    Parameters
    ----------
    filename_abbrev : str
        The full file path to the file that is being written to, WITHOUT the extension
        (since two files are written with default extensions)
    energy : array-like
        The 1-d vector of energy. This should be the same length as all
        the other pointwise quantities given.
    obs : array-like
        The observable, should be the same length as energy
    obserr : array-like
        The uncertainty on the observable (for convenience). Should be the same
        length as energy, and should be the diagonal of the total covariance matrix.
    stat : array-like
        The statistical covariance matrix, can be 2-d or 1-d as only the diagonal is
        non-zero.
    syst : array-like
        The user can provide the systematic covariance matrix here. It should be 
        noted that **care must be taken to ensure the order of the covariance matches
        the order of derivatives given in sys_der_list**. 
    sys_der_list : array-like
        Each vector in this list is the derivative of the function with 
        respect to that variable. **This is the Jacobian, and the order
        of the derivatives must be the same as the order given in the 
        covariance matrix syst**. 
    high_precision : bool, optional
        If the precision written to the file is too low (default 1e-6), write
        the floats with higher precision (1e-12).

    Returns
    -------
    none : ascii file
        Returns a somewhat annotated file with all the information necessary to
        reproduce the point-by-point covariance of the function, in the format
        SAMMY expects.


    """
    def sixe(x):
        return '{:10.6e}'.format(x)
    if high_precision:
        def sixe(x):
            return '{:10.12e}'.format(x)

    # use numpy 
    energy = np.array(energy)
    stat = np.array(stat)
    syst = np.array(syst)
    sys_der_list = np.array(sys_der_list)
    
    # --- Data file ---
    fmt_twenty([energy,obs,obserr],filename_abbrev+'.twenty')

    # --- IDC file ---
    if len(np.shape(stat))>1:
        stat_err = np.sqrt(np.diag(stat))
    else:
        stat_err = np.sqrt(stat)
    with open(filename_abbrev+'.idc','w+') as f:
        f.write('Number of data-reduction parameters = '+str(len(sys_der_list))+'\n')
        f.write('\n')
        f.write('Free-format partial derivatives\n')
        for i in range(len(sys_der_list[0])):
            f.write(sixe(energy[i])+' '+sixe(stat_err[i])+' ')
            for j in range(len(sys_der_list)):
                f.write(sixe(sys_der_list[j][i]))
                if j != len(sys_der_list)-1:
                    f.write(' ')
            f.write('\n')
        f.write('\n')
        f.write('uncertainties on data-reduction parameters\n')
        for i,err in enumerate(np.sqrt(np.diag(syst))):
            f.write(sixe(err))
            if i != len(syst)-1:
                f.write(' ')
        f.write('\n')
        f.write('\n')
        
        off_diagonal = False
        if ( (np.eye(len(syst))*syst - syst).any() > 0):
            off_diagonal = True
            syst_corr = mt.cov_to_corr(syst)
        if off_diagonal:
            f.write('Correlation for data-reduction parameters\n')
            for i,row in enumerate(syst_corr):
                if i==0:
                    continue
                for j,elem in enumerate(row):
                    if j==i:
                        break
                    f.write(sixe(elem))
                    if j != len(row)-1:
                        f.write(' ')
                f.write('\n')
        f.write('\n')

def print_chann(cpts,factor,base_width,t0,FP,uncertainty="0.8"):
    """
    Print the CHANN input for SAMMY resolution function

    Print the "crunch" boundaries in energy [eV], the 
    grouping factor of base bin width, the bin width [ns],
    and the error for for the bin width.

    Parameters
    ----------
    cpts : list
        A list of "crunch" boundaries (tof locations [us]) used
        to group time-of-flight data.
    factor : list
        The factors to group by on each side of the crunch boundaries.
        The len(factor) should be 1 greater than len(cpts).
    base_width : float
        The base width [ns] of the time bins that have been grouped.
    t0 : float
        The time at which neutrons are emitted from a target [us].
    FP : float
        The length of the path neutrons will be traveling, i.e. the
        'flight path' [m]
    uncertainty : str
        The uncertainty on the bin width [ns]

    """

    if (len(factor) < len(cpts)) or (len(factor) == len(cpts)):
        raise ValueError("You should have more grouping factors than cpts")
    if base_width < 0.1:
        print("Warning: Your base width is < 0.1 ns.",
              "Units must be NANOSECONDS.")
    print("--- Paste into SAMMY res. function ---")
    for i,tof in enumerate(cpts[::-1]):
        print("CHANN 0 ",funky.two(tt.tofe(tof-t0,FP)),
              factor[::-1][i]*base_width,uncertainty)
    print("--------------------------------------")
    return 0

def create_sam_inp(inp_list,filename):
    """
    Create a file to automate the interactive input in SAMMY
    which specifies the various input files and energy regions.

    Parameters
    ----------
    inp_list : list
        A list of strings. Each string should be the full line 
        given to SAMMY as input. Some examples of strings are:
        "inp/samndf.inp" 
        "par/samndf.par"
        "data/trans_ta1.dat 200. 300."
    filename : str
        The name of the file that the strings will be written to.
        Each element of inp_list is written to the file and a 
        newline is added after each string.

    Returns
    -------
    no return : ascii file
    """

    with open(filename,'w+') as f:
        for line in inp_list:
            f.write(line+'\n')

def fmt_twenty(data,filename):
    """
    Function designed to take a data array that contains
    transmission or yield data and create a file that is 
    compatible with the SAMMY code TWENTY format

    Parameters
    ----------
    data : numpy array
        Should be in the format of data[0] = energy 1-d
        array, data[1] = trans/yield 1-d array, 
        data[2] = error on t/y 1-d array

    Returns
    -------
    no return : ascii file
    """
    a = np.array(data).T
    np.savetxt(filename,a,fmt="%20.10f",delimiter="") 

def fmt_ten(data,filename,tot=True):
    """
    Function designed to take a data array that contains
    total or capture cross section data and create a file 
    that is compatible with the FITACS data format.

    Parameters
    ----------
    data : numpy array
        Should be in the format of data[0] = energy 1-d
        array, data[1] = trans/yield 1-d array, 
        data[2] = error on t/y 1-d array
    filename : str
        The full file path or relative
    tot : bool
        Whether it is total cross section to be written, 
        if not, capture cross section will be written

    Returns
    -------
    no return : ascii file
    """
    a = np.array(data).T
    np.savetxt(filename,a,fmt="%10.3e",delimiter="")

    with open(filename,'r+') as f:
        file = f.read()
        f.seek(0,0)
        if tot:
            f.write('TOTAL CROSS SECTION\nABSOLUTE UNCERTAINTY\n'+file)
        else:
            f.write('CAPTURE CROSS SECTION\nABSOLUTE UNCERTAINTY\n'+file)

def fmt_par(data,filename):
    """
    Format parameter values into PAR format readable by SAMMY
    
    Parameters
    ----------
    data : array like
        An array with columns of data: energy, \Gamma_\gamma, 
        \Gamma_n, \Gamma_{x1}, \Gamma_{x2}, vary energy, vary
        \Gamma_\gamma, vary \Gamma_n, vary \Gamma_{x1}, vary
        \Gamma_{x2}, spin group. There should be 11 fields
    filename : str
        The file name to write the formatted data into.
        
    Returns
    -------
    no return : ascii file
    """
    
    a = np.array(data).T
    if len(a[0]) < 11:
        ValueError('Requires 11 fields for each resonance energy.')
    with open(filename,'w+') as f:
        for row in a:
            f.write('{:0<11.5f} {:0<10.5f} {:0<10.5f} '.format(row[0],row[1],row[2]))
            f.write('{:0<10.5f} {:0<10.5f} {:1d} '.format(row[3],row[4],int(row[5])))
            f.write('{:1d} {:1d} {:1d} '.format(int(row[6]),int(row[7]),int(row[8])))
            f.write('{:1d} {:1d}\n'.format(int(row[9]),int(row[10])))



def switch_par_flags(par_file,fname_mod,flag_id,flag_val,erange=None):
    """
    Opens par_file and changes the specified flag_id to the
    flag_val. This allows the user to change 0's to 1's or 
    vice versa to turn on/off fitting for each resonce width.

    Parameters
    ----------
    par_file : str
        The address of a PAR file that is in the format typical for
        the code SAMMY.
    fname_mod : str
        String to be added to the end of the par_file name
    flag_id : str
        Either a single str or list of str's that match ve = vary
        energy, vgg = vary gamma gamma, vc1 = vary channel 1, 
        vc2 = vary channel 2, vc3 = vary channel 3.
    flag_val : int
        Either a 0,1,3 for: don't fit, fit, and PUP respectively.
    erange : 1-d array
        An array of length = 2, where erange[0] = lowest energy
        for flags to be changed, and erange[1] = highest energy
        for flags to be changed.

    Returns
    -------
    no return : ascii file
        Creates a new PAR file with only the specified flags changed.
    """
    ve,vgg,vc1,vc2,vc3 = False,False,False,False,False
    if 've' in flag_id:
        ve = True
    if 'vgg' in flag_id:
        vgg = True
    if 'vc1' in flag_id:
        vc1 = True
    if 'vc2' in flag_id:
        vc2 = True
    if 'vc3' in flag_id:
        vc3 = True
    with open(par_file,"r+") as f:
        initial_par = f.readlines()
    with open(par_file+fname_mod,"w+") as f:
        for line in initial_par:
            line_long_enough = False
            inside_erange = False
            if len(line) > 65:
                line_long_enough = True
            if line_long_enough:
                energy = float(''.join(line[0:10]))
                if erange == None:
                    inside_erange = True
                elif energy > erange[0] and energy < erange[1]:
                    inside_erange = True
                if erange == None or inside_erange:
                    line_list = list(line)
                    if ve == True:
                        line_list[56] = str(flag_val)
                    if vgg == True:
                        line_list[58] = str(flag_val)
                    if vc1 == True:
                        line_list[60] = str(flag_val)
                    if vc2 == True:
                        line_list[62] = str(flag_val)
                    if vc3 == True:
                        line_list[64] = str(flag_val)
                    line_join = ''.join(line_list)
                    f.write(line_join)
                else:
                    f.write(line)
            else:
                f.write(line)


def describe_norsum_header(norsum_file,t0,verbose=False,maxE=10000):
    """
    Read the header of a NORSUM file (from RPIXDR) and return a
    description including the crunch boundaries and factors for 
    SAMMY input file.

    Parameters
    ----------
    norsum_file : ascii file
        An ascii file that has the format of NORSUM files that are
        produced by the program RPIXDR.
    t0 : float
        The time zero [us] for the time of flight data give in 
        norsum_file. This along with the flight path that resides
        in norsum_file calculates the energy crunch boundaries.

    Returns 
    -------
    no return : print statement
        A print statement, that includes crunch boundaries and 
        factors that can be copy and pasted into a SAMMY inp 
        file.

    """
    if maxE == 10000:
        print("----\nDefault value of 10 keV for max energy of ")
        print("crunch boundary is being used\n----")
    header_strings = open(norsum_file).readlines()[0:14]
    if verbose == True:
        for string in header_strings:
            print(string.strip())
    print(header_strings[2].strip())
    # collect data
    fp = float(header_strings[11].strip().split()[0])
    bw = float(header_strings[12].strip().split()[0])
    trig = float(header_strings[13].strip().split()[0])
    ncpts = int(header_strings[5].strip().split()[0])
    cpts,cfactors = np.zeros(4),np.zeros(4)
    tofcpts,ecpts = np.zeros(4),np.zeros(4)
    grp_bw = np.zeros(4)
    for i in range(len(cpts)):
        temp_array = header_strings[6+i].strip().split()
        cpts[i],cfactors[i] = float(temp_array[0]),float(temp_array[1])
        print(cfactors[i])
        
        if i==0:
            grp_bw[i] = bw
            tofcpts[i] = cpts[i]*grp_bw[i]
        else:
            grp_bw[i] = grp_bw[i-1]*2**cfactors[i-1]
            tofcpts[i] = (cpts[i]-cpts[i-1])*grp_bw[i]+tofcpts[i-1]
        ecpts[i] = tt.tofe(tofcpts[i]-t0,fp)
    print('FP :         ',fp)
    print('bin width :  ',bw)
    print('triggers :   ',trig)
    if verbose == True:
        print('comp energies : ',ecpts,'\n','comp factors :  ',cfactors)
    print('--------\nSAMMY input: \n--------\n{:9.8f}{:5d}'.format(bw,ncpts))
    print('{:10.2f}{:10.1f}{:10.2f}{:10.1f}{: >10.2f}{:10.1f}{:10.2f}{:10.1f}'.format(ecpts[2],
          grp_bw[3]/bw,ecpts[1],grp_bw[2]/bw,ecpts[0],grp_bw[1]/bw,maxE,grp_bw[0]/bw))

# input: numPars- how many parameters were varied by SAMMY
#        infile- SAMMY.LPT file that includes printed 
#        covariance matrix.
# output: returns string of new covariance matrix file, and 
#         creates new cov. mat. file.
# take LPT file from sammy, that includes covariance matrix
# and return csv file cov matrix
def LPTtoCOV(numPars,infile):
    """
    Uses SAMMY output LPT file to create a covariance matrix
    file.

    SAMMY.LPT file contains correlations and standard deviations
    for the parameters that were allowed to vary. These are used
    to create a covariance matrix in a more usable text file.

    Parameters
    ----------
    numPars : int
        The number of parameters that were allowed to vary in the 
        SAMMY calcualation.
    infile : str
        The file address of the LPT file to be interrogated.

    Returns
    -------
    newfile : str
        File address of the new csv covariance file

    Notes
    -----
    Function is relatively primitive, it will only work with 
    numPars > 15. 
    """
    
    # !! Requires parameters in LPT file to be >15 !!
    if numPars<15:
        print("You have to edit the source to search for")
        print("strings of \"2\" or \"3\" or whichever numbers")
        print("describe the columns of the correlation matrix")
        print("in the LPT file.")
        raise ValueError("Number of parameters < 15.")

    # initialize correlation matrix of zeros
    cormat = np.zeros((numPars,numPars))
    # initialize standard deviation matrix of zeros
    stdmat = np.zeros((numPars,numPars))
    
    
    flag1,flag2,flag3 = False,False,False
    startval,i = 0,0                                 # used to track column position in cor. matrix
    lines = []
    with open(infile) as lpt:
        lines = [line.split() for line in lpt]       # split every line into str's by space
    for line in lines:
        #print(line)
        if "Number" in line and "varied" in line and "parameters" in line:
            number_of_param = int(line[len(line)-1]) # check number of parameters stated in LPT
            print(' '.join(line))
            if number_of_param != numPars:
                raise ValueError("Re-run program with {} parameters".format(number_of_param))
        if flag1 and "MISSING" in line and "LINES" in line:   # if keywords break the loop
            print("Finished.")
            break
        if flag2 and len(line)==0:                   # if end of block of covar's; start new col's
            startval += 15
            flag3 = True
            if i+1 != numPars:
                print("!!!!\nIs the specified number of par's correct?")
                print("Number of parameters found",i+1)
                print("The matrix may be too large!\n!!!!")
            continue
        if flag3 == True:
            flag3 = False                            # need to skip an extra line, continue again
            continue
        if "CORRELATION" and "MATRIX" and  "OUTPUT" in line:          # starts on the output matrix
            print("found it")
            flag1 = True
        if "1" in line and "14" in line and "15" in line and flag1:   # starts on col. definitions
            print("Starting build...")
            flag2 = True
            continue
        if flag2:
            #print(line)
            ll = len(line)
            i = int(line[0])-1                       # i'th row in cormat is line[0]-1=i
            if i>numPars:
                break                                # stop data collection if numPars is exceeded
            for j in range(3,ll):
                cormat[i][startval+j-3] = float(line[j])/100.0 # set proper elem. of cormat to float
            stdmat[i][i] = float(line[1])                      # set proper elem. of stdmat to float
    
    # correlation matrix, so diagonal values=1
    for k in range(numPars):
        cormat[k][k] = 1.0
        
    covmat = np.dot(np.dot(stdmat,cormat),stdmat)    # [cov] = [std]*[cor]*[std]
    print("Matrix dim.: ",np.shape(covmat))
                
    f = open("covmat_{}.dat".format(infile),"w+")
    for line in covmat:
        printline = [str(val) for val in line]
        f.write(",".join(printline)+"\n")

    f.close()
    print("Written in csv to: covmat_{}.dat".format(infile))

    newfile = "covmat_{}.dat".format(infile)

    return newfile                                   # return string filename 
#------------------------------------------------------------------------------
def mini_cov(infile,sel_pars):
    """
    Convert full covariance matrix file to one only made of 
    parameters you choose.

    Takes covariance matrix file generated by LPTtoCOV() and
    converts large matrix into matrix of covariances specific 
    to the user selection so as to reduce the complexity of 
    the matrix. (300x300 matrix --> 10x10 matrix)

    Parameters
    ----------
    infile : text file
        File of covariances, should be formatted the same as
        the output from LPTtoCOV().
    sel_pars : array
        Vector containing the parameter integer ID numbers
        assigned by sammy and listed in the sammy.lpt file 
        that are specific to the varied parameters the user
        is interested in.
    
    Returns
    -------
    select_cov : numpy array
        2-d matrix containing the covariances of the 
        parameters specified by sel_pars. 
    
    """
    
    # check if input filename contains expected string
    if "covmat" not in infile:
        print("Warning, this is not standard output from LPTtoCOV() function.")
    
    # read them into a lower-triang. array
    cov = np.array(pd.read_csv(infile,header=None))
    print("Initial array size: ",np.shape(cov))

    # convert from lower-tri matrix to symmetric full matrix by adding the transpose
    # and subtracting the diagonal elements
    cov = cov + cov.T - np.diag(cov.diagonal())

    # Test each element on whether cov = transposed cov. If symmetric, all elements should
    # be equal, if not, the boolean array will contain a False value
    if f.symmetric2dTest(cov) == False:
        print("!\n!\n!\n!\n!\n!\n!\n")
        print("Warning, information may have been lost during symmetrization of cov. matrix.")
        print("!\n!\n!\n!\n!\n!\n!\n")
    
    # check to see if user has selected a parameter not in matrix
    if True in [np.shape(cov)[0]<val or val<0 for val in sel_pars]:
        raise ValueError("The cov. mat. doesn't contain some of the requested ID's in sel_pars.")
    
    # make a new array of selected parameters
    select_cov = np.array([[cov[row-1][col-1] for row in sel_pars] for col in sel_pars])
    print("New covariance array of {} parameters.".format(np.shape(select_cov)))
    
    return select_cov
#------------------------------------------------------------------------------