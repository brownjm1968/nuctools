
import numpy as np
import pandas as pd

__all__ = ["TAB1",'read_pendf_xs','write_pendf_xs',"write_endf_float","read_endf_float",
           "update_file2","setupdict","update_pardict","strip_line_num",
           "loglog_interp"]

class TAB1:
    """
    A class to mimic an ENDF-6 format TAB1 record

    Parameters
    ----------
    filename : str
        The file path to the ENDF file
    mat : int
        The ENDF format MAT number
    mf : int
        The ENDF format MF number
    mt : int
        The ENDF format MT number
    starting_line_no : int, optional
        The starting line number of the MAT, MF, and MT

    Attributes
    ----------
    interp_table : list
        The interpolation table; a list of 2-element lists with values of number 
        of points and interpolation scheme for those points. The sum of number
        of points should be equal to attribute `np`
    c1 : float
        Control record number 1
    c2 : float
        Control record number 2
    l1 : int
        integer record number 1
    l2 : int
        integer record number 1
    nr : int
        number of interpolation ranges
    np : int
        number of tabulated pairs
    data : DataFrame
        A DataFrame with attributes `e` and `cs` for pointwise cross section

    Methods
    -------
    read_tab1_cr
        Read the control record and interpolation tabel for a TAB1

    """
    def __init__(self,filename,mat,mf,mt,starting_line_no=0):
        self.mat = mat
        self.mf = mf
        self.mt = mt
        self.filename = filename
        self.starting_line_no = starting_line_no
        self.c1 = None
        self.c2 = None
        self.l1 = None
        self.l2 = None
        self.nr = None
        self.np = None
        self.interp_table = []

        self.read_tab1_cr()

        self.data = pd.DataFrame(read_pendf_xs(self.filename,self.starting_line_no,self.ending_line_no).T,
                                  columns=['e','cs'])

    def read_tab1_cr(self):
        """ designed for file 3 """
        start,numlines_nr,numlines_np = 0,0,0
        read_interp_table = False
        with open(self.filename,'r') as f:
            lines = f.readlines()
        for i,line in enumerate(lines):
            if i<self.starting_line_no:
                continue
            # find starting line
            if "{:4d}{:2d}{:3d}".format(self.mat,self.mf,self.mt) in line and start==0:
                start = i+1
                continue
            # read control record
            if i==start and i!=0:
                self.c1 = read_endf_float(line[0:11])
                self.c2 = read_endf_float(line[11:22])
                self.l1 = read_endf_float(line[22:33])
                self.l2 = read_endf_float(line[33:44])
                self.nr = read_endf_float(line[44:55])
                self.np = read_endf_float(line[55:66])
                numlines_nr = int(self.nr//3)
                if self.nr%3 > 0: 
                    numlines_nr +=1
                numlines_np = int(self.np//3)
                if self.np%3 > 0:
                    numlines_np +=1
                read_interp_table = True
            if read_interp_table == True and i<=start+numlines_nr and i > start:
                for j in range(3):
                    npi = read_endf_float(line[j*22:j*22+11])
                    inter = read_endf_float(line[j*22+11:j*22+22])
                    if npi!=0.0 and inter!=0.0:
                        self.interp_table.append([npi,inter])
            if self.nr is not None:
                if i> numlines_nr + start:
                    break
        self.starting_line_no = start+3
        self.ending_line_no = start+2+numlines_np


def read_pendf_xs(file,start,finish):
    """
    Read part of a PENDF format file.
    
    
    PENDF format files have sections of the file where the cross
    section is listed in 10 character exponential format as:
    
    energy1 cross_section1 energy2 cross_section2
    
    format. This function opens the file, grabs the specified lines
    and converts to floats and returns an energy array and a cross 
    section array.
    
    Parameters
    ----------
    file : str
        The full pathname for the file to be read, should be in PENDF
        format.
    start : int
        The line number to begin reading data
    finish : int 
        The line number to stop reading data

    Returns
    -------
	pointwise_cs : numpy array
	    An array that contains pointwise cross section: energies at
	    pointwise_cs[0], and the cross section values at 
	    pointwise_cs[1]
	
	Examples
	--------
	>>> import nuctools as nuc
	>>> file = '/Users/jesse/Desktop/temp.pendf'
	>>> start = 98
	>>> finish = 18984
	>>> energy,cross_section = nuc.read_3col_pendf(file,start,finish)

    
    """
    with open(file) as f:
        e = []
        cs = []

        break_outer = False

        for i,line in enumerate(f):
            # -------------------------------
            # Stop the loop once finish is reached
            # -------------------------------
            if i == finish:
                break
            if i >= start-1:
            	# -------------------------------
            	# Only include first 66 columns, split on space
            	# and convert to an array of strings
            	# -------------------------------
                word_len = 11
                word_start = 0
                for j in range(6):
                    word = line[word_start:word_start+11]

                    if( j%2 == 0 ):
                        # -------------------------------
                        # Grab the energies, convert to readable format
                        # -------------------------------
                        if( word == '           ' ):
                            break_outer = True
                            break # end of TAB1
                        e.append(read_endf_float(word))
                    else:
                        # -------------------------------
                        # Grab cross section, convert to readable format
                        # -------------------------------
                        cs.append(read_endf_float(word))
                    word_start+=word_len

            if( break_outer ):
                break # end of TAB1
    
    # -------------------------------
    # Convert to floats
    # -------------------------------
    e = np.array(e).astype(float)
    cs = np.array(cs).astype(float)

    # -------------------------------
    # Stack them into a numpy array
    # -------------------------------
    pointwise_cs = np.array([e,cs])
    
    return pointwise_cs


def write_pendf_xs(filename,energy,cs,mat_num,file_num,reaction_num):
    """
    Write energy and cross section into the PENDF format to be used 
    mostly for file 3 in the ENDF format.

    Parameters
    ----------
    filename : str
        The full file name of the file which is to be written
    energy : array-like
        The energies for each of the cross section points
    cs : array-like
        The cross section points to be written into the PENDF
        format
    mat_num : int
        The integer number (e.g. 7328 for Ta-181) that describes
        the material for which the data is written. This is 
        specified in the ENDF format 
    file_num : int
        The file number (most often to be 3) as specified by the
        ENDF format
    reaction_num : int
        The reaction number as specified by the ENDF format. For 
        example the total cross section reaction number is 1

    Returns
    -------
    no return : ascii file
        The text in a file in PENDF format that can be placed 
        into the full ENDF file of interest.

    """
    # ----------------------
    # Open new file
    # ----------------------
    with open(filename,'w+') as f:
        # ----------------------
        # Check if cross section or energy are NaN
        # ----------------------
        if (np.isnan(energy).any()) or (np.isnan(cs).any()):
            raise ValueError("Input energy or cross section contains NaN's")
        # ----------------------
        # Check if cross section is negative
        # ----------------------
        cs = np.array(cs)
        if cs[cs<0].any():
            raise ValueError("Input cross section is negative.")
        
        for i in range(len(energy)):
            # ----------------------
            # find exponent on power of 10
            # ----------------------
            exponent_e = int(np.floor(np.log10(energy[i])))
            if exponent_e < 0:
                sign_e = "-"
            else:
                sign_e = "+"

            if( cs[i] == 0.0 ):
                exponent_cs = 0
            else:
                exponent_cs = int(np.floor(np.log10(cs[i])))
            if exponent_cs < 0:
                sign_cs = "-"
            else:
                sign_cs = "+"
            # ----------------------
            # Write the energy and cross section points
            # ----------------------
            f.write(" {:0=7.6f}{}{} {:0=7.6f}{}{}".format(energy[i]/10**exponent_e,sign_e,abs(exponent_e),
                                                        cs[i]/10**exponent_cs,sign_cs,abs(exponent_cs)))
            # ----------------------
            # Tag the last line with ENDF info
            # ----------------------
            if (i+1)%3 == 0:
                f.write("{:d} {:d}  {:d}\n".format(mat_num,file_num,reaction_num))
            # ----------------------
            # Fill in space at the end
            # ----------------------
            if (i == len(energy)-1) and ((i+1)%3 != 0):
                how_many_gaps = 3-(i+1)%3
                spaces_per_gap = 22
                extra_space = " "*spaces_per_gap*how_many_gaps
                f.write(extra_space+"{:d} {:d}  {:d}\n".format(mat_num,file_num,reaction_num))


def write_endf_float(value):
    """
    Return the ENDF format string of a floating point number

    Parameters
    ----------
    value : float
        A single floating point value

    Returns
    -------
    valstring : str
        An ENDF format string of input par `value`
    """
    if(abs(value) < 1e-9 or abs(value) > 9.999e9):
        raise ValueError("value is too small or too big")
    valstring = "{:>13.6e}".format(value).replace('e','').replace('+0','+').replace('-0','-')
    # with AMPX written files we use "-0" instead of "+0" for some reason
    if( '+0' in valstring ):
        valstring = valstring.replace('+0','-0')
    return valstring

def read_endf_float(string):
    """
    Read the ENDF format string of a floating point number and return float

    Parameters
    ----------
    string : str
        A string for a single floating point value

    Returns
    -------
    value : float
        The value for the ENDF string float
    """
    if string.strip() == "":
        return 0.0
    elif "e" in string.lower():
        return float(string)
    elif "." in string:
        strsplit = string.split('.')
        return float(strsplit[0]+"."+strsplit[1].replace("-","e-").replace("+","e+"))
    else:
        return float(string.replace("-","e-").replace("+","e+"))

def update_file2(infile,outfile,energy_dict,mat):
    """
    Replace the res's in the infile with new values
    
    note: resonance must have unique energy

    Parameters
    ----------
    infile : str
        A file that includes ENDF-format File 2. This is the base 
        file that the outfile will match except updated values
    outfile : str
        File written from `infile` and updated resonance values
        from `energy_dict`
    energy_dict : dict
        A dictionary that uses keys based on the energies that can
        be identified **in the base file** `infile`
    mat : int
        The unique ENDF format MAT id

    Returns
    -------
    None : None
        Write a new file to `outfile`
    """

    with open(infile,'r+') as f:
        with open(outfile,'w+') as of:
            matfile = "{} {}".format(mat,2151)
            for line in f:
                #if(we're in file 2: res's)
                if(matfile in line):
                    for ekey in energy_dict:
                        #if( we find the res )
                        if ekey in line:
                            linelist = list(line)
                            for i,channel_pair in enumerate(energy_dict[ekey]):
                                channel_number = channel_pair[0]
                                channel_value = channel_pair[1]
                                start = 11*channel_number
                                end = 11+start
                                linelist[start:end] = endf_float_str(channel_value)
                            line = "".join(linelist)
                of.write(line)

def setupdict(parfile):
    """
    Set up a dictionary of varied parameters

    Only varied parameters are included in the dict
    
    Parameters
    ----------
    parfile : str
        The path to a par file that was used to run a SAMMY-like
        operation (normal, MC, etc.).

    Returns
    -------
    pardict : dict
        A dictionary with keys of res energies. Varied pars are
        added to a list for the appropriate res energy

    Notes
    -----
    ** Must match exact E-values of ENDF file that dict will be used to update **

    """
    pardict = {}
    with open(parfile,'r+') as f:
        for line in f:
            flags = line[56:65].split(' ')
            try:
                flags = [int(f) for f in flags]
            except:
                continue
            # if we found res pars
            if( all(flags) <= 3 ):
                # if any varied pars
                if( any(flags) > 0 ):
                    # energies are dict keys
                    estring = endf_float_str(float(line[0:11]))
                    pardict[estring] = []
                    pars = [float(line[0+11*i:11+11*i]) for i in range(len(flags))]
                    for i,flag in enumerate(flags):
                        if( flag > 0 ):
                            pardict[estring].append((i,pars[i]))
    return pardict

def update_pardict(pardict,mcpars):
    """
    Update the parameter dict using a pandas Series
    of values taken from a MC run 

    Parameters
    ----------
    pardict : dict
        A dict with keys of base res energies and values
        of varied parameters from a SAMMY operation
    mcpars : Series
        A Pandas Series that has pre-set key names from
        a MC SAMMY run (fitAPI)

    Returns
    -------
    none : None
        No return, only updates input parameter `pardict`

    """
    # figure out which pandas key matches dict keys
    colnames = mcpars.keys().tolist()
    keycounter = 0 # where in pardict key are we?
    prevnamenum = 0
    # run over MC names
    for i,name in enumerate(colnames):
        namenum = float(name.split("_")[1])
        if( namenum!= prevnamenum ):
            keycounter = 0
            prevnamenum = namenum
        for key in pardict:
            # get key-energy to see if it matches MC name
            mult = 1
            tempkey=key
            if key[0]=="-":
                tempkey = key[1:len(key)]
                mult = -1
            keynum = mult*float(tempkey.replace('-','e-').replace('+','e+'))
            percdiff = abs((namenum-keynum)/((namenum+keynum)/2))*100
            # found a matching res
            if( percdiff < 0.001 ):
                parname = name.split('_')[0]
                if parname == "Energy":
                    parInd = 0
                elif parname == "gamma":
                    parInd = 1
                elif parname == "IncChan":
                    parInd = 2
                else:
                    raise ValueError("par type: {} is not implemented".format(parname))
                pardict[key][keycounter] = (parInd,mcpars[name])
                keycounter+=1

def strip_line_num(inpfile,outfile):
    """
    Strip out the line numbers in an ENDF file, which occur after
    column number 75

    Parameters
    ----------
    inpfile : str
        The path and name of the input file
    outfile : str
        The path and name of the output file that will have line
        numbers stripped out
    """
    with open(inpfile,"r+") as f:
        inplines = f.readlines()
    with open(outfile,"w+") as f:
        for line in inplines:
            line = line[0:75] + '\n'
            f.write(line)

def loglog_interp(x,xp,fp):
    """
    Translate x,xp,fp into log-log space, linearly interpolate, and
    return them back in the same space as before

    Parameters
    ----------
    x : array-like
        The new grid onto which we want to interpolate
    xp : array-like
        The old grid from which we interpolate
    fp : array-like
        The "y-values" for xp

    Returns
    -------
    y : array-like
        The "y-values" that we linearly interpolated in log-log space
    """
    lx  = np.log(x)
    lxp = np.log(xp)
    lfp = np.log(fp)
    return np.exp(np.interp(lx,lxp,lfp))



