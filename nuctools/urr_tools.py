import os,sys,time
import pandas as pd
import numpy as np
import subprocess as sub
import yaml
from pathlib import Path
from . import funky as fk 
from . import sam_tools as st

__all__ = ["sesh_fitacs","write_pars_fitacs","apply_sesh_corr",
           "write_sesh_pars","read_fitacs_par","create_new_sesh_mcfile"]

def sesh_fitacs(inp_file_name,data):

    """
    This function performs the iteration loop for SESH and FITACS

    !! This is deprecated since A. Golas implemented SESH into SAMMY !!

    SESH and FITACS are used together to correct experimental 
    cross section data and fit that data, respectively. Both 
    programs require average resonance parameters and to use them
    in concert requires that the input to both be nearly the 
    same. Therefore an iteration must be performed feeding the 
    input and output between the SESH and FITACS programs. The
    sesh_fitacs() function performs this.

    Parameters
    ----------
    inp_file_name : str
        The input file must be a yaml format file (.yml or .yaml), with 
        variables defined as in the example file.
    data : list
        A list of pandas DataFrames, which contain very specific column
        names. This list of DataFrames contain the experimental data 
        used to fit in FITACS. The order of this list must be the same
        order as the lists given in the input file inp_file_name.

    Returns
    -------
    none : ascii files
        This function will generate a lot of ascii files output from 
        SESH and FITACS. These will have naming conventions that are
        used in the function. The final fit with FITACS is after 
        the convergence criterion has been met and completes the 
        while loop.

    Examples
    --------
    >>> xs_col = ['e','cs','dcs']
    >>> tr_col = ['e','t','dt']
    >>> totxs_ta1 = pd.read_csv("ta1_sig.dat"  ,names=xs_col,sep=r'\s+')
    >>> trans_ta1 = pd.read_csv("ta1_trans.dat",names=tr_col,sep=r'\s+')
    >>> data1 = pd.concat([totxs_ta1,trans_ta1[['t','dt']]],axis=1)
    >>>
    >>> totxs_ta3 = pd.read_csv("ta3_sig.dat"  ,names=xs_col,sep=r'\s+')
    >>> trans_ta3 = pd.read_csv("ta3_trans.dat",names=tr_col,sep=r'\s+')
    >>> data3 = pd.concat([totxs_ta3,trans_ta3[['t','dt']]],axis=1)
    >>>
    >>> totxs_ta6 = pd.read_csv("ta6_sig.dat"  ,names=xs_col,sep=r'\s+')
    >>> trans_ta6 = pd.read_csv("ta6_trans.dat",names=tr_col,sep=r'\s+')
    >>> data6 = pd.concat([totxs_ta6,trans_ta6[['t','dt']]],axis=1)
    >>>
    >>> capxs_ta1 = pd.read_csv("capxs_ta1.dat",names=xs_col,sep=r'\s+')
    >>> capxs_ta2 = pd.read_csv("capxs_ta2.dat",names=xs_col,sep=r'\s+')
    >>>
    >>> data = [data1,data3,data6,capxs_ta1,capxs_ta2]
    >>> 
    >>> directory = '/Users/jesse/Dropbox/ta_urr_fitting/'
    >>> inp_file = directory+'sesh_fitacs_inp.yml'
    >>>
    >>> nuc.sesh_fitacs(inp_file,data)

    Notes
    -----
    This is to be more fully documented in the greater nuctools documentation.

    """
    
    start_tot = time.time()
    
    # ------------------------------------------------------------------------
    # import all the variables from the input file into a dict
    # ------------------------------------------------------------------------
    with open(inp_file_name,'r+') as f:
        inp = yaml.load(f)
    
    wdir          = inp['workdir']
    fitacs_indir  = wdir+inp['fitacs_indir']
    fitacs_outdir = wdir+inp['fitacs_outdir']
    sesh_indir    = wdir+inp['sesh_indir']
    sesh_outdir   = wdir+inp['sesh_outdir']
    ebd           = inp['e_bounds']
    thickness     = inp['samp_thick']
    cap_bool      = inp['cap_bool']
    data_names    = inp['data_names']
    data_dir      = wdir+inp['data_dir']
    # ------------------------------------------------------------------------
    
    # ------------------------------------------------------------------------
    # iterate until fitted FITACS pars are within the bounds of experimental 
    # uncertainty of the SESH pars used to correct the data input to FITACS
    # ------------------------------------------------------------------------
    #
    # ----------------------
    # set par filename as the starting FITACS par file
    # ----------------------
    temp_filename = fitacs_indir+inp['fitacs_par_name']
    
    # ----------------------
    # Calculate correction factor for first set of pars?
    # Often the beginning set has been calculated by a previous run
    # ----------------------
    calc_first_corr = inp['calc_first_corr']
    
    
    iteration = 0
    sesh_fitacs_pars_not_within_unc = True
    
    while sesh_fitacs_pars_not_within_unc:
        start_it = time.time()
        
        pars = read_fitacs_par(temp_filename,num_e_regions=3)
    
        print('#######################################')
        print('Iteration ',iteration)
        print('#######################################')
        
        # convergence boolean is tested in loop
        if iteration > 0:
            sesh_fitacs_pars_not_within_unc = False

        for e in range(inp['fitacs_numE_regions']):
            start_e = time.time()
            print('----------------')
            print('Energy region',e)
            print('----------------')
            # ----------------------
            # create a sesh file for each E_region from the pars read from fitacs
            # ----------------------
            temp_sesh_filename = sesh_indir+inp['sesh_ifile_name']
            write_sesh_pars(temp_sesh_filename,
                            temp_sesh_filename+'_e{}_it{}'.format(e,iteration),
                            pars[e*3:e*3+3].reset_index(),inp['Rp'])

            # ----------------------
            # write a new interactive input file for each e region
            # ----------------------
            
            temp_int_input = sesh_indir+inp['sesh_int_input']
            with open(temp_int_input,'r+') as f:
                with open(temp_int_input+'{}'.format(e),'w+') as f2:
                    for line in f:
                        # remove quotes, trailing whitespace, and newlines
                        line = line.replace('\"','').strip()
                        f2.write('\"'+line+'_e{}_it{}\"\n'.format(e,iteration))

            # ----------------------
            # run the SESH program
            # ----------------------
            if iteration == 0 and calc_first_corr == False:
                calculate = False
                print("Skipping first SESH calc, as it already exists.")
                print("This can be controlled by input.")
            else:
                calculate = True
            if calculate:
                # -> create a process
                run_sesh = sub.Popen('sesh',stdin=sub.PIPE,stdout=sub.PIPE,stderr=sub.PIPE)
                # -> make a byte string to send to process
                with open(temp_int_input+'{}'.format(e),'r+') as f:
                    byte_str_inp = str.encode(f.read())
                # -> send the byte string to SESH and run
                print("SESH MC calculation running...",time.ctime())
                out, err = run_sesh.communicate(input=byte_str_inp)
                #print(err)
                if (Path(wdir+'fort.99').exists()):
                    Path(wdir+'fort.99').unlink()
                print("SESH run completed.")

            # ----------------------
            # correction file must be interpreted from the results in output
            # ----------------------
            temp_corr_file = sesh_outdir+inp['sesh_cor']+'_e{}_it{}'.format(e,iteration)
            temp_output = sesh_outdir+inp['sesh_output']+'_e{}_it{}'.format(e,iteration)
            create_new_sesh_mcfile(temp_output,temp_corr_file)
            
            # ----------------------
            # Calculate difference between current and previous correction factors
            # ----------------------
            if iteration > 0:
                prev_corr_file = sesh_outdir+inp['sesh_cor']+'_e{}_it{}'.format(e,iteration-1)
                colnames = ['e','avtot','avcap','sscf','davtot','davcap','dsscf','thick','samp_type']
                current = pd.read_csv(temp_corr_file,sep=r'\s+',
                                      names=colnames,header=0)
                previous = pd.read_csv(prev_corr_file,sep=r'\s+',
                                       names=colnames,header=0)
                
                mean_change = np.mean(abs(current.sscf-previous.sscf)/previous.sscf)
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                print(fk.four(mean_change*100),"% change from previous correction factor")
                print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
                
                if mean_change > 0.01:
                    sesh_fitacs_pars_not_within_unc = True

            # ----------------------
            # apply correction to specific energy region of experimental data for 
            # each sample of thickness i,thick
            # ----------------------
            for i,thick in enumerate(thickness):
                # if capture sample
                if cap_bool[i]:
                    new_xs,new_dxs = apply_sesh_corr(temp_corr_file,
                                                     [data[i].e,data[i].cs,data[i].dcs],
                                                     thick,capture=True)
                # if transmission sample
                else:
                    new_xs,new_dxs = apply_sesh_corr(temp_corr_file,
                                                     [data[i].e,data[i].t,data[i].dt],
                                                     thick,capture=False)

                corr = pd.DataFrame({
                    'e': data[i].e,
                    'cs': new_xs,
                    'dcs': new_dxs,
                })
                #print(corr.cs)
                corr = corr[(corr.e > ebd[e]) & (corr.e < ebd[e+1])]
                # set data[thickness i] to corrected data
                data[i].cs[(data[i].e > ebd[e]) & (data[i].e < ebd[e+1])] = corr.cs
                data[i].dcs[(data[i].e > ebd[e]) & (data[i].e < ebd[e+1])] = corr.dcs

            finish_e = time.time()
            print("Time for E region {} = {} sec".format(e,fk.two(finish_e-start_e)))

        
        
        if sesh_fitacs_pars_not_within_unc is False:
            print("\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\")
            print("/////////////////////////////////////////////////////////////")
            print("SESH correction factor changed {:10.4f}% between iterations.".format(mean_change*100))
            print("Convergence achieved.")
            print("Final fit with FITACS...")
            print("\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\")
            print("/////////////////////////////////////////////////////////////")
        
        # ----------------------
        # Create new data files to fit w/ FITACS
        # ----------------------
        for i,thick in enumerate(thickness):
            filename = data_dir+data_names[i]
            st.fmt_ten([data[i].e,data[i].cs,data[i].dcs],
                        filename,tot= not cap_bool[i])

        # ----------------------
        # write a new interactive input file for FITACS
        # ----------------------
        temp_int_input = fitacs_indir+inp['fitacs_int_input']
        with open(temp_int_input,'r+') as f:
            with open(temp_int_input+'{}'.format('_'),'w+') as f2:
                for i,line in enumerate(f):
                    filename = data_dir+data_names[i-2]
                    if i < 2:
                        f2.write(line)
                    else:
                        f2.write(filename+'\n')


        # ----------------------
        # Fit corrected data with FITACS
        # ----------------------
        # -> create a process
        run_fitacs = sub.Popen('sammy',stdin=sub.PIPE,stdout=sub.PIPE,stderr=sub.PIPE)
        # -> make a byte string to send to process
        with open(temp_int_input+'{}'.format('_'),'r+') as f:
            byte_str_inp = str.encode(f.read())
        # -> send the byte string to SESH and run
        out, err = run_fitacs.communicate(input=byte_str_inp)
        #print(err)
        print('\n+++++++++++++++++++++')
        print("FITACS run completed.")
        print('+++++++++++++++++++++')

        # ----------------------
        # Rename FITACS output and delete unnecessary files
        # ----------------------
        fit_out = inp['fitacs_outdir']
        p = Path(wdir)
        for item in list(p.glob('*.lst')):
            item.rename(Path(item.parent, "{}{}".format(fit_out,
                                                        item.stem+item.suffix)))
        # .unlink() deletes the files 
        for item in list(p.glob('*.odf')):
            item.unlink()
        for item in list(p.glob('*.plt')):
            item.unlink()
        for item in list(p.glob('*.DAT')):
            item.unlink()
        Path(wdir+'SAMMY.IO').unlink()

        # rename files
        par_path = Path(wdir+'SAMMY.PAR')
        lpt_path = Path(wdir+'SAMMY.LPT')
        cov_path = Path(wdir+'SAMMY.COV')
        par_path_targ = Path(wdir+fit_out+'fit.par')
        lpt_path_targ = Path(wdir+fit_out+'fit.lpt')
        cov_path_targ = Path(wdir+fit_out+'fit.cov')
        par_path.rename(par_path_targ)
        lpt_path.rename(lpt_path_targ)
        cov_path.rename(cov_path_targ)
        
        # failsafe: only allow 4 iterations
        if iteration > 3:
            print("\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\")
            print("///////////////////////////////////////")
            print("Max iterations reached, no convergence.")
            print("\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\")
            print("///////////////////////////////////////")
            break
        
        # ----------------------
        # If another iteration is necessary, create a new par file
        # ----------------------
        if sesh_fitacs_pars_not_within_unc:
            # ----------------------
            # read the fitted parameters 
            # ----------------------
            pars_new = read_fitacs_par(wdir+fit_out+'fit.par',num_e_regions=3)
            
            p = Path(temp_filename)
            # create new FITACS input par file with new pars
            new_filename = (str(p.parent)+"/"+str(p.stem)+
                            "_it{}".format(iteration)+str(p.suffix))
            
            write_pars_fitacs(temp_filename,new_filename,pars_new)
            
            temp_filename = new_filename
            
        finish_it = time.time()
        it_time = fk.two(finish_it-start_it)
        print("Time for iteration {} = {}".format(iteration,it_time))
        iteration += 1
            
    
    finish_tot = time.time()
    tot_time = finish_tot-start_tot
    print("Total time = {}".format(tot_time))
    
    with open('run_{}_stats.txt'.format(time.ctime()),'w+') as f:
        f.write('Various run statistics:\n')
        f.write('-----------------------\n')
        f.write("Total time = {}\n".format(tot_time))
        f.write('Total number of iterations = {}\n'.format(iteration))

def write_pars_fitacs(init_filename,new_filename,pars):
    """
    Write average parameters from memory into a pre-formatted 
    FITACS par file.
    
    Parameters
    ----------
    init_filename : str
        The initial par file name (full or relative path). This
        should have the exact same number of parameters for each
        partial wave as the pars DataFrame input
    new_filename : str
        The name of the file to be written to (full or relative 
        path). This file will be the same as init_filename but 
        with different parameters
    pars : pandas DataFrame
        The pars DataFrame has average resonance parameters 
        necessary to run the SESH and FITACS codes. This variable
        has columns of data described in read_fitacs_par()
    
    Returns
    -------
    none : ascii file
        Replaces the parameters in init_filename and re-writes the
        file to new_filename.
    
    """
    
    def ten(number):
        return '{:10.6f}'.format(number)
    
    with open(init_filename,'r+') as f: file = f.readlines()
    with open(new_filename,'w+') as f:
        k = -1
        # skip this many lines before writing to new file
        write_line = 0
        for line in file:
            # found the average resonance parameters
            if ('stre' in line) or ('STRE' in line):
                f.write(line)
                k+=1
                # for all ang mom L
                for j in range(int(max(pars.L)+1)):
                    f.write((ten(pars.stren[k+j])+ten(pars.dstren[k+j])+
                      ten(pars.dist_R[k+j])+ten(pars.ddist_R[k+j])+
                      ten(pars.gam_g[k+j])+ten(pars.dgam_g[k+j])))
                    if j==0:
                        f.write(ten(pars.D[0])+'\n')
                    else:
                        f.write('\n')
                write_line = int(max(pars.L)+2)
            if write_line > 0:
                write_line-=1
                continue
            else:
                f.write(line)



def apply_sesh_corr(cor_file,data,N,capture=True):
    """
    Correct experimental data with a Monte Carlo correction from SESH

    Parameters
    ----------
    cor_file : str
        The full file path or relative to the correction factor file
    data : array-like
        2-d array with data[0] as energy; For capture, data[1] as cross section,
        data[2] as the absolute uncertainty in the cross section. If correcting 
        total cross section data[1] should be the observable (transmission), 
        data[2] is the absolute uncert. on the observable. Units for xs are
        barns.
    N : float
        The areal number density [at/barn] of the sample is needed if correcting
        the total cross section, but not for the capture cross section
    capture : bool, optional
        Whether the correction is being applied to capture or total cross
        section. If total, N must be specified

    Returns
    -------
    new_xs,new_dxs : tuple
        A tuple of numpy arrays containing the corrected cross section [barns]
        and the uncertainty on the corrected cross section [barns].

    Examples
    --------
    >>>
    
    Notes
    -----
    Capture

    .. math:: \\alpha = correction
    
    .. math:: \\sigma_{\\gamma,\\alpha} = \\sigma_{\\gamma}/\\alpha
    
    .. math:: \\Delta\\sigma_{\\gamma,\\alpha} = \\sigma_{\\gamma,\\alpha}\\sqrt{\\frac{\\Delta\\sigma}{\\sigma}^2+\\frac{\\Delta \\alpha}{\\alpha}^2}
    
    Trans
    
    .. math:: \\sigma_{t,\\alpha} = \\frac{-1}{N}ln\\left(\\frac{T}{\\alpha}\\right)
    
    .. math:: \\Delta\\sigma_{t,\\alpha} = \\frac{1}{N}\\sqrt{\\frac{\\Delta T}{T}^2+\\frac{\\Delta\\alpha}{\\alpha}^2}
    
    
    """
    # the column names 
    colnames = ['e','totxs','dtotxs','capxs','dcapxs','sscf',
                'dsscf','trans','dtrans','thick','samp_type']

    if (N is None) and (capture is False):
        raise ValueError("If correcting total xs, N must be specified.")
    
    # cast data to numpy 
    data = np.array(data)
    # read in the correction factor file
    sesh = pd.read_csv(cor_file,sep=r'\s+',names=colnames,
                       header=0,float_precision='round_trip') 
    sesh = sesh[sesh.thick == N]
    
    # interpolate the correction onto the data energy grid
    cor = np.interp(data[0],sesh.e*1000,sesh.sscf)
    dcor = np.interp(data[0],sesh.e*1000,sesh.dsscf)/100 # percent
    
    # apply the correction
    if capture:
        new_xs = data[1]/cor
        new_dxs = new_xs*np.sqrt((data[2]/data[1])**2+(dcor)**2)
    else:
        new_xs = -1/N*np.log(data[1]/cor)
        new_dxs = 1/N*np.sqrt((data[2]/data[1])**2 + (dcor)**2)
        
    return new_xs,new_dxs

def write_sesh_pars(ifilename,ofilename,pars,Rp,max_orb_ang_mom=2):
    """
    Re-write the SESH input file with new resonance parameters.

    Parameters
    ----------
    ifilename : str
        The base file you want to re-write. Full path or relative.
    ofilename : str
        The new input file with updated parameters. Full path or relative
    pars : Pandas DataFrame
        DataFrame that is returned from the function read_fitacs_par() which
        is described below. DataFrame attributes (columns) must have the same
        names as described in read_fitacs_par()
    Rp : float
        The effective radius which will be used for the SESH calculation
    max_orb_ang_mom : int, optional
        The maximum orbital quantum angular momentum. SESH accepts avg 
        parameter values for each L value

    Returns
    -------
    None : ascii file
        Returns a re-written SESH ascii input file with new average parameters

    Examples
    --------

    >>> import nuctools as nuc
    >>> directory = '/Users/jesse/Dropbox/ta_urr_fitting/'
    >>> pars = nuc.read_fitacs_par(directory+'multi_fit.par',num_e_regions=3)
    >>> nuc.write_sesh_pars(directory+'input.data',directory+'input.data1',pars[0:3],7.8)

    """
    with open(ifilename,'r+') as f:
        # create a list of line strings
        file = f.readlines()
    with open(ofilename,'w+') as f:
        i = 0
        for line in file:
            # write the pars in file lines 3->(3+max orb ang momentum)
            if (i > 1) and (i <= 2+max_orb_ang_mom):
                # write the avg pars for each partial wave
                f.write(fk.ten(pars.gam_g[i-2])+' '
                        +fk.ten(pars.D[i-2])+' '
                        +fk.ten_exp(pars.stren[i-2])+' '
                        +fk.ten(0.0)+fk.ten(Rp)+fk.ten(1.0)+'\n')
                i += 1
                continue
            # re-write the lines for everywhere else
            f.write(line)
            i += 1


def read_fitacs_par(filename,num_e_regions=1,max_orb_ang_mom=2):
    """
    Function to read the average parameters from a FITACS/SAMMY par file

    The output par file from FITACS has average resonance parameters for
    multiple angular momentums, and multiple energy regions (if user 
    specified). This function reads the file based on the user input 
    for those two values.

    Parameters
    ----------
    filename : str
        The name of the file (full path or relative) to be read
    num_e_regions : int, optional
        The number of energy regions over which SAMMY/FITACS is fitting
        the average resonance parameters
    max_orb_ang_mom : int, optional
        The maximum orbital angular momentum that was considered when 
        running the SAMMY/FITACS program. A typical value is L=2.

    Returns
    -------
    pars : Pandas DataFrame
        A dataframe that contains the following attributes: 

        - L : the orbital angular momentum
        - D : the average level spacing
        - eregion : the number of energy regions fit
        - stren : the neutron strength function
        - dstren : the unc. on the neutron strenth function
        - dist_R : the distant level parameter R-inf
        - ddist_R : the unc. on dist_R
        - gam_g : the average radiation width 
        - dgam_g : the unc. on gam_g

    Examples
    --------

    >>> import nuctools as nuc
    >>> pars = nuc.read_fitacs_par('multi_fit.par',num_e_regions=3)
    >>> print(pars)

    """
    with open(filename,'r+') as f:
        res_pars_next = False
        pars_break = False
        df_len = num_e_regions*(max_orb_ang_mom+1)
        pars = pd.DataFrame({
            'L'      : np.zeros(df_len),
            'D'      : np.zeros(df_len),
            'eregion': np.zeros(df_len),
            'stren'  : np.zeros(df_len),
            'dstren' : np.zeros(df_len),
            'dist_R' : np.zeros(df_len),
            'ddist_R': np.zeros(df_len),
            'gam_g'  : np.zeros(df_len),
            'dgam_g' : np.zeros(df_len),
        })
        region = 0
        L = 0
        
        for i,line in enumerate(f):
            #print(i)
            if len(line) < 4:
                continue
                
            # remove the line endings and lines with only \n
            line = line.strip('\n')
            
            # test if resonance parameters are coming next
            if ('STRE' in line) or ('Stre' in line):
                #print("region = ",region,"\n---------")
                res_pars_next = True
                pars_break = False
                L = 0
                continue
            # test if we've reached the end of the res. pars 
            if '---' in line:
                pars_break = True
            if res_pars_next and not pars_break:
                # grab the level spacing for s-waves
                if L == 0:
                    D = float(line[60:70])
                # choose element to fill
                pos = region*3 + L
                # set the values
                pars.at[pos,'L']       = L
                pars.at[pos,'D']       = D
                pars.at[pos,'eregion'] = region
                pars.at[pos,'stren']   = float(line[0:10])
                pars.at[pos,'dstren']  = float(line[10:20])
                pars.at[pos,'dist_R']  = float(line[20:30])
                pars.at[pos,'ddist_R'] = float(line[30:40])
                pars.at[pos,'gam_g']   = float(line[40:50])
                pars.at[pos,'dgam_g']  = float(line[50:60])
                
                # increment the orbital ang momentum value
                L += 1
            if res_pars_next and pars_break:
                res_pars_next = False
                # increment the energy region
                region += 1
                
    return pars


def create_new_sesh_mcfile(output_filename,clean_mc_filename,runtype='trans'):
    """
    Read the sesh output file defined by output_filename. Collect
    Monte Carlo run data including average capture cross section,
    average total cross section, self-shielding correction factor,
    and the error on all of these.

    Parameters
    ----------
    output_filename : str
        the name of the output file from a sesh run. Usually is
        output.data, but I suggest renaming it as sesh will 
        keep overwriting this file each time you use sesh.
    clean_mc_filename : str
        The name of the file you wish to write the MC info to.    
    runtype : str, optional
        The type of SESH run that was performed. Right now the 
        supported types are transmission (input must include 
        either "trans" or "TRANS" in string) and capture 
        (runtype is left at default None). If both transmission
        and capture are input to the SESH program, this function
        will decide where there is capture or transmission data
        being given in the output file and this input is 
        unnecessary.

    Returns
    -------
    None : ascii file
        Returns an ascii file with a descriptive header for capture.
    """
    # Choose the type of self-shielding that was calculated in SESH
    captype = True
    if 'trans' in runtype or "TRANS" in runtype:
        captype = False

    # open new output sesh file and collect data
    with open(clean_mc_filename,"w+") as fnew:
        # write a header
        fnew.write("Energy[keV] Av.Tot.xs[b] dAv.Tot.xs[b] Av.Cap.xs[b] dAv.Cap.xs[b]")
        fnew.write(" SSCor.Fac.[] dSSCor.Fac.[%] trans[] dtrans[] ")
        fnew.write("thick[at/b] samp_type\n")
        # open the file from sesh
        with open(output_filename,"r+") as f:
            # create array of line strings from sesh file
            f_string_array = f.readlines()
            start_of_mc = False
            type_string = ''
            for line in f_string_array:
                line_long_enough = False
                thick_in_line = False
                if len(line) > 1:
                    line_long_enough = True
                if "TRANSMISSION" in line:
                    captype = False
                    type_string = 'trans'
                if "CAPTURE DATA" in line:
                    captype = True
                    type_string = 'capt'
                if "(NUC./B)" in line:
                    start_of_mc = True
                    continue
                if "M O N T E" in line:
                    start_of_mc = False
                if line_long_enough and start_of_mc and ('.' in line[1:10]):
                    thick_in_line = True
                    thickness = line[0:10]
                if line_long_enough and start_of_mc and ('I N P U T' in line):
                    # break if the end of the MC region is found
                    break
                
                if line_long_enough and start_of_mc:
                    # energy line -- begin new_line to write 
                    if "." == line[13]:
                        # initialize string variables
                        E_n,av_cap_xs,av_tot_xs,ss_cor_fc,t = '','','','',''
                        new_line = ''
                        E_n = line[6:19]
                        # sample thickness takes up extra columns, so read this one differently
                        if thick_in_line:
                            E_n = '    '+line[10:19]
                        av_cap_xs,av_tot_xs = line[19:28],line[29:38]
                        if captype:
                            ss_cor_fc = line[65:74]
                            t = line[119:128]
                        else:
                            ss_cor_fc = line[51:61]
                            t = line[39:48]

                    # collect the uncertainty information from next line and place in new_line
                    if "." != line[13]:
                        # initialize string variables
                        dav_cap_xs,dav_tot_xs,dss_cor_fc,dt = '','','',''
                        dav_cap_xs,dav_tot_xs = line[19:28],line[29:38]
                        if captype:
                            dss_cor_fc = line[65:74]
                            dt = line[119:128]
                        else:
                            dss_cor_fc = line[51:61]
                            dt = line[39:48]
                    
                    
                    if "." != line[13]:
                        new_line += E_n+' '+av_tot_xs+' '+dav_tot_xs+' '+av_cap_xs+' '+dav_cap_xs
                        new_line += ' '+ss_cor_fc+' '+dss_cor_fc
                        new_line += ' '+t+' '+dt+' '+thickness+' '+type_string
                    # write the new_line into the new clean file when the 1st and 2nd line
                    # for the energy has been read
                    if "." != line[13]:
                        fnew.write(new_line+"\n")


