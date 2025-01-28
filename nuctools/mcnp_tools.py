import pandas as pd
import numpy as np
import warnings
import os,sys

__all__ = ['retrieve_tally','modify_sdef']

def retrieve_tally(infile,cells,numbins,skip_extra=0):
    """
    Find the line numbers where tally data starts in mcnp
    output files given the cell numbers being tallied

    Parameters
    ----------
    infile : ascii file
        The mncp output file
    cells : array like
        The cell or surface numbers that are being tallied, 
        if just one, give the single int in an array
    numbins : int
        The number of bins in the tally
    skip_extra : int, optional
        If there are additional comment lines before the 
        tally data, skip this many lines

    Returns
    -------
    ind : numpy array
        An array of the tally data
    """
    file = open(infile).readlines()
    ind = []
    for cell_number in cells:
        for i,line in enumerate(file):
            cell_str = "cell  {} ".format(cell_number)
            surf_str = "surface  {} ".format(cell_number)
            if cell_str in line or surf_str in line:
                ind.append(i+2)

    tally_list = []

    for row in ind:
        temp = np.genfromtxt(infile,skip_header=row+skip_extra,max_rows=numbins,unpack=True)
        tally_list.append(temp)

    return tally_list


def modify_sdef(source,mcnp_file_name,si_string,sp_string):
    """
    Modify the source information and source probability cards in 
    MCNP

    Parameters
    ----------
    source : array-like
        An array of source information and probability. source[0] is
        what will be defined as the "source information" and source[1]
        is what will be defined as the "source probability".
    mcnp_file_name : string
        The full file path including the name of the file
    si_string : str
        The string that exists in the MCNP input file you want to modify,
        this string will be used to identify the line to replace for the 
        source information (gamma energies) Example is \'SI4\'. Only 
        include the card string, it should be unique.
    sp_string : str
        The string that exists in the MCNP input file you want to modify,
        this string will be used to identify the line to replace for the 
        source probability (frequency of gamma energies occuring) Example 
        is \'SP4\'. Only include the card string, it should be unique.

    """
    # ---------------------------------------------------------
    # ----- Write a new input file with the spectrum
    # ---------------------------------------------------------
    found_si_string = False
    found_sp_string = False
    if len(source[0]) != len(source[1]):
        raise ValueError("source energies and probabilities not same length!")

    if '/' in mcnp_file_name:
        filesplit = mcnp_file_name.split('/')
    else:
        filesplit = mcnp_file_name.split('\\')
    filesplit[len(filesplit)-1] = 'mod_'+filesplit[len(filesplit)-1]
    newfilename = '/'.join(filesplit)

    source = np.array(source)

    frequencies = source[1]/source[1].sum()

    with open(mcnp_file_name,'r+') as f:
        start_file = f.readlines()
    with open(newfilename,'w+') as f:
        for line in start_file:
            if si_string in line:
                found_si_string = True
                f.write(si_string+' A')
                for i,info in enumerate(source[0]):
                    if i==0:
                        f.write('      ')
                    f.write('{:1.5e} '.format(info))
                    if i%6==0 and i!=0:
                        f.write('\n      ')
                    if i==len(source[0])-1:
                        f.write('\n')
            elif sp_string in line:
                found_sp_string = True
                f.write(sp_string)
                for i,frequency in enumerate(frequencies):
                    if i==0:
                        f.write('      ')
                    f.write('{:1.5e} '.format(frequency))
                    if i%6==0 and i!=0:
                        f.write('\n      ')
                    if i==len(source[0])-1:
                        f.write('\n')
            else:
                f.write(line)

    if found_si_string == False:
        warnings.warn('si_string was never found')
    if found_sp_string == False:
        warnings.warn('sp_string was never found')

    print("File written: {}".format(newfilename))













