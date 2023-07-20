
import pandas as pd
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt

__all__ = ['plot_h5scale_xs',"get_cross_section","get_std_comp","get_std_comp_nat_abund",
           "get_zaidlist","write_tsl_table"]

def plot_h5scale_xs(filename,scaleid,temp,emin=2.1e7,mt=None):
    """
    Plot xs of every MT (reaction number) in SCALE data H5 format

    Parameters
    ----------
    filename : str
        The path to the H5 library file
    scaleid : int
        The integer representing a single isotope or compound, e.g. 1001 
        for hydrogen, 5008016 for the oxygen in BeO 
    temp : float
        The temperature, assumed to only need one digit after the decimal
    emin : float, optional
        Minimum energy of XS must be below this value to be plotted. 
    mt : int, optional
        The MT ENDF/SCALE reaction number to be plotted. Default is all within 
        energy limit defined

    Returns
    -------
    None

    Examples
    --------
    >>> import nuctools as nuc
    >>> 
    >>> scaleid = 5008016
    >>> temperature = 293.6
    >>> nuc.plot_all_h5scale_xs('n_008016.h5',temperature,scaleid)
    """
    with h5.File(filename,'r+') as f:
        plt.figure()
        if mt is None:
            for key in f['n_{:0>7d}_{:0>6.1f}'.format(scaleid,temp)]:
                try:
                    e  = f['n_{:0>7d}_{:0>6.1f}/{}/energy'.format(scaleid,temp,key)][()]
                    xs = f['n_{:0>7d}_{:0>6.1f}/{}/xs'.format(scaleid,temp,key)][()]
                    if e.min() < emin:
                        plt.plot(e,xs,label='{}'.format(key))
                except:
                    print("key = '{}' cannot be plotted".format(key))
        else:
            try:
                e  = f['n_{:0>7d}_{:0>6.1f}/mt_{:0>4d}/energy'.format(scaleid,temp,mt)][()]
                xs = f['n_{:0>7d}_{:0>6.1f}/mt_{:0>4d}/xs'.format(scaleid,temp,mt)][()]
                if e.min() < emin:
                    plt.plot(e,xs,label='{}'.format("mt={}".format(mt)))
            except:
                print("key = 'mt_{:0>4d}' cannot be plotted".format(mt))
        plt.xlabel('Energy [eV]')
        plt.ylabel("Cross section [b]")
        plt.xscale('log')
        plt.yscale('log')
        plt.legend(ncol=3)

def get_cross_section(filename,scaleid,temp,mt):
    """
    Plot xs of every MT (reaction number) in SCALE data H5 format

    Parameters
    ----------
    filename : str
        The path to the H5 library file
    scaleid : int
        The integer representing a single isotope or compound, e.g. 1001 
        for hydrogen, 5008016 for the oxygen in BeO 
    temp : float
        The temperature, assumed to only need one digit after the decimal
    emin : float, optional
        Minimum energy of XS must be below this value to be plotted. 
    mt : int, optional
        The MT ENDF/SCALE reaction number to be plotted. Default is all within 
        energy limit defined

    Returns
    -------
    data : DataFrame
        A Pandas DataFrame with keys: "e" and "cs" for energy and cross section

    Examples
    --------
    >>> import numpy as np
    >>> 
    >>> scaleid = 5008016
    >>> temperature = 293.6
    >>> data = nuc.get_cross_section('n_008016.h5',temperature,scaleid)

    """
    with h5.File(filename,'r+') as f:
        try:
            data = pd.DataFrame({
                "e" : f['n_{:0>7d}_{:0>6.1f}/mt_{:0>4d}/energy'.format(scaleid,temp,mt)][()],
                "cs" : f['n_{:0>7d}_{:0>6.1f}/mt_{:0>4d}/xs'.format(scaleid,temp,mt)][()]
            })
        except:
            print("key = 'n_{:0>7d}_{:0>6.1f}/mt_{:0>4d}' cannot be found".format(mt))

    return data



def get_std_comp(datadir):
    """
    Return a std composition table pandas DataFrame

    assumes the directory and file names are static in the 
    SCALE data directory

    Parameters
    ----------
    datadir : str
        The directory where the "scale data" exists. This is a 
        specific Git repo

    Returns
    -------
    stdcomp : DataFrame
        A table of the standard composition from SCALE data
    """

    std_composition = datadir+'library_inputs/compoz.inp'
    begin = 0
    max_rows = 3527 # num rows in table not including mixtures
    with open(std_composition,'r+') as f:
        for i,line in enumerate(f):
            if "SCID          ROTH.    ICP      NCZA     ATPM     END" in line:
                begin_natural = i+1
    stdcomp = pd.read_fwf(std_composition,skiprows=begin_natural,nrows=max_rows,
                          colspecs=( (0,14), (15,22), (23,27), (29,39), (40,48), (49,55)),
                          names=['name','rho_theor','compound_flag','scaleid','wtpct','end'])
    stdcomp.fillna(0,inplace=True)
    stdcomp.scaleid = stdcomp.scaleid.astype('int')

    return stdcomp

def get_std_comp_nat_abund(datadir):
    """
    Return a std composition table pandas DataFrame for the 
    natural abundance table

    assumes the directory and file names are static in the 
    SCALE data directory

    Parameters
    ----------
    datadir : str
        The directory where the "scale data" exists. This is a 
        specific Git repo

    Returns
    -------
    stdcomp : DataFrame
        A table of the standard composition of natural abundance
        from SCALE data
    """

    std_composition = datadir+'library_inputs/compoz.inp'
    begin_natural = 0
    max_rows = 361 # num rows with natural abundance
    with open(std_composition,'r+') as f:
        for i,line in enumerate(f):
            if "NZN         ISZA   ABWP      END" in line:
                begin_natural = i+1
    stdcomp = pd.read_fwf(std_composition,skiprows=begin_natural,nrows=max_rows,names=['natid','isoid','frac','end'])
    stdcomp.fillna(0,inplace=True)
    stdcomp.isoid = stdcomp.isoid.astype('int')

    return stdcomp

def get_zaidlist(zlist,stdcomp):
    """
    Get the list of ZAIDs, relevant to ENDF/SCALE ID system, using the
    standard composition to get the natural composition of isotopes

    Parameters
    ----------
    zlist : list
        List of atomic numbers for unique elements in the salt mixture
    stdcomp : DataFrame
        DataFrame object made exactly like the output of `get_std_comp`

    Returns
    -------
    zaidlist : list
        The list of ZAIDs for the mixture of salts
    """
    zaidlist = []
    for z in zlist:
        for za in stdcomp.isoid:
            if( z == za//1000 ):
                zaidlist.append(za)
    return zaidlist

def write_tsl_table(library_master_file,output_file,stdcomp):
    """
    Read library master file to find TSLs and temperatures
    available and write to file

    Parameters
    ----------
    library_master_file : str
        The path for the library master file in HDF5
    output_file : str
        The path where to write a markdown table
    stdcomp : DataFrame
        The SCALE standard composition information which 
        can be placed in a DataFrame with `get_std_comp`

    Returns
    -------
    None


    """
    tsl_table = [['name','scaleid','temps']]
    with h5.File(library_master_file,'r') as f:
        last_scaleid = None
        last_scaleid_pos = 0
        for i,nuclide in enumerate(f):
            if nuclide == 'metadata' or nuclide == "nuclide_md" or "p_" in nuclide:
                continue
            found_tsl = False
            scaleid = nuclide.split('_')[1]
            temp = nuclide.split('_')[2]
            is_water = (int(scaleid) == 1001) or (int(scaleid) == 1002)
            if (int(scaleid) < 100*1000) and not is_water:
                continue
            for mtlist in f[nuclide]['reaction_md']:
                mt = mtlist[0]
                if mt == 1007 or mt == 1008:
                    found_tsl = True
                    break
            if found_tsl:
                if last_scaleid == scaleid:
                    temp = str(float(temp))
                    tsl_table[last_scaleid_pos][2] += ", {}".format(temp)
                else:
                    name = stdcomp.name[stdcomp.scaleid==int(scaleid)].to_numpy()[0]
                    temp = str(float(temp))
                    tsl_table.append([name,str(int(scaleid)),temp])
                    last_scaleid_pos += 1
                last_scaleid = scaleid
    print("Number of TSLs found = ",last_scaleid_pos)
    tsl_df = pd.DataFrame(tsl_table)
    header = tsl_df.iloc[0]
    tsl_df = tsl_df.iloc[1:len(tsl_df)]
    tsl_df.columns = header
    tsl_df.to_markdown(output_file,index=False)









