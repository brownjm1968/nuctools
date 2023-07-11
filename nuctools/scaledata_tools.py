
import pandas as pd
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt

__all__ = ['plot_h5scale_xs',"get_cross_section"]

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












