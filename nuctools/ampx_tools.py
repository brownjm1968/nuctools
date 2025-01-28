import pandas as pd


__all__ = ["read_tab1_ascii"]

def read_tab1_ascii(filename,mat,mf,mt,minus1=False):
    """
    Read Tab1 file in ascii format from AMPX (translated to ascii from charmin)

    Parameters
    ----------
    filename : str 
        The name of the file to be read
    mat : int
        A four digit material MAT id defined by ENDF (e.g. 7328 for Ta-181)
    mf : int
        A two digit file MF id defined by ENDF
    mt : int
        A three digit reaction MT id defined by ENDF
    minus1 : bool
        Whether the length record is wrong and needs to subtract one

    Returns
    -------
    tab1 : DataFrame
        A pandas DataFrame of energy and cross section (1-D)

    Examples
    --------
    >>> import nuctools as nuc
    >>> mat,mf,mt = 7328,3,1
    >>> filename = "mytab1.dat"
    >>> tab1 = nuc.read_tab1_ascii(filename,mat,mf,mt)

    """
    string_to_match = "{:4d}{:2d}{:3d}    5".format(mat,mf,mt)

    skip,end = 0,0
    with open(filename,'r+') as f:
        filelist = f.readlines()
    for i,line in enumerate(filelist):
        if string_to_match in line:
            skip = i+1
            end = int(line.split()[0])
            break
    if end!=0:
        if minus1:
            end -= 1
        # now on line after match
        try:
            tab1 = pd.read_csv(filename,skiprows=skip,nrows=end,names=['e','cs'],sep=r'\s+')
        except pd.errors.ParserError as err:
            print("string to match: \"{}\"".format(string_to_match))
            print("lines to skip: ",skip)
            print("lines to read: ",end)
            print(f"Unexpected {err=}, {type(err)=}")
            raise
    return tab1