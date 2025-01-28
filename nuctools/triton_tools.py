import pandas as pd

__all__ = ["read_depletion_history"]

def read_depletion_history(triton_inp):
    """
    Read the depletion history from a SCALE/TRITON input file
    and put it into a pandas DataFrame

    Parameters
    ----------
    triton_inp : str
        The path and filename for the triton input file

    Returns
    -------
    history : DataFrame
        A DataFrame with columns: 'power','burn','str','down','nlib','end'
        where `power` is the power in MW and `burn` is the days burning
        the fuel, and `down` is the downtime between burns
    """
    skip_lines = 0
    num_lines = 0
    with open(triton_inp,'r') as f:
        lines = f.readlines()
    for i,line in enumerate(lines):
        if "read burndata" in line:
            skip_lines = i+1
        if "end burndata" in line and skip_lines != 0:
            num_lines = i-skip_lines
    if skip_lines == 0 and num_lines == 0:
        raise ValueError("Could not find burndata in input file")

    history = pd.read_csv(triton_inp,skiprows=skip_lines,nrows=num_lines,sep=r'\s+',
                          names=['power','burn','str','down','nlib','end'])
    history['power'] = history.power.str.replace('power=','').astype(float)
    history['burn'] = history.burn.str.replace('burn=','').astype(float)
    print("Total burnup = {:.4f} GWd/t".format((history.power*history.burn).sum()/1000))

    return history

