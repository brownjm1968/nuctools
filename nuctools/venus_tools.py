import tifffile as tif
import numpy as np
import glob

def tif_to_tof(folder,return_mean=False,savefile=None):
    """
    Get sum of counts for each `.tif` file from a VENUS
    experimental campaign and place in order of file name
    (alphabetically). Each file should refer
    to a different time-of-flight, ordered by slice number.

    Parameters
    ----------
    folder : str
        The path to the folder containing the .tif files
    return_mean : bool, optional
        Return the mean instead of sum of counts for a 
        slice (tof bin).
    savefile : str, optional
        The path for the file to write the data to. If the
        default of None is kept, no file is written

    Returns
    -------
    counts : array-like
        A numpy array of counts summed for each file in
        alphabetically sorted file names. Typically the
        file names have a prefix and a number indicating
        "Slice"-number for each time-of-flight bin

    Examples
    --------
    >>> import nuctools as nuc
    >>>
    >>> folder = "data/"
    >>> counts = nuc.tif_to_tof(folder)
    """
    filelist = glob.glob(f'data/{folder}/*.tif')
    filelist = np.sort(filelist)
    counts = np.zeros(len(filelist))

    for i,filename in enumerate(filelist):
        img_array = tif.imread(filename)
        if return_mean:
            counts[i] = img_array.mean()
        else:
            counts[i] = img_array.sum()
    print(f"Read {len(filelist)} files in folder: {folder}")

    if savefile is not None:
        with open(savefile,'w') as f:
            f.write('Slice,sum\n')
            for i,val in enumerate(counts):
                f.write(f"{i+1},{val:.4f}\n")
    return counts


