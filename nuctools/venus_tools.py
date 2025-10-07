import tifffile as tif
import numpy as np
import glob

def tif_to_tof(folder,return_mean=False):
    """
    Get sum of counts for each `.tif` file from a VENUS
    experimental campaign and place in order of file name
    (alphabetically). Each file should refer
    to a different time-of-flight, ordered by slice number.

    Parameters
    ----------
    folder : str
        The path to the folder containing the .tif files

    Returns
    -------
    counts : array-like
        A numpy array of counts summed for each file in
        alphabetically sorted file names. Typically the
        file names have a prefix and a number indicating
        "Slice"-number for each time-of-flight bin
    return_mean : bool, optional
        Return the mean instead of sum of counts for a 
        slice (tof bin).
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
    return counts