import pandas as pd
import numpy as np
import scipy.integrate as intgr
import time
from . import physics_tools as pt
from . import sam_tools as st
from . import math_tools as mt

__all__ = ['tofe','etof','opsum','single_group','comp_group','gelina_group','sgfilter',
           'optstat_group','calc_notch']


def tofe(tof,FP):
    """
    Convert neutron time-of-flight into energy

    Parameters
    ----------
    tof: float
        Time of flight [us] of neutron in micro-seconds
    FP: float
        flight path length [m] of neutron in meters

    Returns
    -------
    E: 
        Energy [eV] of neutron that traveled the specified
        flight path in the specified time. (Similar to
        E=0.5mv^2)

    Notes
    -----
    Accounts for relativistic corrections.

    .. math:: E(TOF,FP) = m_n c^2 \\left( \\frac{1}{\\sqrt{1-\\frac{\\left(\\frac{FP}{TOF}\\right)^2}{c^2}}} - 1 \\right)

    Examples
    --------
    >>> import numpy as np
    >>> import nuctools as nuc
    >>> tof = np.linspace(0.1,2500,1e4)
    >>> FP = 100.141 
    >>> E = nuc.tofe(tof,FP)

    """
    #if type(tof) == float or type(tof) == int:
    #    test_tof = np.array([tof])
    #if type(tof) == list:
    #    test_tof = tof
    #    tof = np.array(tof)
    #if type(tof) == np.ndarray:
    #    test_tof = tof
    #    tof = np.array(tof)
    #if True in [FP/x/c > 1 for x in test_tof]:
    #    raise ValueError("FP/TOF/c > 1, input is incorrect.")

    E = pt.mnc2*(1/np.sqrt(1-(FP/tof/pt.c)**2)-1)
    return E
#------------------------------------------------------------------------------
def etof(E,FP):
    """
    Convert neutron energy into time-of-flight

    Parameters
    ----------
    E: float
        Energy [eV] of neutron
    FP: float
        flight path length [m] of neutron in meters

    Returns
    -------
    tof: 
        Time of flight [us] for neutron of energy E to travel
        the specified distance. (Similar to E=0.5mv^2)

    Notes
    -----
    Accounts for relativistic corrections.
    
    .. math:: TOF(E,FP) = \\frac{FP}{c}\\frac{1}{\\sqrt{1-\\frac{1}{(\\frac{E}{m_nc^2} + 1)^2}}}
        
    Examples
    --------
    >>> import numpy as np
    >>> import nuctools as nuc
    >>> E = np.linspace(10,1e6,1e4)
    >>> FP = 100.141 
    >>> E = nuc.etof(E,FP)

    """
    if type(E) == list:
        E = np.array(E)
    tof = FP/pt.c*1/np.sqrt(1-1/(E/pt.mnc2+1)**2)
    return tof
#------------------------------------------------------------------------------
def opsum(infile_list,outfile,numfiles,hlen):
    """
    Function to sum multiple norsum format open files

    Takes norsum grouped style format files and sums them, the 
    header stays the same as first file except summation of
    the triggers from all files. Keeps the same channel number,
    sums the counts, quadrature sum of the error. 

    .. math:: \\sqrt{a^2+b^2}

    Parameters
    ----------
    infile_list: list
        Strings of filenames, e.g. "C:/Users/example.txt"
    outfile: str
        Filename for newly summed data
    numfiles: int
        Number of files to be summed into outfile
    hlen: int
        Header length: how many lines are in the header
    
    Returns
    -------
    None
        Writes newly summed data to text file.

    """
    file_array = []
    # create array of files where each file is vector of line strings
    for infile in infile_list:
        with open(infile,"r+") as file:
            file_array.append(file.readlines())
    print("Number of files: ",len(file_array))
    # write the header and break after triggers line
    with open(outfile,"w+") as newfile:
        breakFlag = False
        tot_trig = 0
        for i,line in enumerate(file_array[0]):
            if breakFlag:
                break
            # sum the triggers from all files
            if "total number of triggers" in line:
                breakFlag = True
                for file in file_array:
                    print(int(file[i].split()[0]))
                    # get the triggers from each file and sum them
                    tot_trig += int(file[i].split()[0])
                line = "{}                      total number of triggers\n".format(tot_trig)
            newfile.write(line)
    # free memory
    del file_array
    data = []
    for file in infile_list:
        print(file)
        data.append(np.genfromtxt(file,unpack=True,skip_header=hlen))
    # The shape of the data should be ( #files x #cols x #rows ) = 3 dims.
    if len(np.shape(data)) < 3:
        raise ValueError("Files do not match. ( # rows, # cols )")
    # Throw an error if the number of columns in the unpacked file is not 3
    if np.shape(data)[1] != 3:
        raise ValueError("Shape of data is unexpected. File format may be incorrect.")
    data = np.array(data)
    # sum the counts across the files
    counts = data.sum(axis=0)[1]
    # quadrature sum of the errors
    err = np.sqrt((data**2).sum(axis=0)[2])
    sum_data = np.vstack((data[0][0],counts,err)).transpose()
    # append summed data to the outfile
    with open(outfile,"ab") as newfile:
        np.savetxt(newfile,sum_data,fmt=["%-15d","%8.5f","%16.10f"])
    # free memory
    del sum_data
    print("done")
#------------------------------------------------------------------------------
def single_group(data,factor,avg=False,binedge=False,prop_err=False):
    """
    Group time of flight data by a single factor

    Grouping is done using the numpy function resize. Each column of data
    is broken into len(data[0])/factor rows, each row is length factor, 
    to find group tof each row of data[0] is averaged (default), each 
    row of data[1] is summed (default).

    Parameters
    ----------
    data : numpy array 
        First column of data should be tof, second 
        column of data should be counts. If propagation
        of error is desired, 3rd column of data should
        be the error, otherwise 3rd column is optional
    factor : int
        how many bins to group by
    avg : bool, optional
        Return an average of group values instead of
        the sum in the group
    binedge : bool, optional
        Return TOF column for grouped data as the first
        tof in group instead of the mean tof of the
        group.
    prop_err : bool, optional
        IF using the avg option, propagate the error given
        in the data array instead of calculating the 
        standard deviation of the mean of each group
    Returns
    -------
     : 2-d numpy array
        Grouped data with columns: tof,bin_width,counts,error counts

    Notes
    -----
    Data array must have more rows than columns, this is 
    typical of TOF data. 

    Examples
    --------
    >>> import nuctools as nuc
    >>> tof = np.linspace(100,1200,400000)
    >>> counts = nuc.noisy_gauss(tof,170,7)
    >>> data = [tof,counts]
    >>> gdata = nuc.single_group(data,2)
    >>> plt.figure(3)
    >>> plt.plot(tof,gauss,color="k",ls=" ",marker="o")
    >>> plt.plot(gdata[0],gdata[1])
    >>> plt.show()

    """
    # data is array with tof & counts for columns
    ld = max(np.shape(data[0]))
    # col. len. should be greatest integer division possible
    col_len = ld//factor
    # if user asked for leading tof bin edge
    if binedge == True:
        # tof for bin is time of first element in bin vector
        gtof = np.resize(data[0],(col_len,factor))[:,0]
    else:
        # tof for bin is mean of tof bin vector
        gtof = np.resize(data[0],(col_len,factor)).mean(axis=1)
    if avg == True:
        # counts for bin is the mean of all elements in bin vector
        gcounts = np.resize(data[1],(col_len,factor)).mean(axis=1)
        if prop_err == True:
            # sqrt of the sum of the squares
            ergcounts = 1/factor*np.sqrt(np.resize(data[2]**2,(col_len,factor)).sum(axis=1))
        else:
            stdcounts = np.resize(data[1],(col_len,factor)).std(axis=1)
            ergcounts = stdcounts/np.sqrt(factor)
    else:
        # counts for bin is the sum of all elements in bin vector
        gcounts = np.resize(data[1],(col_len,factor)).sum(axis=1)
        ergcounts = np.sqrt(gcounts)
    # difference between time bins to plot cps
    dt = np.diff(gtof)
    # diff is calc. between points, therefore has n-1 pts -- add 1 at 0
    dt = np.insert(dt,0,dt[0])
    
    return np.array([gtof,dt,gcounts,ergcounts])
#------------------------------------------------------------------------------
def comp_group(data,factor,c_pts,FP,avg=False,binedge=False,verbose=False):
    """
    Group counts data with compression points as a function of TOF.

    Group TOF data with different factors i.e compression points according
    to the group selected by the compression points. Function finds indices
    of the tof data (data[0]). These tof indices match closely the TOF's given
    in c_pts. Then the function forms len(c_pts)+1 compression regions, and groups each 
    compression region by the values specified by the user in parameter: factor.

    Parameters
    ----------
    data : 2-d numpy array
        Array with columns "tof","counts", and "error".
    factor : 1-d numpy array
        Compression factors for the tof groups, 1st compression factor
        corresponds to time of flights below c_pt[0], 2nd corresponds to tofs above
        c_pt[0] and below c_pt[1], and so on until end of range is found.
    c_pts : 1-d numpy array
        Compression time of flights that seperate the bins that are to be grouped.
        List points from lowest to highest (should re-write to sort and remove user
        capability to fail.)
    FP : float
        Flight path needed to convert the time of flight values of bins to
        energies for the compression point selection.
    avg : bool, optional
        Return an average of group values instead of
        the sum in the group. group error is 1/N*(sum(error^2))^(1/2)
    binedge : bool, optional
        Return a TOF column for grouped data as the first
        tof in group as well as the mean tof of the group.
    verbose : bool, optional
        Choose whether to print out approximate comp pt energy locations
        and indices. Default is False.

    Returns
    -------
    gdata : 2-d numpy array
        Returns the grouped data with columns: tof,bin_width,counts,error counts
    ascii-file : text file (NOT YET)
        Contains a header of bin widths and indices of where groups end, 
        and the grouped data.

    Notes
    -----
    The number of values in tof data, must be greater than 3 (number of columns 
    defined in the description of parameter: data[tof_i,counts_i,errCounts_i])

    The tof vector no longer contains all the information for bin width, it now
    describes the average tof for each grouped bin. This can cause problems at the
    compression points when converting from counts to cps. Use the output bin_width
    to calculate cps.

    Examples
    --------
    >>> import nuctools as nuc
    >>> import matplotlib as plt
    >>> tof = np.linspace(100,1200,1000)
    >>> counts = nuc.noisy_gauss(tof,170,7) + nuc.noisy_gauss(tof,270,15) 
    >>> counts += nuc.noisy_gauss(tof,400,20)
    >>> counts *= 1000
    >>> ercounts = np.sqrt(counts)
    >>> FP = 45.25
    >>> factor,c_pts = np.array([2,4,8]),np.array([10,100])
    >>> data = [tof,counts,ercounts]
    >>> gdata = nuc.comp_group(data,factor,c_pts,FP)
    >>> plt.figure(1)
    >>> plt.plot(tof,counts,color="blue",label="orig")
    >>> plt.plot(gdata[0],gdata[1],label="grouped")
    >>> plt.legend()
    >>> plt.show()

    """
    # cast data as a numpy array
    data = np.array(data)
    # ensure shape is correct (pandas input would've caused probs.)
    data = np.reshape(data,np.sort(np.shape(data)))
    # if statement for compression pts out of tof range
    for val in c_pts:
        if val > max(data[0]) or val < min(data[0]):
            print("Tof: ",max(data[0]),min(data[0]))
            raise ValueError("Compression pt. is out of time-of-flight range.")
    # should be more than one factor
    if len(factor) < 2:
        raise ValueError("Only one factor was given, use single_group func.")
    # cast factor and c_pts to np.array()
    factor,c_pts = np.array(factor),np.array(c_pts)
    if verbose == True:
        print("Comp. pt at ",tofe(c_pts,FP)," eV")
    if len(c_pts)+1 != len(factor):
        raise ValueError("Length of c_pts should be 1 less than length factor")
    # --- Selecting compression regions -----------------------------
    c_ind,j = [0],0
    ng,ld = len(factor),len(data[0])
    # for each compression point go find the index where it lands in the
    # tof grid, and select an indice close to that tof where the grouping
    # can include all bins (except at the end of tof range)
    for i,pt in enumerate(c_pts):
        # index = (tof-first time of flight)/bin_width
        index = int((pt-data[0][0])/(data[0][1]-data[0][0]))
        while (index-c_ind[i]) % factor[j]:
            index += 1
        j += 1
        c_ind.append(index)
        if verbose == True:
            print("compression indices ",index)
    if c_ind[len(c_ind)-1] + 2*factor[len(factor)-1] > ld:
        print('You\'ve grouped too hard on the last compression')
        print('region. You need at least 2 bins for each, and ')
        print('there is only room for one.')
        raise ValueError("Redefine factor for last comp. region.")
    c_ind = c_ind[1:len(c_ind)]
    # ---------------------------------------------------------------
    gdata = []
    # loop over multiplication factors: mult
    for i, mult in enumerate(factor):
        # if in first compression point group
        if i == 0:
            # shape of temporary array (rows,cols)
            # rows are how many groups we'll have, cols are how 
            # many elements are being binned
            shape = (c_ind[i]//mult,mult)
            # tof for bin is mean time in bin vector
            # resize the data and take the average of the row
            #gtof = [vector.mean() for vector in np.resize(data[0][0:c_ind[0]],shape)]
            gtof = np.resize(data[0][0:c_ind[0]],shape).mean(axis=1)
            # if user selected to return tof of leading bin edge
            if binedge == True:
                # tof for bin is time of first element in bin vector
                edge_gtof = np.resize(data[0][0:c_ind[0]],shape)[:,0]
            # calculate the average counts in data[1]
            if avg == True:
                gcounts = np.resize(data[1][0:c_ind[0]],shape).mean(axis=1)
                ergcounts = 1/mult*np.sqrt(np.resize(data[2][0:c_ind[0]]**2,shape).sum(axis=1))
            # calculate the sum of counts in data[1]
            else:
                # counts for bin is the sum of all elements in bin vector
                # resize the data and take the sum of the row
                gcounts = np.resize(data[1][0:c_ind[0]],shape).sum(axis=1)
                ergcounts = np.sqrt(gcounts)
            # difference between time bins to plot cps
            dt = np.diff(gtof)
            # diff is calc. between points, therefore has n-1 pts -- add 1 at 0
            dt = np.insert(dt,0,dt[0])
            # if user asked for leading tof bin edge append it as an additional 
            # column to gdata
            if binedge == True:
                gdata.append([gtof,dt,gcounts,ergcounts,edge_gtof])
            else:
                gdata.append([gtof,dt,gcounts,ergcounts])
            continue
        # if in the last compression group (in tof)
        if i == ng-1:
            # shape of temporary array
            shape = ((ld-c_ind[i-1])//mult,mult)
            # tof for bin is mean time in bin vector
            gtof = np.resize(data[0][c_ind[i-1]:ld],shape).mean(axis=1)
            # if user selected to return tof of leading bin edge
            if binedge == True:
                # tof for bin is time of first element in bin vector
                edge_gtof = np.resize(data[0][c_ind[i-1]:ld],shape)[:,0]
            # calculate the average counts in data[1]
            if avg == True:
                gcounts = np.resize(data[1][c_ind[i-1]:ld],shape).mean(axis=1)
                ergcounts = 1/mult*np.sqrt(np.resize(data[2][c_ind[i-1]:ld]**2,shape).sum(axis=1))
            # calculate the sum of counts in data[1]
            else:
                # counts for bin is the sum of all elements in bin vector
                gcounts = np.resize(data[1][c_ind[i-1]:ld],shape).sum(axis=1)
                ergcounts = np.sqrt(gcounts)
            # difference between time bins to plot cps
            dt = np.diff(gtof)
            # diff is calc. between points, therefore has n-1 pts -- add 1 at 0
            dt = np.insert(dt,0,dt[0])
            # if user asked for leading tof bin edge append it as an additional 
            # column to gdata
            if binedge == True:
                gdata.append([gtof,dt,gcounts,ergcounts,edge_gtof])
            else:
                gdata.append([gtof,dt,gcounts,ergcounts])
            continue
        else:
            # shape of temporary array
            shape = (((c_ind[i]-(c_ind[i-1]))//mult),mult)
            # tof for bin is avg of elements in bin vector
            # resize the data and take the average of the row
            gtof = np.resize(data[0][c_ind[i-1]:c_ind[i]],shape).mean(axis=1)
            # if user selected to return tof of leading bin edge
            if binedge == True:
                # tof for bin is time of first element in bin vector
                edge_gtof = np.resize(data[0][c_ind[i-1]:c_ind[i]],shape)[:,0]
            # calculate the average counts in data[1]
            if avg == True:
                gcounts = np.resize(data[1][c_ind[i-1]:c_ind[i]],shape).mean(axis=1)
                ergcounts = 1/mult*np.sqrt(np.resize(data[2][c_ind[i-1]:c_ind[i]]**2,shape).sum(axis=1))
            # calculate the sum of counts in data[1]
            else:
                # counts for bin is the sum of all elements in bin vector
                # resize the data and take the sum of the row
                gcounts = np.resize(data[1][c_ind[i-1]:c_ind[i]],shape).sum(axis=1)
                ergcounts = np.sqrt(gcounts)
            # difference between time bins to plot cps
            dt = np.diff(gtof)
            # diff is calc. between points, therefore has n-1 pts -- add 1 at 0
            dt = np.insert(dt,0,dt[0])
            # if user asked for leading tof bin edge append it as an additional 
            # column to gdata
            if binedge == True:
                gdata.append([gtof,dt,gcounts,ergcounts,edge_gtof])
            else:
                gdata.append([gtof,dt,gcounts,ergcounts])
            
    # hstack in this case is appending all the gdata arrays into one array 
    # with columns gtof,gcounts,ergcounts
    gdata = np.hstack([val for val in gdata])
        
    return gdata
#------------------------------------------------------------------------------
def gelina_group(data,factor,c_ind,avg=False,binedge=False,verbose=False):
    """
    Group like AGL for GELINA comparison. This is experimental for now.

    Parameters
    ----------
    data : 2-d numpy array
        Array with columns "tof","counts", and "error".
    factor : 1-d numpy array
        Compression factors for the tof groups, 1st compression factor
        corresponds to time of flights below c_pt[0], 2nd corresponds to tofs above
        c_pt[0] and below c_pt[1], and so on until end of range is found.
    c_pts : 1-d numpy array
        Compression points that define the edges of the bins that are to be grouped.
        List points from lowest to highest (should re-write to sort and remove user
        capability to fail.)
    verbose : bool
        Whether to print out information to STDOUT

    Returns
    -------
    gdata : 2-d numpy array
        Returns the grouped data with columns: tof,bin_width,counts,error counts

    """
    ng,ld = len(factor),len(data[0])
    gdata = []
    # loop over multiplication factors: mult
    for i, mult in enumerate(factor):
        # if in first compression point group
        if i == 0:
            # shape of temporary array (rows,cols)
            # rows are how many groups we'll have, cols are how 
            # many elements are being binned
            shape = (c_ind[i]//mult,mult)
            # tof for bin is mean time in bin vector
            # resize the data and take the average of the row
            #gtof = [vector.mean() for vector in np.resize(data[0][0:c_ind[0]],shape)]
            gtof = np.resize(data[0][0:c_ind[0]],shape).mean(axis=1)
            # if user selected to return tof of leading bin edge
            if binedge == True:
                # tof for bin is time of first element in bin vector
                edge_gtof = np.resize(data[0][0:c_ind[0]],shape)[:,0]
            # calculate the average counts in data[1]
            if avg == True:
                gcounts = np.resize(data[1][0:c_ind[0]],shape).mean(axis=1)
                ergcounts = 1/mult*np.sqrt(np.resize(data[2][0:c_ind[0]]**2,shape).sum(axis=1))
            # calculate the sum of counts in data[1]
            else:
                # counts for bin is the sum of all elements in bin vector
                # resize the data and take the sum of the row
                gcounts = np.resize(data[1][0:c_ind[0]],shape).sum(axis=1)
                ergcounts = np.sqrt(gcounts)
            # difference between time bins to plot cps
            dt = np.diff(gtof)
            # diff is calc. between points, therefore has n-1 pts -- add 1 at 0
            dt = np.insert(dt,0,dt[0])
            # if user asked for leading tof bin edge append it as an additional 
            # column to gdata
            if binedge == True:
                gdata.append([gtof,dt,gcounts,ergcounts,edge_gtof])
            else:
                gdata.append([gtof,dt,gcounts,ergcounts])
            continue
        # if in the last compression group (in tof)
        if i == ng-1:
            # shape of temporary array
            shape = ((ld-c_ind[i-1])//mult,mult)
            # tof for bin is mean time in bin vector
            gtof = np.resize(data[0][c_ind[i-1]:ld],shape).mean(axis=1)
            # if user selected to return tof of leading bin edge
            if binedge == True:
                # tof for bin is time of first element in bin vector
                edge_gtof = np.resize(data[0][c_ind[i-1]:ld],shape)[:,0]
            # calculate the average counts in data[1]
            if avg == True:
                gcounts = np.resize(data[1][c_ind[i-1]:ld],shape).mean(axis=1)
                ergcounts = 1/mult*np.sqrt(np.resize(data[2][c_ind[i-1]:ld]**2,shape).sum(axis=1))
            # calculate the sum of counts in data[1]
            else:
                # counts for bin is the sum of all elements in bin vector
                gcounts = np.resize(data[1][c_ind[i-1]:ld],shape).sum(axis=1)
                ergcounts = np.sqrt(gcounts)
            # difference between time bins to plot cps
            dt = np.diff(gtof)
            # diff is calc. between points, therefore has n-1 pts -- add 1 at 0
            dt = np.insert(dt,0,dt[0])
            # if user asked for leading tof bin edge append it as an additional 
            # column to gdata
            if binedge == True:
                gdata.append([gtof,dt,gcounts,ergcounts,edge_gtof])
            else:
                gdata.append([gtof,dt,gcounts,ergcounts])
            continue
        else:
            # shape of temporary array
            shape = (((c_ind[i]-(c_ind[i-1]))//mult),mult)
            # tof for bin is avg of elements in bin vector
            # resize the data and take the average of the row
            gtof = np.resize(data[0][c_ind[i-1]:c_ind[i]],shape).mean(axis=1)
            # if user selected to return tof of leading bin edge
            if binedge == True:
                # tof for bin is time of first element in bin vector
                edge_gtof = np.resize(data[0][c_ind[i-1]:c_ind[i]],shape)[:,0]
            # calculate the average counts in data[1]
            if avg == True:
                gcounts = np.resize(data[1][c_ind[i-1]:c_ind[i]],shape).mean(axis=1)
                ergcounts = 1/mult*np.sqrt(np.resize(data[2][c_ind[i-1]:c_ind[i]]**2,shape).sum(axis=1))
            # calculate the sum of counts in data[1]
            else:
                # counts for bin is the sum of all elements in bin vector
                # resize the data and take the sum of the row
                gcounts = np.resize(data[1][c_ind[i-1]:c_ind[i]],shape).sum(axis=1)
                ergcounts = np.sqrt(gcounts)
            # difference between time bins to plot cps
            dt = np.diff(gtof)
            # diff is calc. between points, therefore has n-1 pts -- add 1 at 0
            dt = np.insert(dt,0,dt[0])
            # if user asked for leading tof bin edge append it as an additional 
            # column to gdata
            if binedge == True:
                gdata.append([gtof,dt,gcounts,ergcounts,edge_gtof])
            else:
                gdata.append([gtof,dt,gcounts,ergcounts])
            
    # hstack in this case is appending all the gdata arrays into one array 
    # with columns gtof,gcounts,ergcounts
    gdata = np.hstack([val for val in gdata])
        
    return gdata


#------------------------------------------------------------------------------
def optstat_group(data,FP,D,res_per_bin,E_cpts=[1e2,1e3,1e4],verbose=False,
                  t0=0.0,binedge=False):
    """
    Function that finds the grouping necessary to include at least res_per_bin
    resonances in each group.

    Finds the grouping scheme necessary to include at least res_per_bin 
    resonances in each group and then groups the provided data using the 
    nuctools comp_group function. The E_cpts are estimated compression pts
    for the grouper. This is typically helpful in the URR when you want to 
    group with a \"statiscally significant\" number of resonances.

    Parameters
    ----------
    data : array-like
        Data array. Should be a tof vector at data[0], data[0] must
        have equal width time bins. data[1] should be counts, data[1]
        will be grouped and summed.
    FP : float
        Flight path of the experiment, this helps in determining the
        energy to time ratio
    D : float
        nuclear energy level / resonance spacing
    res_per_bin : int
        The number of resonances per bin that are desired
    E_cpts : array-like, optional
        The compression pts desired, given in energy
    verbose : bool, optional
        Whether to spit out extra information about the grouping
    t0 : float
        Use t0 to determine SAMMY input printout
    binedge : bool, optional
        Return a TOF column for grouped data as the first
        tof in group as well as the mean tof of the group.

    Returns
    -------
    gdata : numpy array
        Returns the grouped data with columns: tof,bin_width,counts,error counts
    """
    # min size in energy
    min_bin_size = D*res_per_bin
    # tof compression points
    tof_pts = etof(E_cpts,FP)
    # bin width in time
    bw = data[0][1]-data[0][0]
    
    data_region_start_pts = np.hstack((tof_pts[::-1],data[0][len(data[0])-1]))
    data_region_start_pts_E = tofe(data_region_start_pts,FP)
    # min size in time
    min_tof_bin = data_region_start_pts-etof(data_region_start_pts_E+min_bin_size,FP)
    # what factor is required for min size
    g_factors = np.array(min_tof_bin/bw).astype(int)
    # don't group lowest TOF values so hard
    g_factors[len(g_factors)-1] = g_factors[len(g_factors)-2]
    # group highest TOF values 3 times as hard
    g_factors[0] = g_factors[0]*3
    if verbose:
        print("E cpts:",data_region_start_pts_E)
        print("Min. bin size [us]:",min_tof_bin)
        print("Grouping factors:",g_factors)
        print("TOF comp. pts:",tof_pts[::-1])
        print("--------------")
        st.print_chann(tof_pts[::-1],g_factors,bw*1e3,t0,FP)
    
    gdata = comp_group(data,g_factors,tof_pts[::-1],FP,binedge=binedge)
    
    return gdata
#------------------------------------------------------------------------------
def sgfilter(f,df,window,poly,deriv=0):
    """
    Performs a Savitsky-Golay filter on data, returning the
    smoothed data and the newly weighted errors on the 
    smoothed data.
    
    Parameters
    ----------
    f : array_like
        The function or data to be smoothed
    df : array_like
        The error associated with f
    window : int
        The window over which to weight the data
    poly : int
        The order of polynomial to used to fit the data
    deriv : int, optional
        The order of derivative to compute. Must be non-
        negative. Default is 0.
    
    Returns
    -------
    s : ndarray
        2d array that has smoothed data in column 0, and
        the error weighted by the smoothing in column 1.

    Notes
    -----
    
    The design matrix :math:`\\mathbf{A}` for a smoothing window of 5 points and a polynomial fit of 
    order 2 is:

    .. math:: \\mathbf{A} = \\left[\\begin{array}{ccc} -2^0 & -2^1 & -2^2 \\cr -1^0 & -1^1 & -1^2 \\cr  0^0 &  0^1 &  0^2 \\cr  1^0 &  1^1 &  1^2 \\cr  2^0 &  2^1 &  2^2 \\cr \\end{array}\\right]

    Or for a smoothing window of length :math:`w`, and polynomial of order :math:`M` where 
    :math:`nL = -w\\div2`, (note that :math:`\\div` denotes integer division) and :math:`nR = nL+w`:

    .. math:: \\mathbf{A}= \\left[\\begin{array}{ccc} -nL^0 & ...   & -nL^M \\cr \\vdots & \\vdots & \\vdots \\cr 0^0   & ...   & 0^M   \\cr \\vdots & \\vdots & \\vdots \\cr nR^0  & ...   & nR^M  \\cr \\end{array}\\right]

    The Savitsky-Golay Coefficients :math:`\\mathbf{C}_{sv}` are given by:

    .. math:: \\mathbf{C}_{sv} = \\left(\\mathbf{A}^T\\mathbf{A}\\right)^{-1}\\mathbf{A}^T

    The smoothed signal :math:`Y_{i,j}` comes from the convolution (:math:`\\circledast`) of the original signal :math:`y_i` (:math:`i\\in0...N`) with each column :math:`j` in the matrix :math:`\\mathbf{C}_{sv}`:

    .. math:: Y(t)_j = \\mathbf{C}_{sv,j}\\circledast y_i = \\sum_{k=-nL}^{nR} \\mathbf{C}_{sv,k+nL} \\cdot y_{i+k}

    for (i :math:`\\in nL...N-nR-1`)

    Examples
    --------
    >>> import nuctools as nuc
    >>> import matplotlib.pyplot as plt
    >>> tof = np.linspace(100,250,150)
    >>> counts = nuc.noisy_gauss(tof,170,7)
    >>> counts *= 1000
    >>> ercounts = np.sqrt(counts)
    >>> sm_counts = nuc.sgfilter(counts,ercounts,5,0)
    >>> plt.figure()
    >>> plt.errorbar(x,counts,ercounts,label="orig. data")
    >>> plt.errorbar(x,sm_counts[0],sm_counts[1],label="sm. data")
    >>> plt.legend()
    >>> plt.show()

    .. image:: images/sav_gol_example.png

    """
    
    
    if poly >= window:
        raise ValueError("polyorder must be less than window_length.")

    pos, rem = divmod(window, 2)

    if rem == 0:
        raise ValueError("window_length must be odd.")
    
    nL,nR = pos,window-pos
    
    # The design matrix A
    A = np.array([[i**j for i in range(-nL,nR)] for j in range(poly+1)])
    B = np.linalg.inv(np.dot(A,np.transpose(A)))
    # The Savitsky-Golay coefficients
    SV = np.transpose(np.dot(np.transpose(A),B))

    lf = len(f)
    smooth,sigsmooth = np.array(f),np.array(df)**2
    # ------ leaving this here in case speedy version is confusing -----------
    # for each point to smooth
    #for i in range(nL,lf-nR-1):
        #smooth[i] = 0
        # for each point on the window
        #smooth[i] = (SV[deriv][0:nL+nR]*f[i-nL:i+nR]).sum()
        #sigsmooth[i] = ((SV[deriv][0:nL+nR]*df[i-nL:i+nR])**2).sum()
    smooth[nL:lf-nR+1] = np.convolve(SV[deriv],f,mode='valid')
    sigsmooth[nL:lf-nR+1] = np.convolve(SV[deriv]**2,df**2,mode='valid')
    # ------------------------------------------------------------------------

    dsmooth = np.sqrt(sigsmooth)

    s = np.vstack([smooth,dsmooth])
    
    return s

#------------------------------------------------------------------------------
def calc_notch(spect,lim1,lim2):
    """
    Calculate the average count rate at time-of-flight (TOF) values
    that are within the limits lim1, lim2

    Parameters
    ----------
    spect : Pandas DataFrame
        Must have columns of: "tof","cps","dcps"
    lim1 : float
        lower TOF value limiting average
    lim2 : float
        higher TOF value limiting average

    Returns
    -------
    notch : array
        Numpy array with 3 values: avg TOF, avg count rate, error on avg count rate
    """
    condition = (spect.tof>lim1) & (spect.tof<lim2)
    avg = np.mean(spect.cps[condition])
    avgerr = mt.prop_err_mean(spect.dcps[condition])
    avgtof = np.mean(spect.tof[condition])

    notch = np.array([avgtof,avg,avgerr])
    return notch









