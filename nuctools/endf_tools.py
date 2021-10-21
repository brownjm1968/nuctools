
import numpy as np

__all__ = ['read_3col_pendf','write_pendf_xs']

def read_3col_pendf(file,start,finish):
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
                        e.append(word.replace('-','e-').replace('+','e+'))
                    else:
                        # -------------------------------
                        # Grab cross section, convert to readable format
                        # -------------------------------
                        if( word == '           ' ):
                            break_outer = True
                            break # end of TAB1
                        cs.append(word.replace('-','e-').replace('+','e+'))
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






