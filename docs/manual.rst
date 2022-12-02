============
Introduction
============

Nuctools has been developed to consistently reduce raw data collected from 
the RPI Electron Linear Accelerator to functionals of cross section, i.e.
transmission and yield. Starting with data collected from the Struck
SIS-3305 digitizer operated on the 40m capture system, this manual will
discuss methods of data collection, file format, and more.

==============
HDF5 Structure
==============

In the code ``40mProcess.exe`` HDF5 files are built from the raw data that
have been stored in binary. The HDF5 values being stored in normal operation
(this could change) include:

* TOF : time of flight
* PA : pulse area
* Ch : channel (on the digitizer board)
* TStamp : time stamp
* TrigNo : linac trigger number

When using ``h5py`` to open and manipulate files use ``import h5py as h5``
to import the library and open the files with the ``hdf = h5.File(filename)`` 
command. A very helpful commmand to show the internal file structure of
``hdf`` is the ``visit()`` function, by creating a function such as::


    def printname(name):
        print(name)

and operating this function with::

    hdf.visit(printname)

========================
AGL Input File Structure
========================
This new version of the AGL input file is in JSON format. There is an example template given in the ``tests`` folder. The required keys are:
- 'numadc', 'llduld', 'mcabins', 'unweighted', 'verbose', 'useAGLgrouping', 'ecal', 'wfcoeff', 'bin_width', 'max_tof', 'cfct', 'zones', 'badrundict'

These are defined by:

badrundict : dict
    A python dictionary with keys that match file ID and run number (e.g. 
    ``zr90_fp14a_f01_r01`` or ``zr90_fp14a_f01_r02``, and values that are lists 
    of integers for file numbers
user_max_tof : float
    The maximum time-of-flight value to histogram [us]
bin_width : float
    The base bin width the TOF data were measured with [us].
ecal : array-like
    A 2-dimensional array of N rows and 3 columns. Each column is for 
    parameter ``A0,A1,A2``. The number of rows is determined by the number of
    ADCs in the data file. Accessing `A1` for ADC1, e.g., would be: 
    ``ecal[0][1]``
llduld : array-like
    A 2-dimensional array of N rows and 2 columns. Column 0 represents the 
    lower-level discriminator (LLD) and column 1 the upper-level discriminator
    (ULD). The number of rows is determined by the number of ADCs in the file
wfcoeff : array-like
    A 1-d array for weighting function coefficients defined by the GELINA 
    weighting function (see ``nuctools.agl_tools.read_and_add()``).
mcabins : int
    Number of bins for an MCA spectrum, **if zero no mca is given**
numadc : int, optional
    The number of ADCs listed **in the binary file. This will change the way**
    **the binary file is read.** Often there will be either 2 or 4 ADCs, 4 is
    the default.
unweighted : bool, optional
    Whether to add unweighted histograms to the DataFrame 
verbose : bool, optional
    Whether to increase print output
bin_width : float 
    describing base width TOF bin
cfct : list 
    list of bin-grouping factors
zones : list 
    a list of zones as defined by AGL software

AGL defines "zones" as 1024 bins, so zones of [2,2] with "cfct" of
[1,2] means that 2048 bins have width 2^{1}*bin_width and the 2048 bins
following the first set have width 2^{2}*bin_width. This is GELINA-style
grouping language.


=========================
Installing CUDA on Ubuntu
=========================

Instructions initially found on this Ubuntu `forum`_, and Nvidia provides instructions for 
installing CUDA-9.2 `at this website`_.

.. _forum: https://askubuntu.com/questions/799184/how-can-i-install-cuda-on-ubuntu-16-04
.. _at this website: https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html

Verify you have a CUDA capable GPU, and accepted version of gcc::

    lspci | grep -i nvidia
    gcc --version

Make sure you have the proper kernel headers::

    uname -r
    sudo apt-get install linux-headers-$(uname -r)

Download the `92 toolkit`_, which should be a .run file. Run the md5sum on the file, and check 
if it matches what Nvidia says it should be for the `md5sum`_ , and don't continue unless it's 
correct::

    md5sum <file name>

.. _92 toolkit: http://developer.nvidia.com/cuda-downloads
.. _md5sum: http://developer.nvidia.com/cuda-downloads/checksums

Remove any other installation or drivers::

    sudo apt-get purge nvidia-*

Log out from your GUI. Once at the login screen, open a terminal session (ctrl+alt+f2) and 
stop the lightdm:: 

    sudo service lightdm stop

Navigate into ``/etc/modprobe.d/`` and create the file ``blacklist-nouveau.conf`` , for me ::

    nano blacklist-nouveau.conf

Followed by entering the lines below, (ctrl+o) to save, (ctrl+x) to exit::

    blacklist nouveau
    options nouveau modeset=0

Then regenerate the kernel initramfs::

    sudo update-initramfs -u

Then navigate to the where you've downloaded the cuda toolkit .run file, for me ``Downloads``,
and run the installation file. I used the key word ``--override`` which: ``Ignores compiler, 
third-party library, and toolkit detection checks which would prevent the CUDA Toolkit and CUDA 
Samples from installing.`` Up to you if you want to use it.::

    cd /home/brownj25/Downloads/
    sudo sh cuda_9.2.88_396.26_linux.run --override

Follow the prompts, make sure you say yes to creating a symbolic link. If the installation did not
work for you, I suggest scouring the Nvidia distributed installation guide to make sure everything
is hunky dory. After the toolkit has been installed, restart the lightdm::

    sudo service lightdm start

Then edit your ``.bashrc`` file to include nvidia executables and libraries, add the lines::

    PATH=$PATH:/usr/local/cuda-9.2/bin
    LD_LIBRARY_PATH=/usr/local/cuda-9.2/lib64:$LD_LIBRARY_PATH

and source the ``.bashrc``::

    source ~/.bashrc

Now check if you've got it installed. Type::

    nvcc --version
    nvidia-smi

This should return the version information, and then some device information. Following this,
make sure your installation works, navigate to the ``NVIDIA_CUDA-9.2_Samples/`` folder, and 
compile the sample code::

    make

Then navigate into the folder in ``bin/`` containing ``deviceQuery`` and type::

    ./deviceQuery

This should inform you that the installation and sample files have passed.


=====================
Running SESH & FITACS
=====================


The SESH code uses sample and average resonance parameter dependent Monte Carlo calculations to 
correct transmission and capture cross section for resonance self-shielding and multiple 
scattering. The FITACS code can then fit this corrected cross section for new average resonance 
parameters. Using these two codes to correct and fit experimental cross section data requires 
that they use the same average resonance parameters. To achieve this, the output average 
resonance parameters from FITACS are fed to SESH, which then calculates a new correction for the 
experimental data, which feeds back to FITACS. This process is iterated until the correction is 
no longer changing by more than 1%. 

In the ``nuctools.urr_tools`` module this is process can be operated by the ``sesh_fitacs()`` 
function. This function requires that you have 4 things:

- List of DataFrames containing all data fitted by FITACS (in the same order)
- A yaml format input file
- Operational FITACS and SESH input files
- ``sesh`` and ``sammy`` executables

------------------
List of DataFrames
------------------

The list of DataFrames must be in the same order as the order of data files fit by FITACS. 
Each of the DataFrames must have properly named columns. For total cross section you need:

- e : energy in eV
- cs : cross section in barns
- dcs : absolute uncertainty on the cross section in barns
- t : the transmission corresponding to the cross section for this sample
- dt : the absolute uncertainty on the transmission

For capture cross section you need:

- e : energy in eV
- cs : cross section in barns
- dcs : absolute uncertainty on the cross section in barns

Best practice for now is reading the files that will be fit by FITACS into a DataFrame, and 
organizing it appropriately for the ``sesh_fitacs()`` runner. An example of reading the 
FITACS data files into DataFrames and adding them to a list is given below ::

    >>> folder = '/Users/jesse/data/'
    >>> 
    >>> totxs_ta1 = pd.read_csv(folder+"ta1_sig.dat",skiprows=2,names=['e','cs','dcs'],delim_whitespace=True)
    >>> trans_ta1 = pd.read_csv(folder+"ta1_trans.dat",names=['e','t','dt'],delim_whitespace=True)
    >>> data1 = pd.concat([totxs_ta1,trans_ta1[['t','dt']]],axis=1)
    >>> 
    >>> totxs_ta3 = pd.read_csv(folder+"ta3_sig.dat",skiprows=2,names=['e','cs','dcs'],delim_whitespace=True)
    >>> trans_ta3 = pd.read_csv(folder+"ta3_trans.dat",names=['e','t','dt'],delim_whitespace=True)
    >>> data3 = pd.concat([totxs_ta3,trans_ta3[['t','dt']]],axis=1)
    >>> 
    >>> totxs_ta6 = pd.read_csv(folder+"ta6_sig.dat",skiprows=2,names=['e','cs','dcs'],delim_whitespace=True)
    >>> trans_ta6 = pd.read_csv(folder+"ta6_trans.dat",names=['e','t','dt'],delim_whitespace=True)
    >>> data6 = pd.concat([totxs_ta6,trans_ta6[['t','dt']]],axis=1)
    >>> 
    >>> capxs_ta1 = pd.read_csv(folder+"capxs_ta1.dat",skiprows=2,names=['e','cs','dcs'],delim_whitespace=True)
    >>> 
    >>> capxs_ta2 = pd.read_csv(folder+"capxs_ta2.dat",skiprows=2,names=['e','cs','dcs'],delim_whitespace=True)
    >>> 
    >>> data = [data1,data3,data6,capxs_ta1,capxs_ta2]

---------------
YAML input file
---------------

The yaml format input file contains many of the input variables needed to execute the function 
properly. An example of a  ``sesh_fitacs_inp.yml`` file is given below::

    ############################################################################
    #
    # This is an input file for the iterative operation of 
    # SESH and FITACS/SAMMY. It will be imported into the 
    # sesh_fitacs() function in nuctools.
    #
    ############################################################################

    # --------------------------------------------------------------------------
    # working directory, the fitacs and sesh in and out directories SHOULD BE 
    # IN THE WORKING DIRECTORY. Files will be opened from these directories as 
    # e.g. workdir+fitacs_indir
    # --------------------------------------------------------------------------
    workdir : /Users/jesse/Dropbox/ta_urr_fitting/

    # --------------------------------------------------------------------------
    # Boolean list for capture or trans. The order of this list must match the
    # order of the data list given to sesh_fitacs along with this input file.
    # --------------------------------------------------------------------------
    cap_bool : [False,False,False,True,True]

    # --------------------------------------------------------------------------
    # Sample thickness list for capture or trans [at/barn]
    # --------------------------------------------------------------------------
    samp_thick : [5.66e-3,1.713e-2,3.358e-2,5.631e-3,1.115e-2]

    # --------------------------------------------------------------------------
    # Sample file names that you wish to put corrected data into. (These will
    # be modified and recorded for each iteration)
    # --------------------------------------------------------------------------
    data_names : ['t1.dat','t3.dat','t6.dat','c1.dat','c2.dat']

    # --------------------------------------------------------------------------
    # This is the directory where the the corrected data files will be placed
    # --------------------------------------------------------------------------
    data_dir : corr_data/

    # --------------------------------------------------------------------------
    # This is the directory where the fitacs sesh lives
    # --------------------------------------------------------------------------
    fitacs_indir : fitacs_inp/

    # --------------------------------------------------------------------------
    # This is a directory that the fitacs output will be directed to
    # --------------------------------------------------------------------------
    fitacs_outdir : fitacs_out/

    # --------------------------------------------------------------------------
    # This is file name for the starting FITACS par file in fitacs_indir
    # --------------------------------------------------------------------------
    fitacs_par_name : ta181_urr_mult.par

    # --------------------------------------------------------------------------
    # This is the interactive input strings answering prompts from FITACS
    # --------------------------------------------------------------------------
    fitacs_int_input : fitacs_in

    # --------------------------------------------------------------------------
    # This is the directory that the sesh input files will reside in
    # --------------------------------------------------------------------------
    sesh_indir : sesh_inp/

    # --------------------------------------------------------------------------
    # This is the directory where the sesh output will be placed
    # --------------------------------------------------------------------------
    sesh_outdir : sesh_out/

    # --------------------------------------------------------------------------
    # This is a list of the file names for the base sesh input file. The 
    # parameters on lines 3 through 5 will be changed. (To inlcude pars for > L 
    # python source changes are required.)
    # 
    # -----------
    # -----------
    # - There should be a SESH input file for every sample (e.g. 1mm cap, 3mm trans..)
    # - The order of the inputs needs to match the order of the other lists
    # -----------
    # -----------
    # --------------------------------------------------------------------------
    sesh_ifile_name : ta_sesh.inp

    # --------------------------------------------------------------------------
    # This file lists the interactive input answering prompts by the RPI ver. 
    # of SESH
    # --------------------------------------------------------------------------
    sesh_int_input : sesh_in

    # --------------------------------------------------------------------------
    # This file is the correction factor file output from sesh (same as the 
    # one listed in sesh_int_input file.)
    # --------------------------------------------------------------------------
    sesh_cor : ta_sesh.cor

    # --------------------------------------------------------------------------
    # This file is the output from sesh (same as the one listed in sesh_int_input
    # file.)
    # --------------------------------------------------------------------------
    sesh_output : ta_sesh.out

    # --------------------------------------------------------------------------
    # The R, or the effective nuclear radius
    # --------------------------------------------------------------------------
    Rp : 7.8

    # --------------------------------------------------------------------------
    # The R, or the effective nuclear radius
    # --------------------------------------------------------------------------
    fitacs_numE_regions  : 3

    # --------------------------------------------------------------------------
    # The boundaries separating each of the energy regions. The number of bounds
    # should be one more than fitacs_numE_regions. e.g. if energy reg. 1 is 
    # 200-400 eV, energy reg. 2 is 400-600 eV, and energy reg. 3 600-800 eV then 
    # the boundaries should be: e_bounds = [200,400,600,800]
    # --------------------------------------------------------------------------
    e_bounds  : [2000,10000,45000,12000]

    # ----------------------
    # Calculate correction factor for first set of pars?
    # Often the beginning set has been calculated by a previous run.
    # 
    # If this option is used, the correction files must be named as expected by
    # the code (normally it is named for you.)
    # ----------------------
    calc_first_corr : False



**It should be noted** that the lists provided in the input file are also ordered in the
same order as the data that is being fit by FITACS. In this case, e.g., the sample 
thicknesses are listed as ``[5.66e-3,1.713e-2,3.358e-2,5.631e-3,1.115e-2]``. This order
corresponds to the 1 mm total cross section, 3 mm total, 6 mm total, 1 mm capture cross 
section, 2 mm capture.

-----------------------------------
Operational SESH & FITACS inp files
-----------------------------------

FITACS requires:

- input file
- par file
- data files
- interactive commands

FITACS must be capable of running these files for the ``sesh_fitacs()`` runner to work. The 
par file will be modified with different average resonance parameters throughout. The FITACS
code is a part of the SAMMY program. SAMMY gives prompts asking for each of the files it 
needs to run. To give the program the proper prompts, one can just pipe a file to the SAMMY 
program :: 

    sammy < interactive_cmd_file

This interactive command file looks like the example below (that last \\n is important)::


    fitacs_inp/ta181_urr.inp
    fitacs_inp/ta181_urr_mult.par
    /Users/jesse/Dropbox/ta_urr_fitting/corr_data/t1.dat
    /Users/jesse/Dropbox/ta_urr_fitting/corr_data/t3.dat
    /Users/jesse/Dropbox/ta_urr_fitting/corr_data/t6.dat
    /Users/jesse/Dropbox/ta_urr_fitting/corr_data/c1.dat
    /Users/jesse/Dropbox/ta_urr_fitting/corr_data/c2.dat



The RPI version of SESH also requires prompted inputs. The input can be given the same as sesh::

    sesh < interactive_cmd_file

An example of the command file is given below::


    "sesh_inp/ta_sesh.inp_e0_it2"
    "sesh_out/ta_sesh.out_e0_it2"
    "sesh_out/ta_sesh.ana_e0_it2"
    "sesh_out/ta_sesh.cor_e0_it2"



The input to SAMMY/FITACS is well documented in the SAMMY manual, but an example of a par file 
fitting multiple energy regions is provided for reference below::


    Ta-181 urr par file
    
    --------------------------
    ITERATIONS.=      2
    TOLERANCE. = 0.000005
    RADIUS    =     7.800
    AW.     =  180.947996
    
    --------------------------
    ELASTIC AND INELASTIC STATES
           0.0       3.5       1.0
        6237.0       4.5      -1.0
      136262.0       4.5       1.0
      158554.0       5.5      -1.0
      301622.0       5.5       1.0
      337540.0       6.5      -1.0
      482168.0       2.5       1.0
      495184.0       6.5       1.0
      542510.0       7.5      -1.0
      590060.0       3.5       1.0
      615190.0       0.5       1.0
      618990.0       1.5       1.0
      716659.0       7.5       1.0
      772970.0       8.5      -1.0
      892900.0       5.5       1.0
      965000.0       8.5       1.0
      994200.0       2.5      -1.0
     1022600.0       4.5      -1.0
     1028000.0       9.5      -1.0
     1085600.0       6.5       1.0
     1163600.0       6.5      -1.0
     1205700.0       1.5       1.0
     1239470.0       9.5       1.0
     1278100.0       2.5       1.0
     1304800.0       7.5       1.0
     1307110.0      10.5      -1.0
     
    --------------------------
    BINDING ENERGY (in MeV) = 7.57680000
    PAIRING ENERGY (in MeV) = 0.73000000
    
    --------------------------
    STRENGTH  DEL_S     DISTANT   DEL_D     GAM_WIDTH DEL_G     BETHED
     0.000185 0.0000110 -0.006600 0.0010000 0.0678000 0.0110000 4.1700000
     0.000050 0.0000200 0.0000000 0.0100000 0.0678000 0.0110000
     0.000230 0.0000300 0.0000000 0.0100000 0.0678000 0.0110000
     
    --------------------------
    MINIMUM ENERGY in eV = 2000.0
    ENERGY MAXIMUM in eV = 10000.0
    
    --------------------------
    BINDING ENERGY (in MeV) = 7.57680000
    PAIRING ENERGY (in MeV) = 0.73000000
    
    --------------------------
    STRENGTH  DEL_S     DISTANT   DEL_D     GAM_WIDTH DEL_G     BETHED
     0.000185 0.0000110 -0.006600 0.0010000 0.0678000 0.0110000 4.1700000
     0.000050 0.0000200 0.0000000 0.0100000 0.0678000 0.0110000
     0.000230 0.0000300 0.0000000 0.0100000 0.0678000 0.0110000
     
    --------------------------
    ENERGY MAXIMUM in MeV = 0.045
    
    --------------------------
    BINDING ENERGY (in MeV) = 7.57680000
    PAIRING ENERGY (in MeV) = 0.73000000
    
    --------------------------
    STRENGTH  DEL_S     DISTANT   DEL_D     GAM_WIDTH DEL_G     BETHED
     0.000185 0.0000110 -0.006600 0.0010000 0.0678000 0.0110000 4.1700000
     0.000050 0.0000200 0.0000000 0.0100000 0.0678000 0.0110000
     0.000230 0.0000300 0.0000000 0.0100000 0.0678000 0.0110000
     
    --------------------------
    ENERGY MAXIMUM in MeV = 0.120
    
    END OF RESONANCE PARAMETER DESCRIPTION
    --------------------------
    NORMALIZATIONS 
    TOTAL       1.000000  0.000000  0.000000  0.000000  0.000000  0.000000
    TOTAL       1.000000  0.000000  0.000000  0.000000  0.000000  0.000000
    TOTAL       1.000000  0.000000  0.000000  0.000000  0.000000  0.000000
    CAPTURE     1.000000  0.000000  0.000000  0.000000  0.000000  0.000000
    CAPTURE     1.000000  0.000000  0.000000  0.000000  0.000000  0.000000
    



The SESH manual contains crucial information on SESH, and is a good reference for the
theory basis of the code, but has sparse information on how to run the code. A supplemental
SESH manual may be included in this documentation at a later date. An example of a SESH
input file for 3 partial waves and 5 samples of both transmission and capture yield is 
given below for reference::


    2mm Ta-181 Multiple Scattering & Self Shielding Correction                0
    181.0     1.0       7.57680   0.730     294.      3.50
      0.06780   4.17000 1.850e-04   0.00000  7.80000  1.00000
      0.06780   4.17000 5.000e-05   0.00000  7.80000  1.00000
      0.06780   4.17000 2.300e-04   0.00000  7.80000  1.00000
              5.660E-03 1.713E-02 3.358E-02 5.631E-03 1.115E-02
              0.000E+00 0.000E+00 0.000E+00 3.351E-01 3.351E-01
              0.000E+00 0.000E+00 0.000E+00 0.000E+00 0.000E+00
              1.377E-01 1.377E-01 1.377E-01 1.377E-01 1.377E-01
    17.
    1.50000   15000.    1.75000   15000.    2.00000   15000.    2.25000   15000.
    2.50000   15000.    2.75000   15000.    3.00000   15000.    3.25000   15000.
    3.50000   15000.    4.00000   15000.    4.50000   15000.    5.00000   15000.
    5.50000   15000.    6.00000   15000.    6.50000   15000.    7.00000   15000.
    7.50000   15000.    8.00000   15000.    8.50000   15000.    9.00000   15000.
    9.50000   15000.    10.0000   15000.    10.5000   15000.    11.0000   15000.
    11.5000   15000.    12.0000   15000.    12.5000   15000.    13.0000   15000.
    14.0000   15000.    15.0000   15000.    16.0000   15000.    17.0000   15000.
    18.0000   15000.    19.0000   15000.    20.0000   15000.    21.0000   15000.
    22.0000   15000.    23.0000   15000.    24.0000   15000.    25.0000   15000.
    26.0000   15000.    28.0000   15000.    30.0000   15000.    32.0000   15000.
    36.0000   15000.    40.0000   15000.    44.0000   15000.    48.0000   15000.
    52.0000   15000.    56.0000   15000.    60.0000   15000.    70.0000   15000.
    80.0000   15000.    90.0000   15000.    100.000   15000.    110.000   15000.
    120.000   15000.    130.000   15000.    140.000   15000.    150.000   15000.
    
    
    
    
    
    
    
    
    
    
    


**There should be 26 lines in the input file following the resnonace pair description.** The
resonance pair descriptor is 17 here, seen at line 10. It is also important to note that if 
the number of resonance pairs (17 here) is low (such as 5 or 6) the **SESH program may enter** 
**an infinite loop and never report a problem.** The greater the number of resonance pairs used
the more likely the code is to complete, and the longer it will take the code to complete. This
is an unresolved bug in the code.

--------------------------
SESH and SAMMY executables
--------------------------

The ``sammy`` program can be obtained from the `ORNL sammy website <https://code.ornl.gov/RNSD/SAMMY>`_.
This open-source ``sammy`` program requires SCALE file to compile until the open-source version of 
AMPX is available. The ``sesh`` program can be obtained from
 `Jesses GitHub page <https://github.com/brownjm1968/sesh>`_. ``sesh`` is a stand-alone Fortran code
 and does not require any supporting software. 
 





















