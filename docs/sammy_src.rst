============
Sammy Source
============

---
acs
---

The acs folder contains Fritz Frohner's FITACS code. ``macs0.f`` is
the main program which calls subroutines from ``macs1.f`` through
``macs7.f``. Subroutines include:

macs0.f

* Samacs_0 - main
* Write30

macs1.f

* Begin_Iteration
* Calthe - call the routines to calc. theor. x-sec's and deriv's
* Elastic_Urr - calc. elast. x-sec (total - others)

macs2.f

* Caltot - Calc. theor. tot. x-sec and deriv's
* Find_Urr

macs3.f

* Calpar - Calc. partial x-sec's and their deriv's for given (inel,fis,cap)
* Wav - Computes the width fluct. average and it's deriv's

macs4.f

* Zelast - Elastic called in Calpar
* Zinela - Inelastic called in Calpar
* Zfissi - Fission called in Calpar
* Deriva - Calculates derivatives for partial x-sec's
* Endep - Calculate energy depend. of level spacings and radiation widths
* Find_Www_Yyy 
* Simpx - Simpson's rule
* Eff_Dof - Calculate effective degrees of freedom

macs5.f

* Acsout - prints given and calculated avg. x-sec's
* Plott2 - Writes calculated value plot files
* Make_File_3_Urr

macs6.f

* Define_Constants
* Kwhich

macs7.f

* Dres - Calculate the Dresner Integral

---
amr
---

The amr folder contains the code to change the \"active\" parameter list
in order to accomadate for different experimental situations. Parameters
that may be treated include: broadening parameters (TEMP,THICK,DELL,..)
MISCellaneous parameters, Oak Ridge Resolution function parameters, 
BacKground parameters (ANORM,BACKx) and more.

mamr.f 

* Samamr - main

mamr1.f

* Startx - organize the results

mamr2.f

* Array - fix array sizes
* Kdimen - Keep track of dimensions for character array
* Keep1
* Keepx
* Keep2

mamr3.f

* Ask1 - Learn FILe NAMes and Action to be taken
* Ask2 - Ask questions for what par to change
* Addmix - determine if Action = \'adv\'
* Newrun - Get information about the replacement parameter
* Xx - Par. was not varied before, is now, what will it's number be?
* X - Find number for added varied parameter
* Mixmix - Decide which unused parameter is wanted
* Intkbn 
* Fix_Anorm
* Recovr
* Rematt
* Intatt - Get attenuation factor file

mamr4.f 

* Covfix - Fix up the covariance matrix
* Fixbrd 
* Fixnbk 
* Fixall
* Fixbgf_Old
* Fixorr_Old
* Fixjfl 
* Fixxxx
* Remkbn
* Remat2
* Recatt
* Reckbn
* Intat2
* Intkb2

mamr5.f

* Fixcov
* Fxcadv
* Fxcadd
* Fxcmix - Action=\'mix\'...interchanging two variables
* Fxcrem - Action=\'rem\'...to REMOVE a set of par's without adding any others
* Fxcrec - Action=\'rec\'...to RECOVER a set of par's from UNUSED to ACTIVE
* Fxcint - Action=\'int\'...to introduce a new set of par's
* Fixu
* Fxuadv - Intro. a variable into middle of list 
* Fxuint
* Fxurem - Move Nvp var's from site \'Jold\' to end of list
* Fxurec - Move Nvp var's from site \'Jnew\' to site \'Jold\'

mamr6.f

* Parfix - Calls all the parameter fix functions
* Rdwrtx
* Rdwrt
* Rdwrt_38
* Cblnk
* Rd
* Brdfix
* Mscfix 
* Pmcfix
* Orrfix
* Rpifix
* Nbkfix
* Wrtusd
* Wrtbag - write BAGGAGE cards

mamr7.f

* Wbrd
* Rchn
* Wmsc
* Worr
* Wrpi
* Wnbk
* Wbgf
* Watt
* Bgffix - background function

mamr8.f

* Iffy
* Iffy_old
* Bbbbrd
* Dddist
* Telbrd
* Fixiso

mamr9.f

* Mthbll - learn which parameter is to be mothballed


---
amx
---

Program to change an unvaried paramter to a different value. Use of 
this program is discouraged

---
anf
---

Program writes data from two ODF files for angular distribution data 
into one file. 

manf.f

* Angodf - main
* Finde

---
ang
---

Used to generate the angle-differential cross section at angles 
"Angle", starting form the Legendre coeff's, and the deriv's 
thereof. Also does incident and outgoing neutron attenuation.

mang0.f

* Samang_0 - main
* Estang

mang1.f

* Diffee - generates the cross section's and deriv's
* Polynm - generate Legendre Poly's at the required angles
* Polynm2 - same as Polynm, also generates fudge factor for conv. to scat. energy
* Mult 
* Write_Legendre
* Include_Coulomb
* Only_Norm_Ang
* Reorg_Iso_Pq

------
angodf
------

---
avg
---

Average over energy ranges to give publishable results.

mavg0.f

* Samavg_0 - main
* Estavg - guess size of array needed

mavg1.f

* Tryrng
* Qrange

mavg2.f

* Qprint - Output A values
* Both
* Write_Std_and_Corr - print triang. variance as std. dev. and corr.

mavg3.f

* Test_If_Data - Determine if there are more data
* Again - get ready for another data set

mavg4.f

* Getdel - Decide which averaging method to use
* Engavg - de Saussure's method (bins)
* Timeav - time averaged method?
* Testee
* Testmm
* Testdd
* Xnorm

mavg5.f

* Gen_Avg_Data - avg data to give Datq
* Datfix - Rebin Vardat to Varq

mavg6.f

* Gfix - average partial derivatives G to Gq
* Thefix - multiply ``Emmmq = gq * vrpr * gq``
* Resetgx - Debug printout of partial deriv's

---
blk
---

The blk folder contains files that hold variables, and are called throughout the code using INCLUDE statements

---
ccm
---

Calculates the covariance matrix for cross section.

mccm0.f

* Samccm_0 - main
* Estccm - estimate the size needed for segment CCM

mccm1.f 

* Create_Cov - generate covariance matrix for the theoretical x-sec

---
clm
---

Calculates the absorption cross section by affecting the nuclear cross
section by the phonon distribution probability. Choice is between free
gas and harmonic crystal models. Program is DOPUSH.

---
clq
---

Used to generate cross section that is constant, linear, quadratic,
Dirac-delta, or 1/v for use in debugging the resolution function.

-----
cmake
-----

Contains dependencies of cmake

---
cou
---

Coulomb version of sub pgh. Employs subroutines to calculate the 
coulomb wave function.

---
cpr
---

Read SAM53.dat (which was generated by samnpv or ipq) and create odf
file.

---
cro
---

Contains many functions to control the R-matrix, and generate cross
section and convert total cross section to transmission

---
dat
---

Reads in the experimental data, and determine various coding
requirements such as the energy grid used to calculate the 
theoretical cross section.

---
dbd
---

Perform high-energy gaussian approximation version of Doppler
broadening

---
dex
---

Forms the resolution broadened cross section and derivatives for the
dE/dx resolution

---
dis
---

Samdist program written by Luiz Leal

---
dop
---

Implements Doppler broadening

mdop0.f 

* Samdop_0 - main

mdop1.f

* Dopplh - Implementation of Doppler broadening
* Set_Negative_E - Deal with storage of negative E
* Multiply_By_Coef
* Temp_Deriv
* Transform_to_Cross 

mdop2.f

* Coefgn
* Findst

-------
dopushx
-------

Compilation tools for DOPUSH

---
end
---

Bookkeeping for the master SAMMY program

mend0.f

* Samend_0 - main
* Ridaaa - clean up files before quitting or restarting

mend1.f

* Readxx

mout.f

* Outpar - writes res par's into SAMMY.LPT
* Advzer 
* Header
* Nwrite
* Headpc
* Pwrite
* Averag
* Oooppp - write the res par's when there are more than 3 channels

mout1.f

* Setflg
* Outred - output resonance par's (amplitudes)
* Outre2 - output reduced widths (not amplitudes)
* Outvs - output triangular variance V as std dev plus corr.
* Outvr - output cov matrix for reduced par's
* Outstr
* Outmlb - output equivalent Breit-Wigner par'

mout2.f 

* Outext - output R-external par's
* Outrad - output radii
* Outiso - output isotopic information
* Outdet - output detector efficiencies
* Outbrd - output \'broadening\' par's
* Outmsc - output the miscellaneous par's
* Outpmc - output par's for paramagnetic cross section

mout3.f

* Outorr - output ORELA resolution par's
* Outrpi - output ORELA resolution par's
* Outnbk - output normalization and bkg par's
* Outbgf - output bkg functions into SAMMY.LPT
* Outdtp - output \'data\' par's
* Outusd - output \'unused\' par's
* Outbag - output \'baggage\' par's

mout4.f

* Getvvv - generate the sub - cov. matrix for each width for a spin grp
* Invert - inverting matrices for Getvvv

mout5.f

* Parout - prints cross section par's for URR

mout6.f

* User_Output - write the res par's in user-defined format

msamvv.f

* Initix - interchange names of Filein and Filout
* Initil - set (N,K,J,L)size = (M,K,J,L)size
* Write_Commons_Few 
* Write_Commons_Many
* Write_Commons

msamxx.f - routines here may be machine dependent

* Oldopn - open existing file with 10 char name (SAMxx.DAT)
* Newopn - open a new file with a 10 char name
* Filopn - open a file with a 70 char name (user's files)
* X17opn - open a file with a 17 char name
* Timer 
* Timer_Debug

msamyy.f 

* Idimen - keep track of dimensions for end subroutines
* Fix_Alpha

neuopn.f

* NEUOPEN

---
fdc
---

Calls routines to modify the experimental data, such as normalization
and background subtraction

mfdc0.f

* Samfdc_0 - update the data file via D => (D-b)/a, a=norm, b=bkg
* Estfdc - estimate the size of array needed for fdc

mfdc1.f

* Readwr - call routines to read data, modify values, and rewrite
* Rddat1_Fdc - read E,xs/T,unc. modify. rewrite. format = 3G11.8
* Rddat0_Fdc - read E,xs/T,unc. modify via norm and bkg. rewrite. 3(2e15.8,f7.5)

mfdc2.f

* Stndrd_Fdc - read from standard ODF data file for E,data, diagonal unc
* Anglrd_Fdc - read from \'sort-of STANDARD ODF DATA FILE\' ?
* Angxrd_Fdc - read from \'sort-of STANDARD ODF DATA FILE\' ?

mfdc3.f

* Fixfdc - generate the revised version of data and unc

mfdc4.f

* Paramf - read and copy the parameter file
* Wheren - find final non-blank character in a line

---
fff
---

fff is the input for the implementation of Fritz Frohner's FITACS code in SAMMY

mfff0.f

* Samfff_0 - main
* Begin
* Wrt49
* Adjust_Fff
* Find_Limits - calc. of part. x-sec's and deriv's for given (inel,fiss,cap)

mfff1.f

* Pass - a first pass through FITACS file to get dimensions for array sizing
* Passxx - reorganize now that we've counted the number of data sets
* Passd - read unc. type and cross section data for Pass

mfff2.f

* Parin - read x-sec par's and overwrite w/ info from cov. mat. as needed?

mfff3.f

* Lespac - calc. level spacings for all relevant L and J
* Read_Zta 
* Read_Bp 

mfff4.f

* Annotation
* Annotation_Units
* Units
* Report_Annot

mfff5.f

* Copy_Urr - make a copy of par's from one energy region to use as starting 
values in another region

mfff6.f

* Read_Cov_Urr - read previously existing covariance file

mfff7.f

* Xwrong
* Ywrong
* Zwrong
* Wwrong
* Uwrong

mfff8.f

* Read_Direct - read direct inelastic or cap x-sec values

mfff9.f 

* Acsinp 
* Acsin - reads and prints avg x-sec exp data and checks energy order
* Plott1 - prepare a forodf plot file 
* Setdim

---
fgm
---

Calculates Doppler broadening using the free gas model with auxiliary energy 
grid of varying spacing

mfgm0.f 

* Samfgm_0 - main
* Estfgm - estimate size of array needed for fgm

mfgm1.f

* Velfix 
* Dopfgm - form Doppler broad. x-sec and deriv's

mfgm2.f 

* Modsmp - form the Gaussian Doppler weights via Simpson's rule (sort of)
* Resets
* Modfpl - form the Doppler weights via modified four-point Lagrange rule
* Reset
* Wtxxd - Start integrals for different conditions
* Wtxxcd - Start integrals for different conditions
* Wtxbcd - Start integrals for different conditions
* Wtabcd - Start integrals for different conditions
* Wtabcx - Start integrals for different conditions
* Wtabxx - Start integrals for different conditions
* Wtaxxx - Start integrals for different conditions
* Start4 
* Start3 
* Start2 
* Start1 
* WtZxxd - Quit integrals for different conditions
* WtZxcd - Quit integrals for different conditions
* WtZbcd - Quit integrals for different conditions
* WtabcZ - Quit integrals for different conditions
* WtabxZ - Quit integrals for different conditions
* WtaxxZ - Quit integrals for different conditions
* Quit4
* Quit3
* Quit2
* Quit1

mfgm3.f

* Kountd - Find Kc and Iup, and Ipnts
* Which
* Stetd - copy data and derivatives b/c there are too few pts to broaden

mfgm4.f

* Xdofgm - perform integration for Dopp. broadening from Kc to Iup
* Funfgm

mfgm5.f

* Add_Nuclide - add the nuclides...sum over isotopes
* Store_S - store in Sigxxx and Dasigx and Dbsigx temporarily
* Store_Nuclide - Copy from Sigxxx,Dasigx,Dbsigx to W-arrays
* Store_W
* Store_W_Iso

---
fin
---

Purpose is to convert to physical parameters and output results

---
fit
---

Purpose is to convert to physical parameters and output results when 
operating FITACS portion of SAMMY

---
fnc
---

Contains functions (log, exponential, erfc...etc)

---
ftz
---

Samftz program to \'fix\' TZero

---
grp
---

Generates group averages (Bondarenko and others)

---
gy2
---

Reads the MC_2.DAT file from samsmc and smooths the Y2 values

---
idc
---

Reads and sorts the implicit data covariance information

---
inp
---

Read the INPut file and get information organized

minp0.f

* Saminp_0 - main

minp01.f

* Files - read file names, open files
* Estinp
* Set_Dimensions
* Filesx - read SSM (self-shielding) and extra file names

minp02.f

* Inp2 - read input file to learn how many groups and other info..
* Set_Kip
* Firstr - Determine values by reading through the PAR file
* Wrongf
* Wrongi

minp03.f

* Inpfil - read input file for control messages, non-variable par's, quantum numbers
    
    * Card Set 11 
        
        * Pread
        * 
* Huhx

minp04.f

* Huh - read input verbage (card set 3 in input file)
* Setdef - set default values for flags to be set in Huh

minp05.f

* Fixdef - adjust values of flags
* Set_Spninc

minp06.f

* Rdone_0 - read first 2 lines of input file, write it out
* Rdone_1 - re-read first 2 lines of input file
* Rdbrd - read broadening information
* Nrdbrd - assign values to broadening parameters when they are not read
* Fissil - set Nfissl=0 if non-Fissile =1 if Fissile
* Ssm_Read - read card set 11 but ignore b/c already read
* Lmaxxx - modified from subroutine Clbsch to give more accur. answer

minp07.f

* Rd_Pp - read the particle pair definitions
* Rd_Pp_Key_Word - Card set 4
* Report_3 - write particle pair information onto unit 21

minp09.f

* Set_Nam - extract alphanumeric \'name\' from large array \'A\'
* Get_Parity - find plus or minus sign
* Get_Ps - extract yes or no (Jppairx=1,0)
* Get_Particle - determine part. name for one of the part. pair
* Which_Particle 
* Set_Wp 
* Get_Mass
* Get_Charge
* Get_Spin
* Get_Value - extract assigned val. for \'Value\' following equal sign
* Find_Iab

minp10.f

* Adjust - adjust constants
* Wrcros - interpret info about type of x-sec's
* Only - learn about type of x-sec if there is only one
* Zernst

minp11.f

* Rdspin - read res quantum numbers and other info
* Testex
* Test_Emmms

minp12.f

* Report_1 - report quantum number info when Card set 10.1 is used

minp13.f

* Set_Dum_Ver
* Set_Nam_Pp
* Final_States - determine which part. pair numbers are to be included in final
state
* Copy_Nam_Pp_8 - reorganize
* Test_Angula

minp14.f

* Rdspno - read res quantum numbers and other info using "ORIGINAL SPIN FORMAT"
* Set_Dum_Ver_X
* Rdspnx - Write "Card set 10.1", i.e. new spin group format

minp15.f

* Rdspi2 - read resonance quantum numbers and other info using card set 10.2, pp 
def's
* Get_Lab_Cm_Conversion - Find fudge factors for conversion of center-of-mass to
laboratory for angular distribution
* Reorg_Cmlab - reorganize storage of Cmlab and Iso_Qv

minp16.f 

* Report_2 - report quantum number info etc when part. pair def's are given

minp17.f

* Organize_Bound_Etc - generate Goj and Bound and other things

minp18.f

* Qresp 
* Qfudge
* Qextr - read external R-function parameters (Card set 3 in PAR file)
* Qrext - read R-external parameters (altern. to card set 3 in PAR file)
* Qradi - \'RADII PARAMETERS ARE ON THE NEXT CARDS\'
* Find_Key_Word
* Qisot - isotopic abundance and mass
* Qdetec - detector efficiencies
* Qbroa - \'BROADENING PARAMETERS MAY BE VARIED\'
* Qmisc - Miscellaneous parameters are given here

minp19.f

* Qpmcs - paramagnetic x-sec parameters are given
* Qorre - '\ORRESOLUTION\' Oak ridge resolution function Card set 9
* Qnorm - \'NORMALIZATION AND CONSTANT BACKGROUND\' Card set 6
* Qbackg - background functions re RPI specifications
* Qdatp - \'DATA PARAMETERS ARE GIVEN HERE\'
* Qunus - \'UNUSED BUT CORRELATED VARIABLES\'
* Qbaga - baggage par's
* Qcova - \'COVARIANCE MATRIX IS IN BINARY FORM\'
* Qrelu - \'RELATIVE UNCERTAINTIES\'
* Qexpl - explicit covariance matrix elements

minp20.f

* Qrpire - RPI resolution function
* Qrpitc - RPI res. func.; transmission or capture defaults
* Qgeeln - Geel or nTOF resolution function

minp21.f

* Quserd - user-defined resolution func. 

minp25.f

* Pread - reads message line from .par or .inp
* Set_Talk

minp26.f

* Gg1234 - skim through Card sets 1 and 2 of INPut file. read card set 3 and 4
* Gg4xxx - Read card set 4.5 of INPut file - energies to plot res. func.
* Gggg56 - read card sets 5 and 6 of INPut file to learn Ncf
* Gggg78 - read card sets 7 and 8 of INPut file
* Gg10p0 - read card sets 9 and 10 of INPut file
* Gg10p1 - read card set 10.1 of INPut file
* Gg10p2 - read card set 10.2 of INPut file
* Gggg11 - read card set 11 of INPut file

minp27.f

* Last_Character
* Convert_To_Caps

---
int
---

Code for output of cross sections

mint0.f

* Samint_0 - main
* Estint - estimate array size for int
* Intermediate - prepare to write even though there's more broad. and others
* Finished - no more broadening, result ready to write

mint1.f

* Interp - interp. between equally spaced pts to experimental grid (Lagrange)
* Mix_V_to_W 
* Flip
* Flip2

mint2.f

* Outthr - print theory as a func. of energy (not for angular output)
* Outtwo - print theory as a func. of energy when Outthr doesn't work
* Outg
* Out_Deriv - output derivatives as a func. of energy
* Chzero 
* Wrxxx
* Wascii - write theor. values into ascii file SAMTHE.dat
* Set_Title

mint3.f

* Thodf - fill section 4,5 (and maybe 8,9) of \'regular\' od file
* Undo 
* Plotun - generate plot file for unbroadened theor. x-sec

mint4.f

* Pdwrit - write the partial derivatives

mint5.f 

* Outddd - output theory as a function of energy, in cols, for Na > 1

mint6.f

* Leal_Hwang - write x-sec's and part. deriv's, interpolated onto 
Leal-Hwang energy grid

---
ipq
---

mipq0.f

* Samipq_0 - solve Bayes' equations to obtain results using I+Q inversion scheme
* Estipq - Estimate the array sized needed for IPQ
* Estx 

mipq1.f

* Newpar_Ipq - solve Bayes' equations via (I+Q) inversion scheme
* Stqp1 - multiply G * V^-1 * G * M, set P1 = G * V^-1 * (D-T-..) - G * V^-1 * X * Wpxvxi * X * V^-1 * (D-T-..)

mipq2.f

* Mul 
* Chise2 - generate chi-squared and Bayesian weighted residuals
* Chngx2 - generate updated values of par's: parm is updated
* Chngy2 - update par values via method of steepest descent
* Chngc2 - update cov. mat. generate change in par. val's: Delpar, Em, Vrpr
* Add - add a=b+c (new par = old par + del par)
* Stpadd - find new par for method of steepest descent

mipq3.f

* Fff111 - aprox. inv. of a real mat. found by solveing AX=I using Crout's method
* Fff333 - ``THE UNSYMMETRIC MATRIX, A, IS STORED IN THE N*N ARRAY A(I,J), I=1,N, J=1,N. THE DECOMPOSITION A=LU, WHERE L IS A LOWER TRIANGULAR MATRIX AND U IS A UNIT UPPER TRIANGULAR MATRIX, IS PERFORMED AND OVERWRITTEN ON A, OMITTING THE UNIT DIAGONAL OF U. A RECORD OF ANY INTERCHANGES MADE TO THE ROWS OF A IS KEPT IN P(I), I=1,N, SUCH THAT THE I-TH ROW AND THE P(I)-TH ROW WERE INTERCHANGED AT THE I-TH STEP. THE SUBROUTINE WILL FAIL IF A, MODIFIED BY THE ROUNDING ERRORS, IS SINGULAR OR ALMOST SINGULAR. SETS Ifail = 0 IF SUCCESSFUL ELSE Ifail = 1.``
* Fff444 - SOLVES AX=B, WHERE A IS AN UNSYMMETRIC MATRIX AND B IS AN N*IR MATRIX OF IR RIGHT-HAND SIDES. THE SUBROUTINE Fff444 MUST BY PRECEDED BY Fff333 IN WHICH L AND U ARE PRODUCED IN A(I,J), FROM A, AND THE RECORD OF THE INTERCHANGES IS PRODUCED IN P(I). AX=B IS SOLVED IN THREE STEPS, INTERCHANGE THE ELEMENTS OF B, LY=B AND UX=Y. THE MATRICES Y AND THEN X ARE OVERWRITTEN ON B.
* Pp1111 - returns the value of error or terminates the program
* Xxx333 - calc. the val. of a scalar product using basic or additional precision and adds it to a basic or additional precision initial value
* Xxx444 - changes the error message for Nerr

mipq4.f

* Print_Array - output I vs A1(I) in columns
* Setfor 
* Readu
* Readiu

---
lru
---

mlru0.f

* Samlru_0 - write ENDF file 2 and 32 for the URR

mlru1.f

* Urr
* Read_Eee
* Read_Ndf_Urr
* Write_Urr_1
* Write_Urr_2
* Write_Urr_3
* Write_Urr_4

mlru2.f

* Parout_Endf - find and print par's into ENDF file
* Zurr
* Urr_Dd - generate level density and uncertainty
* Urr_Gn - generate values and uncertainties for elastic width
* Urr_Gg_Cov - generate uncertainties for gamma width
* Urr_Gf - generate uncertainties for fission width
* Urr_Gx - generate uncertainties for inelastic width
* Write_Urr_5

mlru3.f

* Read_Covar - Read Vrpr = covariance mat. for SAMMY par's
* Fix_Covar - Determ. covar. = cov. mat. for ENDF par's
* Write_Cmpc_Cov_2 - write cov. mat. in compact form
* Ijxxx
* Write_Urr_Lcomp_1 - write cov. mat. in Lcomp="undef." format

---
mas
---

Master control program for running SAMMY

mmas0.f

* Sammas_0 - main program

mmas1.f

* Inppar - takes user input for inp and par file
* Finpx - reads the input file for: Kywywy, Kwywyw, Krdspn, Kpntws, Ndfinp, Kompci, Nretro, Kedepu, Mcy2
* Fpar_If_Cov - read parameter file to determine if covariace file exists
* Finp - read INPut file, obtain info needed to create other INP files for each pass
* Zeroxx 
* Whatis 
* Sumstr - prepare to evaluate summed strengths
* Get_Stop_Segment - determine the segment before which the code should write files and stop for debugging

mmas2.f

* Fixnam_2 - Read line from Unit 5, return two values (2E) and Filnam
* Fixnam - Read line from Unit 5, return Filnam
* Whaten
* Odfdat - Write data from ODF file onto dummy file named Fdatax for use in samdat segment
* Stdodf - Read standard data file to be sure data exist within limits given
* Ifind

mmas3.f

* Datacov - Learn DATa file names and energy ranges
* Fdat - Check whether data in the specified range are in the DATa file
* F3dat - Read ENDF data file to be sure there are data within limits given
* Ascidt - Read ASCII data file to be sure there are data within limits given
* Ascxdt - Read ASCII data file for differential data, to be sure there are data within limits given
* Fdat1 - Divide into regions if needed
* Parcov - Learn name of par cov. file, create PAR file if needed

mmas4.f

* Wrt16x - open and write SAM16.DAT, which contains name of FITACS input file
* Writ16 - open and write SAM16.DAT, which is the equivalent of the TTY file for the first and second pass thru SAMMY
* Writ19 - set up file (SAM19) to be used as TTY input the third pass thru SAMMY
* File2x - find out whether there exist further data

mmas5.f 

* New_Input_File - look at INPut file, create new INPut files to be read in SAMMY-INP on the three passes
* Sett
* Set
* Finis1 
* Write_Alpha

mmas6.f

* Endf1 - copy most of first part of the INPut file to use ENDF-file info along with user-supplied information
* Newinp - read alphanumeric information from original input file and write new input file
* Mina
* Endf2 - Create rest of INPut file using ENDF-file info, also make PARameter file using ENDF-file info
* Quantu - generate all possible combinations of l,s,J for given I
* Smlbw 
* Rmoore
* Rcontx
* Rcont
* Rtab1
* Rlist
* Findsp - find group number (Is) for given Llll and Spinjj (Spinjj is ENDF value Aj)
* Fixjvl - find "correct" value of aJ if J-value in ENDF/B is not known (ENDF sets J=I+L)
* Newin2 - write remainder of INPut file
* Wrtrad - write radii at end of parameter file
* Read_Cov_Mpar - open and read ENDF covariance file to see how many parameters per resonance are to be flagged

mmas6a.f

* Read_Pp7 - read and write particle-pair definitions
* Endf7 - Make rest of INPut file using ENDF-file information, also make rest of PAR file

mmas7.f

* Setcns - set the constant values (e.g. mass of neutron, speed of light, etc.)

mmas9.f

* Rewrite_12 - put paramter file into better format
* Fix_30 - reduce 30 characters to Num_Col (10 or 11)
* Asterisk_30 - reduce 12 characters to 11 
* Blankx - blank characters, replace Zero
* Fix_Extension - put "extension" command into Beta
* Copy_X - copy Kountr characters from Alpha to Beta
* Round_Alpha - round up the numbers if needed
* Special_Case - test for special cases where Jfirst must be increased by 1
* Kase 
* Kasex
* Values

mmasa.f

* Find_Matnum_In_ENDF_File - read thru the ENDF file until Matnum, Ifile, and Mtx are found
* 


---
ssm
---


Self-Shielding and multiple scattering corrections. This version contains both single
and double scattering as well as self-shielding

m012.f

* Fix_Sam012 - read ascii file (pieces of multiple-scattering correction) and create comparable odf file

mssm00.f

* Samssm_0 - main program
* Estss1 
* Estssm - estimate the size of array needed for samSSM
* Set_Logic_Ssm - e.g. linear interp. or quadratic interp.
* Ssm_Get_Organized 

mssm01.f

* Qqqxxx - read the dimensions of Xtpt_V and Xtpt_W etc, both are used for interpolation on Sqfb(V,W,mu), where mu = cos(theta)
* Qqqyyy - read the Xtpt's 
* Thtget - initialize array Ftheta, when there are no edge corrections
* Qqqget - Read the arrays Ftheta, and Sqfb. Also generate the rest of Ftheta
* Sam_Initialize 
* Getcrs - find the next cross section, store values in appropriate places
* Getem - determine where we are, energy-wise; i.e. which pieces of the calculation need to happen now
* Get_Angles - define angle grids

mssm02.f

* Ssssds - generates the self-shielded + multiple scattered capture yield for an infinite slab, and calls Mulsca to generate the finite-slab results
* Zero0_1f
* Zero0_2i
* Zero0_2f
* Tell_Finite
* X_Trpths_Lin
* X_Trpths_Quad

mssm03.f 

* Ssssds_0x - generates the self-shielded capture yield with no scattering corrections
* Ssssds_1il - generates the self-shielded capture with single scatter yield using linear interpolation (infinite slab)
* Ssssds_1iq - generates the self-shielded capture with single scatter yield using quadratic interpolation (infinite slab)
* Ssssds_1fl - generates the self-shielded capture with single scatter yield using linear interpolation (finite slab)
* Ssssds_1fq - generates the self-shielded capture with single scatter yield using quadratic interpolation (finite slab)

mssm04.f

* Ssssds_2il - generates the self-shielded capture with multiple scatter yield using linear interpolation (infinite slab)
* Ssssds_2iq - generates the self-shielded capture with multiple scatter yield using quadratic interpolation (infinite slab)
* Ssssds_2fl - generates the self-shielded capture with multiple scatter yield using linear interpolation (finite slab)
* Ssssds_2fq - generates the self-shielded capture with multiple scatter yield using quadratic interpolation (finite slab)

mssm05.f - functions to apply the multiple scattering correction

mssm06.f - functions to apply the multiple scattering correction

mssm07.f - functions to apply the multiple scattering correction

mssm08.f - functions to apply the multiple scattering correction

mssm19.f

* Ynrm_0 - normalizes to the yield for no scattering yield, and stores the yield which is output to plot files in sigxxx
* Ynrm_0_Derivative - derivative for the case of Ynrm_0

mssm20.f

* Fixy_1i
* Fixy_1f
* Non_Uniform_Thickness
* Get_Ratio_Sensin
* Ynrm_1 - normalizes to the yield for singly scattered yield, and stores the yield which is output to plot files in sigxxx
* Ynrm_1_Derivative - derivative for the case of Ynrm_1
* Selfin - self indication ratio

mssm21.f

* Fixy_2
* Ynrm_2 - normalizes to the yield for multiply scattered yield, and stores the yield which is output to plot files in sigxxx
* Ynrm_2_Derivative - derivative for the case of Ynrm_2
* Gamma_attenuation_corr - applies a correction to the yield (primary, and singly,multiply scattered) to account for gamma attenuation in the capture sample

mssm22.f

* Finish_01
* Finish_02

---
sta
---






























