import nuctools.urr_tools as ut
import numpy as np
import pandas as pd
import os,sys

def test_apply_sesh_corr():

    # ------- Total xs correction -------------------------------------------
    directory = os.getcwd() + '/tests/test_files/'

    totxs_ta1 = pd.read_csv(directory+"test_tot_xs.dat",skiprows=2,
                            names=['e','cs','dcs'],sep=r'\s+')
    trans_ta1 = pd.read_csv(directory+"test_trans.dat",names=['e','t','dT'],sep=r'\s+')

    N = 0.1

    new_xs,new_dxs = ut.apply_sesh_corr(directory+'test_sscf.dat',
                                         [totxs_ta1.e,trans_ta1.t,trans_ta1.dT],
                                         N,capture=False)

    answer_xs = np.repeat(-1/N*np.log(0.5),3)
    answer_dxs = np.repeat(1/N*np.sqrt((0.1/0.5)**2+(0.1/100)**2),3)

    np.testing.assert_array_equal(answer_xs,new_xs)
    np.testing.assert_array_equal(answer_dxs,new_dxs)
    # -----------------------------------------------------------------------

def test_read_fitacs_par():

    

    test = np.array([
        [0.00015753,0.00000313,-0.0065603, 0.0009997, 0.0629051, 0.0017586, 4.17],
        [0.00007308,0.00002404, 0.0000197, 0.0100000, 0.0692753, 0.0112199, 4.17],
        [0.00023035,0.00003005, 0.0000002, 0.0100000, 0.0629051, 0.0000000, 4.17],

        [0.00016361,0.00000559,-0.0065196, 0.0009991, 0.0651041, 0.0027096, 4.17],
        [0.00006400,0.00001313, 0.0000600, 0.0099997, 0.0705997, 0.0109313, 4.17],
        [0.00023132,0.00003002, 0.0000004, 0.0100000, 0.0651041, 0.0000000, 4.17],
        
        [0.00019104,0.00001031,-0.0066895, 0.0009987, 0.0569646, 0.0051531, 4.17],
        [0.00006823,0.00001808,-0.0006987, 0.0099972, 0.0574557, 0.0078683, 4.17],
        [0.00021568,0.00002653,-0.0001183, 0.0099997, 0.0569646, 0.0000000, 4.17]])


    directory = os.getcwd() + '/tests/'
    a = ut.read_fitacs_par(directory+'test_urr.par',num_e_regions=3)
    b = a[['stren','dstren','dist_R','ddist_R','gam_g','dgam_g','D']].to_numpy()

    np.testing.assert_array_almost_equal(b, test)
