import nuctools.rpixdr_tools as rxt
import numpy as np
from numpy.testing import assert_array_equal,assert_allclose
import os

def test_Rpy_xdr():

    mon_sample_0 = np.array([[200000,300000],[200000,300000],
                             [200000,300000],[200000,300000],
                             [200000,300000],[200000,300000],
                             [200000,300000],[200000,300000],
                             [200000,300000],[200000,300000],
                             [22,55]]) # <- detector sums
    mon_sample_1 = np.array([[200000,300000],[200000,300000],
                             [200000,300000],[200000,300000],
                             [200000,300000],[200000,300000],
                             [200000,300000],[200000,300000],
                             [200000,300000],[200000,300000],
                             [22,55]]) # <- detector sums

    det_sum_0 = np.repeat(2,11) + np.repeat(5,11)
    det_sum_1 = np.repeat(2,11) + np.repeat(5,11)


    # define inputs to the class
    folder = os.getcwd() + '/tests/test_files/'
    number_of_samples = 2

    rpixdr_class = rxt.Rpy_xdr(folder,num_samp=number_of_samples,raw_bw=0.0064)

    assert_array_equal(rpixdr_class.monitors[0],mon_sample_0)
    assert_array_equal(rpixdr_class.monitors[1],mon_sample_1)

def test_Rpy_xdr_std():

    N = 2
    mon = np.array((200000,300000))
    det = np.array((22,55))
    mr = mon/det

    std = np.repeat(np.sqrt(np.sum( (mr-np.mean(mr))**2 )/N)/np.mean(mr), 10)
    std = np.array((std,std))


    # define inputs to the class
    folder = os.getcwd() + '/tests/test_files/'
    number_of_samples = 2

    rpixdr_class = rxt.Rpy_xdr(folder,num_samp=number_of_samples,raw_bw=0.0064)

    assert_array_equal(std,np.array(rpixdr_class.mon_std)[:,:10])



