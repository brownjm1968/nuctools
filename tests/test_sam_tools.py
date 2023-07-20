import nuctools.sam_tools as st
import numpy as np
from numpy.testing import assert_array_equal,assert_allclose
import os


def test_read_pds():
    pds_file = os.getcwd() + '/tests/test_files/sammy.pds'

    df = st.read_pds(pds_file)

    gold = np.array([0.794262,4.302552E-02,10.3511,-0.200283,3.46130,-3.926995E-02,-6.50670, 38.5814,5.76615]) 

    assert_allclose(gold,df.iloc[0])