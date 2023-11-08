from nuctools import *
import numpy as np
import pandas as pd
from matplotlib.pyplot import *
from numpy.testing import assert_array_equal,assert_allclose
import os


def test_sum_tof():
    directory = os.getcwd() + '/tests/test_files/'
    pa_answer = pd.read_csv(directory+'pa_la.csv').to_numpy()
    tof_answer = pd.read_csv(directory+'tof_la.csv')[['uwadc1','uwadc2']].to_numpy()

    filelist = [directory+'la_fp4a_t03_r03_0022.lst']
    tof_df,pa_df = sum_tof(directory+'la-agl.json',filelist)
    tof_uw_df = tof_df[['uwadc1','uwadc2']]
    
    assert_allclose(pa_answer,pa_df.to_numpy())
    assert_allclose(tof_answer,tof_uw_df.to_numpy())

