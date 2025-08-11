from nuctools import *
import numpy as np
from numpy.testing import assert_array_equal,assert_allclose
import os
import filecmp

def test_TAB1():

    file = os.getcwd() + '/tests/test_files/test_pendf.dat'
    mat,mf,mt = 7328,3,1

    ta = TAB1(file,mat,mf,mt)

    gold_array = np.array([[1. , 1.5],
                           [1.1, 1.5],
                           [1.2, 1.5],
                           [1.3, 1.5],
                           [1.4, 1.5],
                           [1.5, 1.5],
                           [1.6, 1.5],
                           [1.7, 1.5],
                           [1.8, 1.5]])

    assert_array_equal(ta.data.to_numpy(),gold_array)

def test_read_endf_float():

    gold_float = 7e4

    assert read_endf_float(" 7.000000+4") == gold_float

def test_write_endf_float():

    # cases represent some common failures, like the integer 10088013 which cannot have 
    # exponential format or it loses precision
    gold_string_list = ["-7.000000+4",
                        " 1.024200-5",
                        " 4.217507-2",
                        " 0.42175079",
                        " 1.234560+9",
                        " 1.23400-14",
                        " 1.1000-100",
                        "-9.990090-0",
                        "-2.438383+8",
                        "-2.438383+2",
                        "-2.440000+7",
                        " 10088013.0"]
    inp_value_list   = [-7e4,
                        0.000010241999999999999517471920007505303829020704142749309539794921875,
                        4.217507e-2,
                        0.421750794,
                        1.234560e+9,
                        1.234e-14,
                        1.1e-100,
                        -9.990090-0,
                        -2.438383e8,
                        -2.438383e2,
                        -2.44000e7,
                        10088013]

    for i,val in enumerate(inp_value_list):
        assert write_endf_float(val) == gold_string_list[i]

def test_read_pendf_xs():
    answer = np.array([[1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8],
                       [1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5]])

    file = os.getcwd() + '/tests/test_files/test_pendf.dat'
    start = 98
    finish = 100
    ptwise_cross_sec = read_pendf_xs(file,start,finish)

    assert_array_equal(answer,ptwise_cross_sec)


def test_write_pendf_xs():
    energy = [1e2,2e2,3e3,4e3,5e3,6e3,7e3,8e3,9e3,1e4]
    cs = [2,2,2,2,2,2,2,2,2,2]
    material_number = 7328
    file_number = 3
    reaction_number = 1
    test_file = os.getcwd() + '/tests/test_files/temp_pendf_xs.dat'
    ans_file = os.getcwd() + '/tests/test_files/test_pendf_write.dat'

    write_pendf_xs(test_file,energy,cs,material_number,file_number,reaction_number)

    assert filecmp.cmp(test_file,ans_file) == True

