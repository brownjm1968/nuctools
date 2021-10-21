import numpy as np
from nuctools import *


def test_sys_cov():
    a = np.array([2,2])
    b = np.array([4,4])
    c = 1.0
    da = np.array([2,2])
    db = np.array([2,2])
    dc = 2.0

    f = a*c/b

    dfdc = a/b

    answer = np.ones((2,2))

    cy = sys_cov([dc],[dfdc])

    np.testing.assert_array_equal(answer,cy)

def test_stat_cov():
    a = np.array([2,2])
    b = np.array([4,4])
    da = np.array([2,2])
    db = np.array([2,2])

    f = a/b

    dfda = 1/b
    dfdb = a/b**2

    answer = np.array([[5/16,0],[0,5/16]])

    cy = stat_cov([da,db],[dfda,dfdb])

    np.testing.assert_array_equal(answer,cy)


def test_symmetric2dTest():
    true = [[1,2,3],
            [2,1,2],
            [3,2,1]]
    false = [[1,2,3],
             [1,2,3],
             [1,2,3]]

    assert symmetric2dTest(true) == True
    assert symmetric2dTest(false) == False

def test_covar():

    fx = np.array([[1,2],[1,2]])
    cx = np.array([[2,0],[0,2]])
    
    cy = covar(fx,cx)
    
    answer = np.array([[10,10],[10,10]])

    np.testing.assert_array_equal(cy,answer)


def test_cov_to_corr():

    fx = np.array([[1,2],[1,2]])
    cx = np.array([[2,0],[0,2]])
    
    cy = covar(fx,cx)

    corr = cov_to_corr(cy)

    answer = np.array([[1.0,1.0],[1.0,1.0]])

    np.testing.assert_array_equal(corr,answer)


def test_prop_err_mean():

    x = [1,2,3,4]
    dx = [0.2,0.2,0.2,0.2]

    d_avg_x = prop_err_mean(dx)

    answer = 0.1

    assert d_avg_x == answer



