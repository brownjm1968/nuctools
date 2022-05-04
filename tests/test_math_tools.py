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

def test_flat_triang_to_full():

    gold1 = np.array([[1,3,5],
                      [3,4,9],
                      [5,9,2]])
    gold2 = np.array([[1,3,5,7],
                      [3,4,9,6],
                      [5,9,2,3],
                      [7,6,3,4]])
    gold3 = np.array([[ 1, 2, 3, 4, 5, 6, 7, 8, 9,10],
                      [ 2,11,12,13,14,15,16,17,18,19],
                      [ 3,12,20,21,22,23,24,25,26,27],
                      [ 4,13,21,28,29,30,31,32,33,34],
                      [ 5,14,22,29,35,36,37,38,39,40],
                      [ 6,15,23,30,36,41,42,43,44,45],
                      [ 7,16,24,31,37,42,46,47,48,49],
                      [ 8,17,25,32,38,43,47,50,51,52],
                      [ 9,18,26,33,39,44,48,51,53,54],
                      [10,19,27,34,40,45,49,52,54,55]])

    lt1 = np.array([1,3,5,4,9,2])
    lt2 = np.array([1,3,5,7,4,9,6,2,3,4])
    lt3 = np.arange(1,56)

    np.testing.assert_array_equal( flat_triang_to_full(lt1), gold1 )
    np.testing.assert_array_equal( flat_triang_to_full(lt2), gold2 )
    np.testing.assert_array_equal( flat_triang_to_full(lt3), gold3 )

