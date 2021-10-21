from nuctools import *
import numpy as np
from numpy.testing import assert_array_equal,assert_allclose
import os
import filecmp

def test_read_3col_pendf():
	answer = np.array([[1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8],
		               [1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5]])

	file = os.getcwd() + '/tests/test_files/test_pendf.dat'
	start = 98
	finish = 100
	ptwise_cross_sec = read_3col_pendf(file,start,finish)

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

