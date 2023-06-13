from nuctools import *
import numpy as np
import h5py as h5
from numpy.testing import assert_array_equal,assert_allclose
import os,sys


def test_inc():
	assert inc(4) == 5

def test_get_montrig():
	sample = "test"
	i = np.zeros(8)
	j = np.zeros(2)
	with h5.File("tests/{}_master.h5".format(sample),"w") as hdf:
	    temp_mon = np.array([100.0,200.0,300.0,0.0,0.0,0.0,0.0,400.0])
	    temp_trig = np.array([[3e8,3.1e8],[0,0]])
	    temp_data = np.array([1,2,3,4,5])
	    mon = hdf.create_dataset("test/CYCLE0001/MON",data=temp_mon)
	    trig = hdf.create_dataset("test/CYCLE0001/TRIG",data=temp_trig)
	    data_g = hdf.create_dataset("test/CYCLE0001/DATA/GOOD",data=temp_data)
	    data_t = hdf.create_dataset("test/CYCLE0001/DATA/TRUNC",data=temp_data)
	    mon2 = hdf.create_dataset("test/CYCLE0002/MON",data=temp_mon)
	    trig2 = hdf.create_dataset("test/CYCLE0002/TRIG",data=temp_trig)
	    data_g2 = hdf.create_dataset("test/CYCLE0002/DATA/GOOD",data=temp_data)
	    data_t2 = hdf.create_dataset("test/CYCLE0002/DATA/TRUNC",data=temp_data)
    
	with h5.File("tests/{}_master.h5".format(sample),"r") as hdf:
	    for cycle in hdf[sample+"/"]:
	        print(cycle)
	        mon_vector = hdf[sample+"/"+cycle+"/MON"]
	        trig_vector = hdf[sample+"/"+cycle+"/TRIG"][0]
	        datag_array = hdf[sample+"/"+cycle+"/DATA/GOOD"]
	        datat_array = hdf[sample+"/"+cycle+"/DATA/TRUNC"]
	        sum_detcts = len(datag_array)+len(datat_array)
	        i+=mon_vector
	        j+=trig_vector

	manual_montrig = np.concatenate((i,j))
	
	fnc_montrig = get_montrig("tests/",sample)

	assert_array_equal(manual_montrig,fnc_montrig)

def test_count_rate():
	tof_arr = np.array([[1,2,3],[20,10,5]])
	trig = 3e8
	tau = .6
	bw = 1.0

	# convert to cps
	cps = tof_arr[1]/bw/1e-6/trig
	# error cps
	e_cps = np.sqrt((np.sqrt(tof_arr[1])/bw/trig/1e-6)**2)
	# dead-time correction
	dead_time_corr = dtc(tof_arr[1],tau,trig,bw)
	# dead-time corrected cps
	dc_cps = cps*dead_time_corr[0]
	# error in dc cps
	e_dc_cps = np.sqrt((e_cps/cps)**2+(dead_time_corr[1]/dead_time_corr[0])**2)*dc_cps
	a = np.array([tof_arr[0],dc_cps,e_dc_cps])

	fcn_cps = count_rate(tof_arr,trig,tau) 

	assert_allclose(a,fcn_cps)

def test_dtc():
	counts = np.linspace(10,1,10)
	dt = 1.125
	w = 0.5
	trigs = 100

	lc = len(counts)
	test_dtcf = np.zeros(lc)
	test_ddtcf = np.zeros(lc)
	# SUM from Dr. Danon's dead_time.pdf (maybe his thesis...)
	SUM = 0.75*counts[0:lc-2] + counts[1:lc-1] + 0.5*counts[2:lc]
	# the error on the sum
	dSUM = np.sqrt(0.75**2*counts[0:lc-2]+counts[1:lc-1]+0.5**2*counts[2:lc])

	test_dtcf[2:lc] = 1/(1-SUM/trigs)
	test_ddtcf[2:lc] = dSUM/trigs/(1-SUM/trigs)**2

	dtcf,ddtcf = dtc(counts,dt,trigs,w)

	assert_allclose(test_dtcf,dtcf)
	assert_allclose(test_ddtcf,ddtcf)

