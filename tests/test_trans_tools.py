from nuctools import *
import numpy as np
import pandas as pd
from numpy.testing import assert_array_equal,assert_allclose
import os

def test_Trans():
	open_file_name = os.getcwd() + '/tests/test_files/test_Ogrp_open.dat'
	samp_file_name = os.getcwd() + '/tests/test_files/test_Ogrp_samp.dat'
	oMon_file_name = os.getcwd() + '/tests/test_files/test_Mon_open.dat'
	sMon_file_name = os.getcwd() + '/tests/test_files/test_Mon_samp.dat'

	odat = pd.read_csv(open_file_name,sep=r'\s+',
		               names=['tof','bin_width','counts'])
	sdat = pd.read_csv(samp_file_name,sep=r'\s+',
		               names=['tof','bin_width','counts'])
	omon = pd.read_csv(oMon_file_name,names=['mon'])
	smon = pd.read_csv(sMon_file_name,names=['mon'])

	# default behavior of Trans
	odat['cps'] = odat.counts/odat.bin_width/omon.mon[0]/1e-6
	sdat['cps'] = sdat.counts/sdat.bin_width/smon.mon[0]/1e-6

	odat['bkg'] = np.repeat(0.5,len(odat))/odat.bin_width/omon.mon[0]/1e-6
	sdat['bkg'] = np.repeat(0.5,len(sdat))/sdat.bin_width/smon.mon[0]/1e-6

	trans = ( sdat.cps - sdat.bkg ) / ( odat.cps - odat.bkg )

	Trans_test = Trans(open_file_name,oMon_file_name,
		               samp_file_name,sMon_file_name,
		               o_bkgcoeff=[0.0,0.0,0.5/odat.bin_width[0]/omon.mon[0]/1e-6],
		               s_bkgcoeff=[0.0,0.0,0.5/odat.bin_width[0]/omon.mon[0]/1e-6])

	assert_array_equal(Trans_test.trans,trans)

def test_calc_count_rate():
	"""
	This test may be superfluous
	"""

	counts = 2.0
	bin_width = 0.5
	triggers = 10

	# counts/bin_width/triggers/1e-6
	cps = 0.4/1e-6
	# np.sqrt(counts)/bin_width/triggers/1e-6
	dcps = np.sqrt(counts)/5/1e-6

	assert cps,dcps == Trans.calc_count_rate(counts,bin_width,triggers)






















