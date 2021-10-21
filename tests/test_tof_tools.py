import numpy as np
import pandas as pd
import nuctools.tof_tools as tt

mnc2 = 939565413.3  # NIST neutron mass [eV]
c = 299.792458      # NIST speed of light [m/us]

def test_tofe():
	tof = 100.0
	FP = 100.0

	E = mnc2*(1/np.sqrt(1-(FP/tof/c)**2)-1)
	assert tt.tofe(tof,FP) == E

def test_etof():
	E = 100.0
	FP = 100.0

	tof = FP/c*1/np.sqrt(1-1/(E/mnc2+1)**2)
	assert tt.etof(E,FP) == tof

def test_single_group():
	# x-axis
	x = np.arange(8)
	# random y
	answer = np.arange(len(x))*np.random.random()
	
	# group size
	gs = 2
	# data to be grouped
	xdata = np.repeat(x,gs)
	data = np.repeat(answer,gs)
	
	# grouped data
	gdata = tt.single_group([xdata,data],gs)

	
	np.testing.assert_allclose(gdata[2]/gs,answer)