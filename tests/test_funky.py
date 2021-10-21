import numpy as np
import pandas as pd
import nuctools.funky as f

mnc2 = 939565413.3  # NIST neutron mass [eV]
c = 299.792458      # NIST speed of light [m/us]




def test_threshFinder():
	pulse = [0,1,2,4,6,8,10,9,8,7,6,5,4,3,2,1,0]
	threshold = 2.0
	answer = 3 #---^

	assert f.threshFinder(pulse,threshold) == answer

def test_decayFinder():
	pulse = [0,1,2,4,6,8,10,9,8,7,6,5,4,3,2,1,0]
	threshold = 2.0
	answer = 15 #---------------------------^

	assert f.decayFinder(pulse,threshold) == answer

def test_giveIt():
	assert f.giveIt() == 900

