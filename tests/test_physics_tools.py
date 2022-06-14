import numpy as np
import nuctools.physics_tools as pt

def test_spin_stat_factor():
    J = 2.0
    I = 3.5

    g = (2*J+1)/(2*(2*I+1))

    assert g == pt.spin_stat_factor(J,I)


def test_r_inf():
    r_prime = 9.0
    a_c = 10.0
    Rinf = 1.0 - r_prime/a_c

    r_inf = pt.r_inf(r_prime,a_c=a_c)

    assert Rinf == r_inf

def test_r_prime():
    a = 10.0
    Rinf = 0.1
    Rp_answer = 9.0

    Rp = pt.r_prime(Rinf,a_c=a)

    assert Rp == Rp_answer

def test_wave_number():
    A = 2.0
    E = 2.0

    per_Fermi = 1e-15 # how many meters per Fermi

    # sqrt{2mE}/hbar * A/(A+1) 
    k = 2*np.sqrt(pt.mnc2)/(pt.c*1e6)/pt.hbar * 2/3 * per_Fermi

    assert k == pt.wave_number(A,E)


def test_rho():
    A = 2.0
    E = 2.0
    a_c = 5.0

    per_Fermi = 1e-15 # how many meters per Fermi

    # sqrt{2mE}/hbar * A/(A+1) 
    k = 2*np.sqrt(pt.mnc2)/(pt.c*1e6)/pt.hbar * 2/3 * per_Fermi

    p = k*a_c

    assert p == pt.rho(a_c,A,E)

def test_penetrability():
    A = 2.0
    E = 2.0
    a_c = 5.0

    per_Fermi = 1e-15 # how many meters per Fermi

    # sqrt{2mE}/hbar * A/(A+1) 
    k = 2*np.sqrt(pt.mnc2)/(pt.c*1e6)/pt.hbar * 2/3 * per_Fermi

    p = k*a_c 

    assert pt.penetrability(0,a_c,A,E) == p
    
    assert pt.penetrability(1,a_c,A,E) == p**3/(1+p**2)
    
    assert pt.penetrability(2,a_c,A,E) == p**5/(9+3*p**2+p**4)
    
    assert pt.penetrability(3,a_c,A,E) == p**7/(225+45*p**2+6*p**4+p**6)
    
    assert pt.penetrability(4,a_c,A,E) == p**9/(11025+1575*p**2+135*p**4+10*p**6+p**8)

def test_wigner_spacing():
    x = np.array([2.0])  # spacing
    D = 2.0              # average spacing

    assert pt.wigner_spacing(x,D) == np.pi/4*np.exp(-np.pi/4)
