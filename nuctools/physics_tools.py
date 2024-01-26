import numpy as np
import pandas as pd

# -----------------------------------------------------------------------------
# Define some constants that can be used anywhere
# -----------------------------------------------------------------------------
hbar = 6.582119514e-16 # [eV-s] NIST Planck's over 2 pi
mnc2 = 939565413.3  # NIST neutron mass [eV]
c = 299.792458      # NIST speed of light [m/us]
Na = 6.02214076e23  # NIST Avogadro's num. [1/mol]
# -----------------------------------------------------------------------------

__all__ = ['pole_strength',
           'n_strength_func_frohner','n_strength_func','spin_stat_factor',
           'r_inf','r_prime','calc_ac','penetrability','rho','wave_number',
           'wigner_spacing']

def pole_strength(S,A,a_c):
    """
    Pole strength as defined by Frohner
    
    Parameters
    ----------
    S : float
        Neutron strength function (can be calculated with 
        physics_tools). Typical value is 1e-4. [unitless]
    A : float
        The mass number of the nucleus.
    a_c : float,optional
        The nuclear scattering radius of the isotope [fm].
    
    Returns
    -------
    s : float
        The pole strength as defined by Frohner

    Notes
    -----
    .. math:: s_c = \\frac{S_l}{2k_ca_c}\\sqrt{\\frac{E}{1 eV}} 

    .. math:: s_c = \\frac{S_l}{2\\left(\\frac{\\sqrt{2m_nE}}{\\hslash/2\\pi}\\frac{A}{A+1}\\right)a_c}\\sqrt{\\frac{E}{1 eV}} 

    .. math:: s_c = \\frac{S_l(h/(2\\pi))(A+1)}{a_cA\\sqrt{8m_n}}
    """
    
    s = S*hbar*(A+1)/a_c/1e-15/A/np.sqrt(8*mnc2/(c*1e6)**2)
    return s

def n_strength_func_frohner(e_vec,gn_vec,J,L,I,A=None,a_c=None,e_high=None):
    """
    This should really be the same as Mughabghab. This does not need 
    to be used.

    """

    # Is A defined?
    if A is not None:
        if a_c is None:
            a_c = calc_ac(A)
    else:
        raise ValueError('You have not specified A')
    # Vectors match?
    if len(e_vec) != len(gn_vec):
        raise ValueError('e_vec and gn_vec must be same length')

    if e_high is not None:
        vals = pd.DataFrame({
            'e' : e_vec,
            'gn': gn_vec,
            })
        e_vec = np.array(vals.e[vals.e < e_high])
        gn_vec = np.array(vals.gn[vals.e < e_high])
    else:
        e_vec = np.array(e_vec)
        gn_vec = np.array(gn_vec)

    # channel penetrability
    P_c = penetrability(L,a_c,A,e_vec)

    # average reduced neutron width amplitude squared, 1/1000 for meV->eV
    avg_red_width = np.mean(gn_vec/1000/2/P_c)

    # average level spacing
    D = np.mean(np.diff(e_vec))

    # strength function
    g = spin_stat_factor(J,I)
    k = wave_number(A,e_vec)
    S_l = np.mean(2*g*k*a_c/np.sqrt(e_vec)*avg_red_width)/D/(2*L+1)

    return S_l






def n_strength_func(e_vec,gn_vec,J,L,I,A=None,a_c=None,e_high=None):
    """
    Calculate the neutron strength function

    Calculate the neutron strength function as defined
    in Mughabghab's Atlas with the reduced neutron 
    amplitudes.

    Parameters
    ----------
    e_vec : array-like
        A 1-d vector of energy eigenvalues associated with
        each of the resonance amplitudes [eV]
    gn_vec : array-like
        A 1-d vector of the neutron widths [meV]
    J : float or array-like
        J is the compound nucleus spin, and takes on values
        within :math:`\\mid I-L+1/2\\mid \\leq J\\leq\\mid I+L+1/2\\mid`
        where L is the quantum angular momentum.
    L : float
        The quantum angular momentum. Can have only integer values.
    I : float
        I is the spin of the target nuclei and can be integer or half
        integer
    A : float
        The mass number of the nucleus. If a_c is not specified
        the mass number must be specified to calculate a_c
    a_c : float,optional
        The nuclear scattering radius of the isotope [fm]. a_c must
        be specified if A is not specified.
    e_high : float
        The upper energy in which to sum of the resonance values. Sum
        will go from lowest to max(e_vec) or to e_high

    Returns
    -------
    S_l : float
        The neutron strength function [dimensionless]

    Examples
    --------
    >>> import nuctools as nuc
    >>> E  = np.array([5,10,15,20])
    >>> gn = np.array([3,12, 8,21])
    >>> J,L,I,A = 4.0,0,3.5,181
    >>> print(nuc.n_strength_func(E,gn,J,L,I,A=A))

    Notes
    -----
    .. math:: S_l = \\frac{1}{(2l+1)\\Delta E}\\sum_j g_j\\Gamma^l_{n,j}

    where,

    .. math:: \\Gamma^l_{n,j} = \\sqrt{\\frac{1 eV}{E_0}}\\frac{\\Gamma_{n,j}}{V_l}

    where,

    .. math:: V_l = \\frac{P_l}{k_ca_c}

    :math:`k_c` is the wave number and :math:`a_c` is the channel radius.

    Mughabghab, S.F., Atlas of Neutron Resonances, 2006, 5th
    Edition

    """
    # Is A defined?
    if (A is not None) and (a_c is None):
        a_c = calc_ac(A)
    # Is a_c defined?
    elif( a_c is None ):
        raise ValueError('You have not specified a_c or A')
    if (A is None):
        raise ValueError('You have not specified A')
    # Vectors match?
    if len(e_vec) != len(gn_vec):
        raise ValueError('e_vec and gn_vec must be same length')

    if e_high is not None:
        vals = pd.DataFrame({
            'e' : e_vec,
            'gn': gn_vec,
            'J' : J,
            })
        e_vec = np.array(vals.e[vals.e < e_high])
        gn_vec = np.array(vals.gn[vals.e < e_high])
        J = np.array(vals.J[vals.e < e_high])
    else:
        e_vec = np.array(e_vec)
        gn_vec = np.array(gn_vec)
        J = np.array(J)

    #Mughabghab's V_l value
    V_l = penetrability(L,a_c,A,e_vec)/rho(a_c,A,e_vec)
    # Reduced neutron width amplitudes (with extra factors), 1/1000 for meV->eV
    reduced_gn = (gn_vec/1000/e_vec**(1/2.)*
        spin_stat_factor(J,I)/V_l)
    
    S_l = 1/((2*L+1)*(max(e_vec)-min(e_vec)))*reduced_gn.sum()

    N = len(e_vec)

    dS_l = (2/N)**0.5*S_l

    return S_l,dS_l

def spin_stat_factor(J,I):
    """
    The spin statistical factor 
    
    The spin statistical factor is the probability that
    a particular compound state will form for a given J
    and I when neutrons interact with the nucleus.

    Parameters
    ----------
    J : float or array-like
        J is the compound nucleus spin, and takes on values
        within :math:`\\mid I-L+1/2\\mid \\leq J\\leq\\mid I+L+1/2\\mid`
        where L is the quantum angular momentum.
    I : float
        I is the spin of the target nuclei and can be integer or half
        integer

    Returns
    -------
    g : float
        The probability of forming a state with spin compound nucleus
        spin J

    Notes
    -----
    .. math:: g = \\frac{2J+1}{2(2I+1)}
    
    """
    
    g = (2*J+1)/(2*(2*I+1))
    return g

def r_inf(Rprime,a_c=0.0,A=None):
    """
    Calculate the distant level contribution to 
    nuclear scattering radius.

    Calculate :math:`R^{\\infty}` from the neutron channel
    radius (or mass number A) and the effective nuclear 
    radius. This is for s-waves, using the R-matrix 
    formalism. **It is recommended that the user input
    parameter A, and let the program calculate :math:`a_c`.**

    Parameters
    ----------
    Rprime : float
        The effective potential scattering radius [fm] 
        including the distant levels contributions from 
        :math:`R^{\\infty}`.
    a_c : float,optional
        The nuclear scattering radius of the isotope [fm]
    A : float, optional
        The mass number of the nucleus. If a_c is not specified
        the mass number must be specified to calculate a_c

    Returns
    -------
    Rinf : float
        The distant levels contribution to the scattering
        radius [fm]

    Notes
    -----
    .. math:: R' = a_c(1-R^\\infty)

    .. math:: R^\\infty = 1-\\frac{R'}{a_c}

    Mughabghab, S.F., Atlas of Neutron Resonances, 2006, 5th
    Edition

    Examples
    --------
    >>> import nuctools as nuc
    >>> A = 181
    >>> Rprime = 7.8
    >>> Rinf = nuc.r_inf(Rprime,A=181)
    
    """

    if a_c == 0.0:
        if A is None:
            raise ValueError("If a_c is to be default, specify mass number A")
        a_c = calc_ac(A)

    Rinf = 1-Rprime/a_c

    return Rinf

def r_prime(Rinf,a_c=0.0,A=None):
    """
    Calculate the effective nuclear radius.

    Calculate the effective nuclear radius from
    the input scattering radius and the contributions
    from distant levels: :math:`R^{\\infty}`. This is 
    for s-waves, using the R-matrix formalism. **It is
    recommended that the user input parameter A, and 
    let the program calculate :math:`a_c`.**

    Parameters
    ----------
    Rinf : float
        The :math:`R^{\\infty}` term which is the contributions
        of many distant levels to the potential scattering [fm].
    a_c : float,optional
        The nuclear scattering radius of the isotope [fm]
    A : float, optional
        The mass number of the nucleus. If a_c is not specified
        the mass number must be specified to calculate a_c

    Returns
    -------
    R_prime : float
        The effective nuclear scattering radius [fm]

    Notes
    -----
    .. math:: R' = a_c(1-R^\\infty)

    Mughabghab, S.F., Atlas of Neutron Resonances, 2006, 5th
    Edition

    Examples
    --------
    >>> import nuctools as nuc
    >>> A = 181
    >>> Rinf = 0.1
    >>> R_prime = nuc.r_prime(Rinf,A=A)
    
    """

    if a_c == 0.0:
        if A is None:
            raise ValueError("If a_c is to be default, specify mass number A")
        a_c = calc_ac(A)

    R_prime = a_c*(1-Rinf)

    return R_prime

def calc_ac(A):
    """
    Calculate good guess of neutron channel radius.

    Parameters
    ----------
    A : int
        The mass number of the nucleus e.g. Ta has 
        one isotope of mass number 181 (essentially)

    Returns
    -------
    a_c : float
        The channel radius [fm] or [Fermi]

    Notes
    -----
    .. math:: a_c = 1.23A^{1/3} + 0.8

    F.H. Frohner, JEFF Report 18, p. 52

    """

    a_c = 1.23*A**(1/3.) + 0.8
    return a_c

def penetrability(L,a_c,A,E):
    """
    Hard sphere penetrability for neutron.

    Parameters
    ----------
    L : int
        The quantum angular momentum
    a_c : float
        The channel radius, which is the value seperating
        the internal and external wave functions in R-matrix,
        or the nuclear radius [fm] or [Fermi].
    A : float
        The mass number of the target isotope, e.g. neutron
        incident on Ta-181 has A=181.
    E : float
        The kinetic energy of the neutron [eV]

    Returns
    -------
    P_c : float
        The hard sphere penetrability factor

    Examples
    --------
    >>> import nuctools as nuc
    >>> L = 0
    >>> a_c = 7.8  # fm
    >>> A   = 181  # mass number
    >>> E   = 100  # eV
    >>> P_c = nuc.penetrability(L,a_c,A,E)

    """

    # calculate rho
    p = rho(a_c,A,E)


    if (L==0):
        P_c = p
    if L==1:
        P_c = p**3/(1+p**2)
    if L==2:
        P_c = p**5/(9+3*p**2+p**4)
    if L==3:
        P_c = p**7/(225+45*p**2+6*p**4+p**6)
    if L==4:
        P_c = p**9/(11025+1575*p**2+135*p**4+10*p**6+p**8)

    return P_c


def rho(a_c,A,E):
    """
    The dimensionless quantity rho

    Parameters
    ----------
    a_c : float
        The channel radius, which is the value seperating
        the internal and external wave functions in R-matrix,
        or the nuclear radius [fm] or [Fermi].
    A : float
        The mass number of the target isotope, e.g. neutron
        incident on Ta-181 has A=181.
    E : float
        The kinetic energy of the neutron [eV]

    Returns
    -------
    p : float
        The dimensionless quantity rho

    Examples
    --------
    >>> import nuctools as nuc
    >>> p = nuc.rho(7.8,181,100)
    >>> print(p)

    Notes
    -----
    .. math:: \\rho = k_ca_c

    where :math:`k_c` is the wave number and :math:`a_c` is 
    the channel radius
    """

    p = wave_number(A,E)*a_c

    return p

def wave_number(A,E):
    """
    Calculate the neutron wave number.

    Calculate the neutron wave number with the reduced
    mass of the neutron and it's target. Useful for 
    equations involving neutron interactions.

    Parameters
    ----------
    A : float
        The mass number of the target isotope, e.g. neutron
        incident on Ta-181 has A=181.
    E : float
        The kinetic energy of the neutron [eV]

    Returns
    -------
    k : float
        The neutron reduced mass wave number [1/Fermi]

    Examples
    --------

    >>> import nuctools as nuc
    >>> k = nuc.wave_number(181,100)
    >>> print(k)

    Notes
    -----
    .. math:: k_c = \\frac{\\sqrt{2m_nE}}{\\hslash/(2\\pi)}\\frac{A}{A+1}

    where :math:`h` is Planck's constant and :math:`m_n` is neutron
    mass, and A is the mass number.
    
    """

    k = np.sqrt(2*mnc2*E)/(c*1e6)/hbar * A/(A+1) * 1e-15 # 1e-15 [m/Fermi]
    return k

def wigner_spacing(x,avD):
    """
    Calculate the Wigner level-spacing distribution from the average 
    spacing D (`avD`)

    Parameters
    ----------
    x : array-like
        The abscissae for the distribution
    avD : float
        The average level spacing

    Returns
    -------
    wig : array-like
        The y-values of the probability distribution

    Examples
    --------
    >>> import nuctools as nuc
    >>> import numpy as np
    >>> x = np.linspace(0,100,100)
    >>> Px = nuc.wigner_spacing(x,4.0)


    Notes
    -----
    .. math:: P(S) = (\\pi S/2D^2)e^{-\\pi/4 S^2/D^2} dS
    """

    wig = np.pi/2/avD**2*x*np.exp(-1/4*np.pi*x**2/avD**2)
    return wig





