import numpy as np



__all__ = ['stat_cov','sys_cov','build_sys_cx',
           'symmetric2dTest','covar','cov_to_corr','prop_err_mean',
           'flat_triang_to_full']


def stat_cov(stat_err,stat_der):
    """
    Statistical covariance from the given errors and derviatives,
    which will always be diagonal as they are uncorrelated. This is
    meant to be used along with sys_cov() to get the full covariance
    matrix.

    Parameters
    ----------
    stat_err : array-like
        The statistical error in the variables. This should be 2-d,
        where each position stat_err[j] is a vector, usually TOF 
        or energy.
    stat_der : array-like
        The derivative of the variable in whichever equation you 
        are building the covariance for. This is 2-d, where each
        position stat_der[j] is a vector with length equal to the
        length of the function you are building covariance for.
        **The order of stat_der should be same as stat_err**.

    Returns
    -------
    cy : numpy array
        The covariance given the input error and derivatives. This
        is an NxN array where N = len(derivative vector)

    Notes
    -----
    This will return a diagonal matrix with each element ii as:

    .. math:: C_{y,ii} = \\sum_j \\left(\\frac{\\partial f}{\\partial j}\\Delta j\\right)^2

    where the derivatives of f are taken with respect to each of the consitituents
    j of the f equation. i goes from 1..N where N = len(f).

    Examples
    --------
    >>> a = np.array([1,2,3,4])
    >>> b = np.array([4,3,2,1])
    >>> da = np.array([0.1,0.2,0.1,0.2])
    >>> db = np.array([0.1,0.2,0.1,0.2])
    >>> f = a/b
    >>> dfda = 1/b
    >>> dfdb = a/b**2
    >>> cy = nuc.stat_cov([da,db],[dfda,dfdb])


    """
    # cast as numpy arrays
    stat_err = np.array(stat_err)
    stat_der = np.array(stat_der)
    # square the error for the variance
    stat_err2 = np.power(stat_err,2)
    # square the derivatives
    stat_der2 = np.power(stat_der,2)
    # sum of (derivative * error)^2
    cy = np.diag((stat_err2*stat_der2).sum(axis=0))
    
    # return the diagonal matrix
    return cy


def sys_cov(sys_err,sys_der,cov_mat=None):
    """
    Systematic covariance from the given errors and derviatives.

    Parameters
    ----------
    sys_err : array-like
        The error associated with each of the systematic variables. We 
        assume that these parameters are not already correlated. The 
        sys_err array is expected to be 1-d, with len(sys_err) = the 
        number of systematic variables are in the equation. As these
        variable do not change over the function, they are not vectors.
    sys_der : array-like
        The derivatives of the function of interest with repsect to 
        each of the variables. Each sys_der[i] should be the length
        of the function of interest. **The order of sys_der should be 
        same as sys_err, with any cov_mat variables plugged into the
        leading rows **.
    cov_mat : array-like
        This is a covariance matrix that can be plugged into the 
        leading rows and columns. This allows already correlated data
        to be used in the systematic error. (sys_err still has to 
        contain variables for now.)

    Returns
    -------
    cy : numpy array
        The covariance given the input error and derivatives.

    Notes
    -----
    The output covariance will be NxN where N = len(function) of interest.
    The algebra to get the covariance is shown below.

    .. math:: \\mathbf{C_y} = \\mathbf{F_x}\\mathbf{C_x}\\mathbf{F_x}^T

    Here the dimensions of F_x are NxM where M is the number of systematic 
    variables (and the len(sys_der)), and N is the len(function) (and the 
    len(sys_der[0])). The matrix C_x is MxM.

    Examples
    --------
    >>> a = np.array([1,2,3,4])
    >>> b = np.array([4,3,2,1])
    >>> c = 3.0
    >>> da = np.array([0.1,0.2,0.1,0.2])
    >>> db = np.array([0.1,0.2,0.1,0.2])
    >>> dc = 0.1
    >>> f = a*(b*c)**-1
    >>> dfdc = a*(b*c)**-2*(-1)*(b)
    >>> cy = sys_cov([dc],[dfdc])

    """

    if cov_mat is not None:
        lcm = len(cov_mat)
        N = lcm + len(sys_err)
        cx = np.zeros((N,N))
        # plug in matrix at the leading rows and cols
        cx[0:lcm,0:lcm] = cov_mat
        
    # cast as numpy arrays
    sys_err = np.array(sys_err)
    sys_der = np.array(sys_der)
    # square the error for the variance and diagonalize
    sys_err2 = np.diag(np.power(sys_err,2))
    # determine the final input covariance matrix
    if cov_mat is None:
        cx = sys_err2
    else:
        cx[lcm:,lcm:] = sys_err2
    # sandwich rule, matrix algebra to get outgoing covariance
    # (@ symbol is matrix dot product)
    cy = sys_der.T@cx@sys_der
    
    # return NxN matrix (N = len(derivatives))
    return cy

def build_sys_cx(sys_err,cov_mat=None):
    """
    Systematic covariance input matrix from the given errors.

    Parameters
    ----------
    sys_err : array-like
        The error associated with each of the systematic variables. We 
        assume that these parameters are not already correlated. The 
        sys_err array is expected to be 1-d, with len(sys_err) = the 
        number of systematic variables are in the equation. As these
        variable do not change over the function, they are not vectors.
    cov_mat : array-like
        This is a covariance matrix that can be plugged into the 
        leading rows and columns. This allows already correlated data
        to be used in the systematic error. (sys_err still has to 
        contain variables for now.)

    Returns
    -------
    cy : numpy array
        The covariance given the input error and derivatives.

    Notes
    -----
    Same as sys_cov(), but only builds the input matrix. Initial purpose
    is for printing the input matrix to a file.

    """
    if cov_mat is not None:
        lcm = len(cov_mat)
        N = lcm + len(sys_err)
        cx = np.zeros((N,N))
        # plug in matrix at the leading rows and cols
        cx[0:lcm,0:lcm] = cov_mat
        
    # cast as numpy arrays
    sys_err = np.array(sys_err)
    # square the error for the variance and diagonalize
    sys_err2 = np.diag(np.power(sys_err,2))
    # determine the final input covariance matrix
    if cov_mat is None:
        cx = sys_err2
    else:
        cx[lcm:,lcm:] = sys_err2

    return cx

def symmetric2dTest(matrix2d):
    """
    Test 2-dimensional matrix on whether it is symmetric

    Parameters
    ----------
    matrix2d: numpy array
        2-dimensional array to be tested for symmetry

    Returns
    -------
    symmBool: Boolean 
        True or False value whether matrix is symmetric
    """
    
    # is the matrix 2-d?
    if len(np.shape(matrix2d)) != 2:
        raise ValueError("Matrix dimensions are not equal to 2.")
    matrix2d = np.array(matrix2d)

    # create boolean for whether 2-d matrix = its transpose
    symmBool = (matrix2d == matrix2d.T).all()
            

    if symmBool == False:
        print("Matrix not symmetric.")
        print("Max assymetry = ",np.max(matrix2d-matrix2d.T))

    return symmBool




def covar(fx,cx):
    """
    Find covariance of a function whose inputs covariance matrix
    are in cx, and the Jacobian of the function with respect to
    it's input variables is in fx.

    Parameters
    ----------
    fx : numpy array
        The sensitivity matrix (or Jacobian) that has
        N rows (corresponding to the number of function
        values) and M columns (the number of input values
        to the function). NxM matrix.
    cx : numpy array
        The input covariance matrix which corresponds to 
        fx. This should include every input (a_1,b_1,a_2,
        b_2...) which is differentiated in fx. MxM matrix.


    Returns
    -------
    cy : numpy array
        The resultant covarianc matrix due to the input 
        covariance matrix cx, and sensitivities of the 
        given function in fx. NxN matrix.

    Examples
    --------
    >>> import nuctools as nuc
    >>> fx = np.array([[1,2],[1,2],[1,2]])
    >>> cx = np.array([[2,0],[0,2]])
    >>> print(fx.shape,"x",cx.shape,"x",fx.T.shape)
    >>> cy = nuc.covar(fx,cx)
    >>> print(cy.shape)

    Notes
    -----
    For independent variables at x_i:

    .. math:: a_i, b_i

    .. math:: C_y = F_xC_xF_x^T

    .. math:: f(x_i) = \\frac{a_i}{b_i}

    The covariance between two points in the function are:

    .. math:: covar(f(x_1),f(x_2)) = F_xC_xF_x^T

    .. math:: F_x = \\left[ \\begin{array}{cccc} \\frac{\\partial f(x_1)}{\\partial a_1} & \\frac{\\partial f(x_1)}{\\partial b_1} & \\frac{\\partial f(x_1)}{\\partial a_2} & \\frac{\\partial f(x_1)}{\\partial b_2} \\cr \\frac{\\partial f(x_2)}{\\partial a_1} & \\frac{\\partial f(x_2)}{\\partial b_1} & \\frac{\\partial f(x_2)}{\\partial a_2} & \\frac{\\partial f(x_2)}{\\partial b_2} \\end{array} \\right]
    
    .. math:: C_x = \\left[ \\begin{array}{cccc} \\Delta a_1^2 & 0 & 0 & 0 \\cr 0 & \\Delta b_1^2 & 0 & 0 \\cr 0 & 0 & \\Delta a_2^2 & 0 \\cr 0 & 0 & 0 & \\Delta b_2^2  \\end{array} \\right]

    """
    
    fx = np.array(fx)
    cx = np.array(cx)
    
    shape_fx = fx.shape
    shape_cx = cx.shape
    
    
    if shape_fx[1] != shape_cx[0]:
        print('-----------------------------------------')
        print("Shapes of fx and cx cannot be multiplied:")
        print(shape_fx,"x",shape_cx)
        print('-----------------------------------------')
        raise ValueError('Input matrices are not compliant')
    
    cy = np.dot(np.dot(fx,cx),fx.T)
    
    print("Size of Cy matrix: ",np.shape(cy))
    
    return cy

def cov_to_corr(cy):
    """
    Convert covariance matrix to it's correlation matrix

    Parameters
    ----------
    cy : numpy array
        The input covariance matrix. Should be square and
        symmetric.

    Returns
    -------
    corr : numpy array
        The corresponding correlation matrix

    Examples
    --------
    >>> import nuctools as nuc
    >>> fx = np.array([[1,2],[1,2],[1,2]])
    >>> cx = np.array([[2,0],[0,2]])
    >>> print(fx.shape,"x",cx.shape,"x",fx.T.shape)
    >>> cy = nuc.covar(fx,cx)
    >>> print(cy.shape)
    >>> corr = nuc.cov_to_corr(cy)

    Notes
    -----
    The relationship of correlation to covariance is:

    .. math:: \\frac{cov(X,Y)}{sd(X)sd(Y)}

    """
    
    N = len(cy)
    
    corr = np.zeros((N,N))
    
    sd = np.sqrt(np.diag(cy))
    
    sdinv = np.diag(1/sd)
    
    corr = np.dot(np.dot(sdinv,cy),sdinv)
    
    #print(np.shape(corr))
    return corr


def prop_err_mean(dx,axis=None):
    """
    Calculate the propagated error for the mean of values
    whose error is known. Only observed to work for 1-d 
    2-d arrays.


    Parameters
    ----------
    dx : 1-d array
        The error on the values which are input for the mean

    Returns
    -------
    dmean : numpy array
        The propagated error on the mean given the error on 
        the averaged values.

    Examples
    --------
    >>> import nuctools as nuc
    >>> vals = np.array([1,2,3])
    >>> vals_error = np.array([0.1,0.1,0.1])
    >>> mean = np.mean(vals)
    >>> error_of_mean = nuc.prop_err_mean(vals_error)

    Notes
    -----
    If axis is specified, the average can be used for 2-d 
    arrays in axis-0,1.

    .. math:: \\Delta\\langle x\\rangle = \\frac{1}{N}\\sqrt{\\sum_i \\Delta x_i^2}

    """
    dx = np.array(dx)
    if axis is not None:
        if axis == 0:
            N = len(dx)
        if axis == 1:
            N = len(dx[0])
    else:
        N = len(dx)
    dx2 = np.power(dx,2)
    dmean = 1/N*np.sqrt(dx2.sum(axis=axis))
    return dmean

def flat_triang_to_full(triang):
    """
    Take a flat 1-dim array that represents upper-triang
    matrix and return the full symmetric matrix

    Parameters
    ----------
    triang : array-like
        Should be 1-dimensional, upper-triangular row-major order

    Returns
    -------
    full : array-like
        The full symmetric array

    Examples
    --------
    >>> import nuctools as nuc
    >>> import numpy as np
    >>> lt = np.array([1,3,5,4,9,2])
    >>> nuc.flat_triang_to_full(lt)
    array([[1., 3., 5.],
          [3., 4., 9.],
          [5., 9., 2.]])

    Notes
    -----
    Future changes may be needed to allow options for upper/lower 
    and row-/column-major ordering of the flat array

    """
    rows = int(1/2*(np.sqrt(8*len(triang)+1)-1))
    mask = np.tri(rows,dtype=bool).T # transpose: lower to upper
    full = np.zeros((rows,rows))
    full[mask] = triang
    full += full.T*np.tri(rows,k=-1)
    return full




