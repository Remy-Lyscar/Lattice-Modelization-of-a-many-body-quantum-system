import numpy as np 
from scipy.sparse import csc_matrix 
from scipy.sparse.linalg import eigsh # Find a few eigenvalues for Hermitian matrices


def diagonalize_py(H, k=1, sigma=0, which='LM'):
    """
    Diagonalize a Hermitian matrix H and return the k lowest eigenvalues and eigenvectors.
    
    Parameters
    ----------
    H : csc_matrix
        The Hamiltonian matrix to diagonalize.
    k : int, optional
        The number of eigenvalues and eigenvectors to return. Default is 1.
    sigma : float, optional
        Find eigenvalues near sigma. Default is 0.
    which : str, optional
        Which k eigenvectors and eigenvalues to find. Default is 'LM'.
        
    Returns
    -------
    evals : np.ndarray
        The k lowest eigenvalues.
    evecs : np.ndarray
        The k corresponding eigenvectors.
    """

    evals, evecs = eigsh(H, k=k, sigma=sigma, which=which)
    return evals



