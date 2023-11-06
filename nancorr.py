def nancorr(A,B):
    """
    NANCORR - Pearson correlation coefficient
    
    [coef, t, n] = NANCORR(A, B)
    INPUT:
      A, B - input matrices, single or double, with equal number of rows
    
    OUTPUT:
      coef - matrix of Pearson correlation coefficients
      t    - matrix of t-statistics
      n    - matrix containing the number of defined values
    
    NOTES  
    pvalue can be calculated as 2*tcdf(-abs(t), n - 2)
    """
    
    import numpy as np

    Am=~np.isfinite(A); Bm=~np.isfinite(B)

    if issubclass(A.dtype.type,np.float32): # np.floating
        Ap=np.single(~Am); Bp=np.single(~Bm); # float32    
    else:
        Ap=np.double(~Am); Bp=np.double(~Bm); # float64

    # zero out nan elements
    A[Am]=0; B[Bm]=0

    #   code one of the formulas from https://en.wikipedia.org/wiki/Pearson_correlation_coefficient
    #   this procedure might be numericaly unstable for large values,
    #   it might be reasonable to center each column before calling nancorr.

    xy = A.T @ B         # sum x_i y_i
    n  = Ap.T @ Bp        # number of items defined both in x and y
    
    mx = A.T @ Bp / n    # mean values in x, calculated across items defined both in x and y
    my = Ap.T @ B / n    # mean values in y, calculated across items defined both in x and y
    
    x2 = (A*A).T @ Bp    # sum x^2_i, calculated across items defined both in x and y
    y2 = Ap.T @ (B*B)    # sum y^2_i, calculated across items defined both in x and y
    
    sx   = np.sqrt(x2 - n * (mx**2));  # sx, sy - standard deviations 
    sy   = np.sqrt(y2 - n * (my**2));
    
    coef = (xy - n * mx * my) / (sx * sy);      # correlation coefficient
    t    = coef * np.sqrt((n - 2) / (1 - coef**2));  # t-test statistic

    return coef, t, n
