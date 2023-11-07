import numpy as np
import gwas_tools
A = np.random.rand(10, 3)
B = np.random.rand(10, 2)

A[1, 1] = np.nan

r,t,n = gwas_tools.nancorr(A, B)

for i in range(A.shape[1]):
    for j in range(B.shape[1]):
        idx = np.isfinite(A[:, i] + B[:, j])
        assert np.isclose(r[i, j], (np.corrcoef(A[idx, i], B[idx, j])[0, 1]))
print('ok')
