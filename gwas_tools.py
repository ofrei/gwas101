import numpy as np
import pandas as pd
import sys

def get_byte_geno_map():
    """
    Construct mapping between bytes 0..255 and 4-element arrays of a1 genotypes
    from plink bed file.
    Return 256 x 4 array A, where A[i] = [a1, a2, a3, a4], each ai from {2, -1, 1, 0}.
    """
    genotype_codes = np.array([2, -1, 1, 0],dtype=np.int8)
    byte_map = np.empty((256,4), dtype=np.int8)
    for b in range(256):
        for a in range(4):
            byte_map[b,a] = genotype_codes[(b >> 2*a) & 3]
    return byte_map


def get_geno_byte_map():
    genotype_codes = np.array([2, -1, 1, 0],dtype=np.int8)
    geno_map = {}
    for b in range(256):
        four_geno = tuple(genotype_codes[(b >> 2*a) & 3] for a in range(4))
        geno_map[four_geno] = b
    return geno_map


def get_snp_geno(bed, i_snp, n_samples, byte_map=None):
    """
    Get {2, -1, 1, 0} genotype for i-th SNP from bed file.
    """
    if byte_map is None:
        byte_map = get_byte_geno_map()
    i_snp_geno = np.empty(bed.shape[1]*4, dtype=np.int8)
    j = 0
    for b in bed[i_snp]:
        i_snp_geno[j:j+4] = byte_map[b]
        j += 4
    return i_snp_geno[:n_samples]


def get_snp_bytes(snp_geno, geno_map=None):
    if geno_map is None:
        geno_map = get_geno_byte_map()
    n_samples = len(snp_geno)
    n_samples_mod4 = n_samples%4
    i_end = n_samples - n_samples_mod4
    bb = bytearray(geno_map[tuple(snp_geno[i:i+4])] for i in range(0,i_end,4))
    if n_samples_mod4 != 0:
        four_geno = tuple(list(snp_geno[-n_samples_mod4:]) + [0]*(4 - n_samples_mod4))
        bb.append(geno_map[four_geno])
    return bytes(bb)


def count_lines(filename):
    with open(filename, 'r', encoding='utf-8', errors='ignore') as f:
        return sum(1 for line in f)


BYTE_GENO_MAP = get_byte_geno_map()


def read_bed(bedprefix, snps=None, num_samples=None):
    snp_list = list(range(count_lines(f'{bedprefix}.bim'))) if snps is None else snps
    nsubj = count_lines(f'{bedprefix}.fam') if num_samples is None else num_samples
    n_bytes = nsubj//4
    if 4*n_bytes != nsubj:
        n_bytes += 1
    print(f'    {len(snp_list)} SNPs')
    print(f'    {nsubj} samples')

    genomat = np.zeros(shape=(len(snp_list), nsubj), dtype=np.int8)
    with open(f'{bedprefix}.bed', "rb") as bedid:
        for i, snp_index in enumerate(snp_list):
            bedid.seek(3+n_bytes*snp_index, 0)
            bed = np.frombuffer(bedid.read(n_bytes), dtype=np.uint8)
            genomat[i, :] = BYTE_GENO_MAP[bed].ravel()[:nsubj]
    return genomat

def write_bed(geno, bed_perm_path):
    n_snps = geno.shape[0]
    n_samples = geno.shape[1]
    magic = bytes([0x6c, 0x1b, 0x01])    
    geno_map = get_geno_byte_map()
    with open(bed_perm_path, 'wb') as bed_perm:
        bed_perm.write(magic)
        for i_snp in range(n_snps):
            snp_bytes = get_snp_bytes(geno[i_snp, :], geno_map)
            bed_perm.write(snp_bytes)
            if i_snp%1000 == 0:
                print(f'    {i_snp+1} SNPs processed')


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
    A = A.copy()
    B = B.copy()

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
