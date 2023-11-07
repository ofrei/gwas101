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


def read_bed(bfile_orig, num_snps=None, num_samples=None):
    n_snps = count_lines(f'{bfile_orig}.bim') if num_snps is None else num_snps
    n_samples = count_lines(f'{bfile_orig}.fam') if num_samples is None else num_samples
    n_cols = n_samples//4
    if 4*n_cols != n_samples:
        n_cols += 1
    print(f'    {n_snps} SNPs')
    print(f'    {n_samples} samples')
    bed = np.fromfile(f'{bfile_orig}.bed', dtype=np.uint8)[3:]
    geno = BYTE_GENO_MAP[bed].reshape(n_snps, 4*n_cols)
    if (n_samples % 4) != 0:
        geno = geno[:, :n_samples]
    return geno

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

