import math
import matplotlib.pyplot as plt
import numpy as np
from os import listdir
import pandas as pd
import random
import re
from scipy.stats import norm as std_norm
import seaborn as sns
import tqdm


# plots a heat map of chromosome interaction intensities
def interaction_heat_map(all_files, data_path):
    biggest_num = max([max(chrm_nums(file)) for file in all_files])
    heat_matrix = np.zeros((biggest_num, biggest_num))

    # intensity = { (chr_x, chr_y): interaction intensity }
    intensity = count_intermingling(all_files, data_path)
    for pair in intensity:
        chr_x, chr_y = pair
        heat_matrix[chr_x-1, chr_y-1] = intensity[pair]
        heat_matrix[chr_y-1, chr_x-1] = intensity[pair]
    sns.heatmap(heat_matrix, cbar_kws={'label': 'Number of Intermingling Regions'})
    plt.title('Heat Map of Inter-Chromosome Interaction Counts')
    plt.xlabel('Chromosome number')
    plt.ylabel('Chromosome number')
    plt.show()


# maps a pair of chrms to the intensity of their interaction
def count_intermingling(all_files, data_path):
    n_intermingle = dict()
    mu, sigma = mean_std(all_files, data_path)
    print('\nCounting interminglings...')
    for file in tqdm.tqdm(all_files):
        x, y = chrm_nums(file)
        if x != y:
            # intensity score = sum of the rectangular areas with significant interaction freq
            n_intermingle[(int(x), int(y))] = sum([w*h for _, w, h in monte_carlo(x, y, mu, sigma, data_path)])
    return n_intermingle


# returns the chromosome numbers in the file name
def chrm_nums(file):
    chrx = re.search('chr(.+?)_', file).group(1)
    chry = re.search('_chr(.+?).txt', file).group(1)
    return int(chrx), int(chry)


# runs a greedy search 50 times to find the regions with lowest p-values
def monte_carlo(x, y, mu, sigma, data_path):
    file_name = f'{data_path}chr{x}_chr{y}.txt'
    Z = t_interaction_matrix(file_name)

    interacting = []
    while True:
        witnessed = [greedy_search(Z, mu, sigma) for _ in range(50)]
        # return the trial with the lowest p-value
        best_M, best_p, best_row, best_col = min(witnessed, key=lambda wit: wit[1])
        # the extended matrix is no longer significant
        if best_p > 0.01:
            return interacting
        else:
            h, w = best_M.shape
            # sub-matrix is an area with significant interaction
            interacting.append(((best_col, best_row), w, h))

            mean = np.nanmean(best_M)
            Z[best_row: best_row + h, best_col: best_col + w] -= mean


# fills a freq matrix, Z, where Z[i, j] = scaled_freq btwn chrm_i and chrm_j
def t_interaction_matrix(file):
    df = read(file)
    n_rows, n_cols = df['xloc'].max()+1, df['yloc'].max()+1
    freq_matrix = np.zeros((n_rows, n_cols))
    freq_matrix[df['xloc'], df['yloc']] = df['freq']
    return freq_matrix


# finds regions with small p-values
# uses adjusted p-value to perform a greedy search--instead of testing every sub-matrix
def greedy_search(freq_matrix, mu, sigma):
    n_rows, n_cols = freq_matrix.shape

    def adjusted_p_value(M):
        # returns the number of sub-matrices in an m x n matrix
        def count_submatrices(m, n):
            return m * (m + 1) * n * (n + 1) / 4
        top = (np.nanmean(M) - mu) * math.sqrt(M.size)
        n_sub = count_submatrices(n_rows, n_cols)
        # adjust by the number of possible sub-matrices
        return std_norm.sf(top / sigma) * n_sub

    # the region begins as a random position
    row, col = random.randint(0, n_rows - 1), random.randint(0, n_cols - 1)
    rand_choice = freq_matrix[row, col]
    M_0 = np.array([[rand_choice]])

    # continue extending the region until the p-value starts increasing
    while True:
        p_value = adjusted_p_value(M_0)
        p_values = []
        for sub_m in extend(row, col, M_0, freq_matrix):
            p_values.append((adjusted_p_value(sub_m[0]), sub_m))
        best_p, best_m = min(p_values, key=lambda x: x[0])
        if best_p >= p_value:
            return M_0, p_value, row, col
        else:
            M_0, row, col = best_m


# extends the sub-matrix M to the right, left, top, and bottom
def extend(r, c, M, matrix):
    h, w = M.shape
    M_r, rr, cr = matrix[r: r + h, c: c + w + 1].copy(), r, c
    M_l, rl, cl = matrix[r: r + h, c - 1: c + w].copy(), r, c - 1
    M_t, rt, ct = matrix[r - 1: r + h, c: c + w].copy(), r - 1, c
    M_b, rb, cb = matrix[r: r + h + 1, c: c + w].copy(), r, c
    return (M_r, rr, cr), (M_l, rl, cl), (M_t, rt, ct), (M_b, rb, cb)


# aggregates the interaction frequencies among all files, and computes the mean and std
def mean_std(all_files, data_path):
    all_freqs = pd.Series([])
    # files are sparse, so must account for the extra zeros
    total_size = 0
    print('\nCalculating mean and std...')
    for file in tqdm.tqdm(all_files):
        chrx, chry = chrm_nums(file)
        if chrx != chry:
            chr_df = read(data_path + file)
            file_freqs = chr_df.freq
            all_freqs = pd.concat([all_freqs, file_freqs])
            total_size += max(chr_df.xloc) * max(chr_df.yloc)

    # calculate mean and variance of sparse data
    mean = all_freqs.sum() / total_size

    se_present = ((all_freqs - mean)**2).sum()
    n_sparse = total_size - len(all_freqs)
    se_missing = mean**2 * n_sparse
    variance = (se_present + se_missing) / (total_size-1)

    print(f'# data points = {len(all_freqs)}')
    return mean, variance**0.5


# reads and scales the data
def read(file):
    def transform(data):
        return np.log(1 + data)

    df = pd.read_csv(file, header=None, sep='\t')
    df.columns = ['xloc', 'yloc', 'freq']
    # scale by the data's resolution
    norm = 250000
    df[['xloc', 'yloc']] = df[['xloc', 'yloc']] // norm
    df['freq'] = transform(df['freq'])
    return df


if __name__ == '__main__':
    data_path = './data/'
    all_files = listdir(data_path)
    interaction_heat_map(all_files, data_path)
