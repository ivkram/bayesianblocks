import pandas as pd
import numpy as np
import csv
from astropy.stats import bayesian_blocks
from astropy.time import Time
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'font.size': 22})


def block_av(lc_, edges):
    """
    :param lc_: original light curve
    :param edges: edges of bayesian blocks
    :return: array bayesian blocks fluxes
    """
    block_flux = np.zeros(len(edges) - 1)
    block_start = block_end = 0
    for i in range(1, len(edges)):
        while times[block_end] < edgesBB[i]:
            block_end += 1
        if i == len(edges) - 1:
            block_end += 1
        block_flux[i - 1] = np.average(lc_.loc[block_start:(block_end - 1), 'flux'])
        block_start = block_end
    return block_flux


if __name__ == '__main__':
    print('--- Bayesian blocks ---\n')

    try:
        p0 = float(input("False alarm probability (0.05 by default): "))
    except ValueError:
        p0 = 0.05
    else:
        if not (0 < p0 < 1):
            p0 = 0.05

    try:
        ts = float(input("TS threshold (4 by default): "))
    except ValueError:
        ts = 4

    try:
        n_pr = float(input("N predicted threshold (10 by default): "))
    except ValueError:
        n_pr = 10

    plot = input('Draw plot? (y/n): ')
    print('\n')

    title = '_p0_' + str(p0)
    title += '_ts' + str(ts)
    title += '_npr' + str(n_pr)

    p = Path('.')
    for lc_file_ext in list(p.glob('*.txt')):
        lc_file = str(lc_file_ext.parts[-1])[4:16]
        print(lc_file)
        try:
            lc = pd.read_csv(lc_file_ext,
                             usecols=[1, 2, 3, 4, 5, 7],
                             names=['start', 'end', 'TS', 'flux', 'flux_err', 'n_pred'],
                             dtype=np.float64,
                             delim_whitespace=True,
                             skiprows=1)
        except TypeError:
            print('Something went wrong while reading input file')
            exit()

        # apply TS filter
        lc = lc[lc.TS > ts].reset_index(drop=True)

        # apply N_predicted filter
        lc = lc[lc.n_pred > n_pr].reset_index(drop=True)

        # bayesian blocks algorithm
        times = (lc.start + lc.end) / 2
        edgesBB = bayesian_blocks(times, lc.flux, lc.flux_err, fitness='measures', p0=p0)
        fluxBB = block_av(lc, edgesBB)

        # write to file
        with open(lc_file + title, 'w', newline='') as bb:
            bb_writer = csv.writer(bb, delimiter=',')
            for i in range(len(fluxBB)):
                bb_writer.writerow([edgesBB[i], edgesBB[i + 1], fluxBB[i]])

        # draw plot
        if plot == 'y' or plot == '':
            plt.figure(figsize=(16, 9))
            times_yr = Time((times - 239557417) / 86400 + 54682, format='mjd').decimalyear
            edgesBB_yr = Time((edgesBB - 239557417) / 86400 + 54682, format='mjd').decimalyear
            plt.errorbar(times_yr, lc.flux * 1E8, yerr=lc.flux_err * 1E8, label='original LC', fmt='o',
                         elinewidth=1, markersize=2, alpha=0.3)
            plt.title(lc_file)
            plt.step(edgesBB_yr, np.append(fluxBB, fluxBB[-1]) * 1E8, where='post',
                     label='bayesian blocks, ' + r'p0=' + str(p0), linewidth=2)
            plt.legend()
            plt.xlabel('Epoch (yr)')
            plt.ylabel("Photon flux (10$^{-8}$ ph cm$^{-2}$ s$^{-1}$)")
            plt.savefig((lc_file + title + '.pdf'))
            plt.close()

    print('Success\n')
