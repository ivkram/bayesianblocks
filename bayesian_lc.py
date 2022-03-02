import pandas as pd
import numpy as np
import csv
from astropy.stats import bayesian_blocks
from astropy.time import Time
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'font.size': 22})


def block_av(flux, epoch_, edges_bb_):
    """
    :param flux: initial fluxes
    :param epoch_: initial epochs
    :param edges_bb_: edges of bayesian blocks
    :return: array of BB fluxes as a result of averaging fluxes in each block
    """
    block_flux = np.zeros(len(edges_bb_) - 1)
    block_start = block_end = 0
    for j in range(1, len(edges_bb_)):
        while epoch_[block_end] < edges_bb_[j]:
            block_end += 1
        if j == len(edges_bb_) - 1:
            block_end += 1
        block_flux[j - 1] = np.average(flux[block_start:block_end])
        block_start = block_end
    return block_flux


def bayesian_plot(source_name_, params_, flux, flux_err, epoch_, flux_bb_, edges_bb_):
    plt.figure(figsize=(16, 9))
    plt.title(source_name_)
    plt.errorbar(epoch_, flux * 1E8, yerr=flux_err * 1E8, label='original LC', fmt='o',
                 elinewidth=1, markersize=2, alpha=0.3)
    plt.step(edges_bb_, np.append(flux_bb_, flux_bb_[-1]) * 1E8, where='post',
             label='bayesian blocks, ' + r'p0=' + str(p0), linewidth=2)
    plt.legend()
    plt.xlabel('Epoch (yr)')
    plt.ylabel("Photon flux (10$^{-8}$ ph cm$^{-2}$ s$^{-1}$)")
    plt.savefig((source_name_ + params_ + '.pdf'))
    plt.close()


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

    params = '_p0_' + str(p0)
    params += '_ts' + str(ts)
    params += '_npr' + str(n_pr)

    p = Path('.')
    for lc_file in sorted(p.glob('*.txt')):
        source_name = str(lc_file.parts[-1])[4:16]
        print(source_name)
        try:
            lc = pd.read_csv(lc_file,
                             usecols=[1, 2, 3, 4, 5, 7],
                             names=['start', 'end', 'TS', 'flux', 'flux_err', 'n_pred'],
                             dtype=np.float64,
                             delim_whitespace=True,
                             skiprows=1)
        except TypeError:
            print('Something went wrong while reading input file')
            quit()

        # apply TS filter
        lc = lc[lc.TS > ts].reset_index(drop=True)

        # apply N_predicted filter
        lc = lc[lc.n_pred > n_pr].reset_index(drop=True)

        # bayesian blocks algorithm
        epoch = (lc.start + lc.end) / 2
        edges_bb = bayesian_blocks(epoch, lc.flux, lc.flux_err, fitness='measures', p0=p0)
        flux_bb = block_av(lc.flux, epoch, edges_bb)

        # write to file
        with open(source_name + params, 'w', newline='') as bb:
            bb_writer = csv.writer(bb, delimiter=',')
            for i in range(len(flux_bb)):
                bb_writer.writerow([edges_bb[i], edges_bb[i + 1], flux_bb[i]])

        # draw plot
        if plot == 'y' or plot == '':
            # rescale epochs from mjd to decimal year
            epoch_dec = Time((epoch - 239557417) / 86400 + 54682, format='mjd').decimalyear
            edges_bb_dec = Time((edges_bb - 239557417) / 86400 + 54682, format='mjd').decimalyear

            # rescale flux
            lc.flux, lc.flux_err, flux_bb = map(lambda x: x * 1E8, (lc.flux, lc.flux_err, flux_bb))

            bayesian_plot(source_name, params, lc.flux, lc.flux_err, epoch_dec, flux_bb, edges_bb_dec)

    print('Success\n')
