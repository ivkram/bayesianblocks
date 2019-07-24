from astropy import stats
import numpy as np
import matplotlib.pyplot as plt
import glob

SNR = 0 #signal-to-noise ratio
output_format = 'png'

path = '*.txt'
files = glob.glob(path)

for src in files:    
    curve = open(src, 'r')

    print('\n'+ src[:-4])
    t, f, ferr = ([] for i in range(3)) #time, flux, flux error
    for line in curve:
        scan = np.array(line.split())
        scan = scan.astype('float')
        if (scan[1] > SNR * scan[2] and scan[2] > 0):
            t.append(scan[0]); f.append(scan[1]); ferr.append(scan[2])
        if scan[2] == 0:
            print('zero f_err! {} {}'.format(scan[0], scan[1]))

    if len(t) < 2:
        continue
    
    #bayesian blocks algorithm
    edgesBB = stats.bayesian_blocks(t, f, ferr, fitness='measures', p0=5.7e-7)

    #average flux in block
    fluxBB = []
    for i in range(len(edgesBB)-1):
        fluxBBi = 0; nbbi = 0
        for j in range(0,len(f)):
            if (edgesBB[i] <= t[j] < edgesBB[i+1]) or (j == len(f)-1):
                fluxBBi += f[j]; nbbi += 1
        try:
            fluxBB.append(fluxBBi/nbbi)
        except ZeroDivisionError:
            print('block {:d} is somehow empty'.format(i+1))
            fluxBB.append(0)
            
        #here we try to calculate the weighted mean
        '''
        fluxBBi = 0; sigma = 0
        for j in range(0,len(f)):
            if (edgesBB[i] <= t[j] and t[j] < edgesBB[i+1]) or (j == len(f)-1):
                sigma += 1.0 / (ferr[j] ** 2)
        for j in range(0,len(f)):
            if (edgesBB[i] <= t[j] and t[j] < edgesBB[i+1]) or (j == len(f)-1):
                fluxBBi += 1.0 / (ferr[j] ** 2) / sigma * f[j]
        fluxBB.append(fluxBBi)
        '''
    
    print(edgesBB, fluxBB)
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    
    ax1.errorbar(t, f, yerr = ferr, fmt = 'ok', ecolor='black', fillstyle = 'none',
                 markersize = 5, markeredgecolor = 'black', linewidth = 1,
                 alpha = 0.3, zorder = 0)
    
    ax1.step(edgesBB, [fluxBB[0]] + fluxBB, color = 'steelblue', linewidth = 2, zorder = 10)
    
    plt.xlabel("Epoch (yr)")
    plt.ylabel("F$_\mathrm{0.1-300 GeV}$ (10$^{-8}$ ph cm$^{-2}$ s$^{-1}$)")
    ax1.set_title(src[:-4], loc = 'right', fontsize = 10)
    #for express view use plt.show()
    plt.savefig(src[:-3] + output_format, bbox_inches = 'tight')
    plt.close()
    
    output = open(src[:-4] + '_bb', 'w')
    for i in range(len(fluxBB)):
        output.write('{:.4f} {:.4f} {:9.5f}\n'.format(edgesBB[i], edgesBB[i+1], fluxBB[i]))
    output.close()
    
    curve.close()
