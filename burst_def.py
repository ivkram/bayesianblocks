import numpy as np
import math

def burst_def(edgesBB, fluxBB):
    #global average
    global_average = 0
    for counter, flux in enumerate(fluxBB):
        global_average += flux * (edgesBB[counter+1] - edgesBB[counter])
    global_average /= (edgesBB[-1] - edgesBB[0])
    
    flux_lowstate = fluxBB
    time_lowstate = edgesBB[1:] - edgesBB[:-1]
    time = time_lowstate
    for i in range(len(fluxBB)):
        average = 0
        for counter, flux in enumerate(flux_lowstate):
            average += flux * time_lowstate[counter]
        average /= sum(time_lowstate)

        anti_av = 0
        summ = 0
        for i in range(len(fluxBB)):
            if (fluxBB[i] not in flux_lowstate) and (fluxBB[i] > average):
                anti_av += fluxBB[i] * time[i]
                summ += time[i]
        if anti_av != 0:
            anti_av = anti_av / summ
        else:
            anti_av = average

        disp = 0
        m = 0
        n = 0
        while n < len(flux_lowstate):
            if abs(flux_lowstate[n] - average) > disp:
                disp = abs(flux_lowstate[n] - average)
                m = n
            n += 1
        if disp/anti_av < 0.25:
            break
        flux_lowstate = np.delete(flux_lowstate, m)
        time_lowstate = np.delete(time_lowstate, m)

    #dispersion of low state
    true_disp = 0.0
    for j in range(len(flux_lowstate)):
        true_disp += ((flux_lowstate[j] - average) ** 2) * time_lowstate[j]
    true_disp = math.sqrt(true_disp/sum(time_lowstate))

    return average, true_disp
