import matplotlib.pyplot as plt
import numpy as np

import os
import glob

from astropy.time import Time
from astropy import stats

import find
import csv
import math

'''
Для успешного выполнения скрипта надо прописать:
1) SNR, TS пороги
2) Добавлять ли OVRO (если да, то прописать название папки -- по умолчанию 'OVRO')
3) Сохранять ли картинки (если да, то прописать название папки сохранения)
4) Относительное расположение radio и gamma данных, файла с ассоциациями
'''

#порог на SNR
SNR = 0.0
#порог на TS
TS  = 0.0

#добавляем ли OVRO
DRAW_OVRO = False
#сохранять ли картинки (иначе -- потоковый вывод)
SAVE = False
#папка для сохранения
SAVE_FOLDER = 'r vs g light curves (bayesian)'

fileDir = os.path.dirname(os.path.realpath('__file__'))

#файлик с радио
radio_data = open(fileDir + '/Radio_Gamma/MOJAVE_LC/source_epoch_flux', 'r')

#файлики с гамма
path = fileDir + '/Radio_Gamma/Fermi_LC/weekly/' + '*.txt'
files = glob.glob(path)
files.sort()

#файл с просписанными ассоциациями
assoc_file = open(fileDir + '/Radio_Gamma/Fermi_LC/weekly/radio_gamma_names', 'r')

#папка с ovro
files_ovro = glob.glob(fileDir + '/OVRO/' + '*.csv')
files_ovro.sort()


'''
сборка данных OVRO
'''
def draw_ovro(name, start, end):
    ovro_src = find.get_J2000(assoc_file, name)
    print(ovro_src)
    for src in files_ovro:
        if ovro_src in src:
            print('OVRO:  ' + ovro_src)
            ovro_data = open(src, 'r')
            ovro_data.readline()

            source = []
            for line in ovro_data:
                obs = np.array(line.split(','))
                obs = obs.astype(float)
                obs_time = Time(obs[0], format = 'mjd')
                obs[0] = obs_time.decimalyear
                if (start < obs[0] < end):
                    source.append(obs)
            ovro_data.close()
            return np.array(source)
    return []


'''
главная функция отрисовки
'''
def draw_lc(name, z, r_source):

    #find association
    print(name)
    gamma_src = find.find_gamma(assoc_file, name)

    #finally plot gamma vs radio
    for src in files:
        if (gamma_src in src):
            print('GAMMA: ' + gamma_src)
            fig, (ax1, ax2) = plt.subplots(2, sharex=True)
            fig.set_size_inches(11.69, 8.27)
            
            #get gamma data
            gamma_data = open(src, 'r')
            gamma_data.readline()
            
            g_source = []
            for counter, line in enumerate(gamma_data):
                obs = np.array(line.split())[1:6]
                obs = obs.astype(float)
                obs[0] = start_0[counter]
                obs[1] = end_0[counter]
                obs[3] = obs[3] * (10**8)
                obs[4] = obs[4] * (10**8)
                if (obs[2] > TS) and (obs[3]/obs[4] > SNR):
                    g_source.append(obs)
            g_source = np.array(g_source)
            gamma_data.close()
            
            g_flux    = g_source[:,3]
            g_flux_err= g_source[:,4]
            start = g_source[:,0]
            end = g_source[:,1]

            #add gamma
            snr = g_source[:, 3]/g_source[:, 4]
            g_source = g_source[snr > 0.05, :]
            true_flux = g_source[g_source[:, 2] > 4, :]
            upper_lim = g_source[g_source[:, 2] <= 4, :]

            ax2.errorbar( (true_flux[:,0]+true_flux[:,1])/2, true_flux[:,3], yerr = true_flux[:,4], fmt = 'ok', fillstyle = 'full',
                          markersize = 4, markeredgecolor = 'black', linewidth = 1, alpha = 0.8)
            ax2.errorbar( (upper_lim[:,0]+upper_lim[:,1])/2, upper_lim[:,3], yerr = upper_lim[:,4], fmt = 'vk', fillstyle = 'none',
                          markersize = 4, markeredgecolor = 'black', linewidth = 1, alpha = 0.8)
            ax2.set_xlabel('Epoch (yr)')
            ax2.set_ylabel("F$_\mathrm{0.1-300 GeV}$ (10$^{-8}$ ph cm$^{-2}$ s$^{-1}$)")
            ax2.set_title(src[-38:-22], loc = 'right', fontsize = 10)
            
            #bayesian
            
            t = (start+end)/2
            edgesBB = stats.bayesian_blocks(t, g_flux, g_flux_err, fitness='measures', p0 = 1.0e-3)
            
            fluxBB = []     #'average' flux in block
            BB_start = 0
            for i in range(len(edgesBB)-1):
                BB_end = BB_start
                while t[BB_end] <= edgesBB[i + 1]:
                       BB_end += 1
                       if BB_end == len(g_flux):
                           break
                fluxBB.append(np.median(g_flux[BB_start:BB_end]))
                BB_start = BB_end
            
            #define burst

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

            #draw threshold
            limit = average + 3 * true_disp
            ax2.plot([t[0], t[-1]], [limit, limit], color = '#da1ba5', zorder = 9)
            fluxBB_max = max(fluxBB)
            #draw bursts
            for counter, flux in enumerate(fluxBB):
                if (flux > limit):
                    ax1.axvspan(edgesBB[counter], edgesBB[counter+1], alpha=0.4*fluxBB[counter]/fluxBB_max, color='steelblue', zorder = 2)

            #add bayesian
            fluxBB = np.array([fluxBB[0]] + fluxBB)
            ax2.step(edgesBB, fluxBB, color = 'steelblue', linewidth = 2, zorder = 10)
            
            #moving radio borders
            left = 0; right = len(r_source)
            while r_source[left][0] < start[0] - 1:
                left += 1
                if left == len(r_source):
                    break
            while r_source[right-1][0] > end[-1] + 1:
                right -= 1
                if (right == 0):
                    break

            r_epoch = r_source[left:right, 0]
            core  = r_source[left:right, 1]
            jet   = r_source[left:right, 2]

            #add radio total and core
            ax1.scatter(r_epoch, jet + core, c = 'black', s = 7, zorder = 4, label='total')
            ax1.scatter(r_epoch, core, c = 'red', s = 7, zorder = 4, label='core')
            ax1.set_title('MOJAVE ' + name, loc = 'right', fontsize = 10)
            ax1.set_ylabel('F$_{15 GHz}$ (Jy)')
            ax1.legend()

            #add OVRO
            if DRAW_OVRO == True: 
                o_source = draw_ovro(name, start[0] - 1, end[-1] + 1)
                if (len(o_source) > 0):
                    ax1.plot(o_source[:,0], o_source[:,1], linewidth = 1, color = 'grey', alpha = 0.4, label = 'OVRO', zorder = 3)

            #start radio from zero
            bottom, top = ax1.get_ylim()
            ax1.set_ylim(0, top)
            
            if SAVE == True:
                plt.savefig(SAVE_FOLDER + '/' + name + ' vs ' + src[-38:-22] + ' (bayesian).png', dpi = 300, bbox_inches = 'tight')
                print('saved successfully')
            else:
                plt.show()
            plt.close(fig)

'''
посокльку некоторые времена 'битые', то полагаемся на  на самый широкий диапазон времен -- 1й источник!!!!
make basis start end times in gamma
'''
gamma_data = open(files[0], 'r')
gamma_data.readline()
start_0, end_0 = ([] for i in range(2))
for line in gamma_data:
    obs = line.split()

    #time to decimal year
    obs_time1 = Time((float(obs[1])/86400.0+2451544.5+365), format='jd')
    obs_time2 = Time((float(obs[2])/86400.0+2451544.5+365), format='jd')
    start_0.append(obs_time1.decimalyear)
    end_0.append(obs_time2.decimalyear)

gamma_data.close()


'''
reading radio and making plots
'''
name = ''
for line in radio_data:
    if not line.startswith('#'):
        obs = np.array(line.split())
        
        #change date to decimal
        obs_time = Time(obs[1])
        obs[1] = obs_time.decimalyear

        if obs[0] != name:
            if name != '':
                if len(source) > 4:
                    if name != '2037+511': #and (int(name[:2]) > 10):
                        source = source.astype(float)
                        draw_lc(name, z, source)
                    print()
            source = np.array([obs[1:-1]])
            name = obs[0]
            z = obs[-1]
        else:
            source = np.concatenate((source, [obs[1:-1]]))
# для последнего источника
source = source.astype(float)
draw_lc(name, z, source)

radio_data.close()
assoc_file.close()
