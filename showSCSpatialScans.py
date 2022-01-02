import numpy as np
import matplotlib.pyplot as plt
import cantrips as can
import sys
import AstronomicalParameterArchive as apa
import scipy

def readInSCandPDDataFile(data_file, data_dir, ignore_str = 'DELETED'):
    with open (data_dir + data_file, 'r') as open_file:
        data = open_file.readlines()
    data = [row.split('|') for row in data]
    positions = [float(row[0])  for row in data]
    SC_illum = [[ ( float(elem)) for elem in can.recursiveStrToListOfLists(row[1]) if elem != ignore_str] for row in data]

    return [positions, SC_illum ]

if __name__=="__main__":
    dir_root = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/solarCell/spatialUniformity/'
    #data_files = [ stem + date_str + '.txt' for stem in ['SC_vs_PD_QE_with_B2987A_cellPosA_inten2_10_', 'SC_vs_PD_QE_with_B2987A_cellPosA_inten4_10_', 'SC_vs_PD_QE_with_B2987A_cellPosB_inten2_10_', 'SC_vs_PD_QE_with_B2987A_cellPosB_inten4_10_'] ]

    stems = ['SC_AX_scanPerpElectrodes_10umPinhole_', 'SC_HG_scanPerpElectrodes_10umPinhole_'  ]
    suffixes = ['', '']
    date_strs = ['20210531', '20210531']
    #extra_title_str = 'Spatial'
    cell_ids = ['Cell AX (old generation) - ' + date_strs[0], 'Cell HG (new generation) - ' + date_strs[0] ]
    shunt_resistances = [610, 610]
    pos_sep = 0.1

    SC_neg_current = 0
    PD_neg_current = 1
    deleted_element_id = 'DELETED'
    data_files = [ stems[i] + date_strs[i] + suffixes[i] + '.txt' for i in range(len(stems)) ]
    file_roots = [f[0:-len('.txt')] for f in data_files]
    save_file_names = [root + '.pdf' for root in file_roots]
    xlims = [-5, 105]
    ylims = [0.89, 1.01]
    recorded_bias_voltages = 1
    title_font_size = 24
    labelsize = 16
    ticksize = 14
    plt.rc('xtick',labelsize=ticksize)
    plt.rc('ytick',labelsize=ticksize)
    data_dirs = [dir_root + date_str + '/' for date_str in date_strs]


    data_sets = [readInSCandPDDataFile(data_files[i], data_dirs[i],  ignore_str = deleted_element_id) for i in range(len(data_dirs))]
    print ('len(data_sets[0][0]) = ' + str(len(data_sets[0][0])))

    f, axarr = plt.subplots(len(data_sets), 2, figsize = [14, 5 * len(data_sets)], squeeze = False)
    plt.subplots_adjust(hspace=0.00)
    for i in range(len(data_sets)):
        data_set = data_sets[i]
        positions = data_set[0]
        currents = data_set[1]
        positions = positions[60:1060]
        currents = currents[60:1060]
        mean_currents = np.mean(currents, axis = 1) * 10.0 ** 6.0
        normalized_photocurrents = mean_currents / np.max(mean_currents)
        std_currents = np.std(currents, axis = 1) * 10.0 ** 6.0
        current_fft = scipy.fft.fft(normalized_photocurrents - np.median(normalized_photocurrents))
        current_abs_fft = np.abs( current_fft[1:len(normalized_photocurrents)//2] )
        pos_freqs = scipy.fft.fftfreq(len(mean_currents), pos_sep)[1:len(normalized_photocurrents)//2]

        print ('len(mean_currents) = ' + str(len(mean_currents)))
        print ('[len(positions), len(mean_currents)] = ' + str([len(positions), len(normalized_photocurrents)]))
        photocurrent = axarr[i,0].plot(np.array(positions) - positions[0], normalized_photocurrents , color = 'k')[0]
        #axarr[i,0].errorbar(positions, mean_currents, yerr = std_currents, fmt = 'none')
        axarr[i,0].set_xlabel('Relative spot position (mm)', fontsize = labelsize)
        axarr[i,0].set_ylabel(r'SC Normalized Photocurrent', fontsize = labelsize)
        axarr[i,0].set_xlim(xlims)
        axarr[i,0].set_ylim(ylims)
        axarr[i,0].legend([photocurrent], [cell_ids[i]], fontsize = labelsize,  bbox_to_anchor=[0.5, 0.95], loc='center')

        axarr[i,1].plot(pos_freqs, current_abs_fft, color = 'k')
        axarr[i,1].set_xlabel(r'Wave-number (mm$^{-1}$)', fontsize = labelsize)
        axarr[i,1].set_ylabel(r'Fourier Amplitude of Photocurrent', fontsize = labelsize)
        axarr[i,1].legend([photocurrent], [cell_ids[i]], fontsize = labelsize,  bbox_to_anchor=[0.5, 0.95], loc='center' )

        first_peak_index = np.argmax(current_abs_fft[50:200])
        first_peak_loc = pos_freqs[50:200][first_peak_index]
        first_peak = current_abs_fft[50:200][first_peak_index]
        first_peak_precision = pos_freqs[50:200][first_peak_index + 1] - pos_freqs[50:200][first_peak_index]
        print (r'[first_peak_loc $\pm$ first_peak_precision, first_peak] = ' + str([first_peak_loc, first_peak_precision, first_peak]))


    plt.tight_layout()
    plt.subplots_adjust(hspace=0.02)
    plt.savefig(data_dirs[-1] + file_roots[-1] + '.pdf')
    plt.show()
