import numpy as np
import matplotlib.pyplot as plt
import cantrips as can
import sys
import AstronomicalParameterArchive as apa
import scipy

def readInSCandPDDataFile(data_file, data_dir, recorded_bias_voltages = 1 ):
    with open (data_dir + data_file, 'r') as open_file:
        data = open_file.readlines()
    data = [row.split('|') for row in data]
    LEDs = [row[0]  for row in data]
    PD_illum = [[ (float(elem)) for elem in can.recursiveStrToListOfLists(row[1])] for row in data]
    PD_dark = [[(float(elem)) for elem in can.recursiveStrToListOfLists(row[2])] for row in data]
    SC_illum = [[(float(elem)) for elem in can.recursiveStrToListOfLists(row[3])] for row in data]
    SC_dark = [[(float(elem)) for elem in can.recursiveStrToListOfLists(row[4])] for row in data]
    if recorded_bias_voltages:
        PD_illum_bv = [[ (float(elem)) for elem in can.recursiveStrToListOfLists(row[5])] for row in data]
        PD_dark_bv = [[(float(elem)) for elem in can.recursiveStrToListOfLists(row[6])] for row in data]
        SC_illum_bv = [[(float(elem)) for elem in can.recursiveStrToListOfLists(row[7])] for row in data]
        SC_dark_bv = [[(float(elem)) for elem in can.recursiveStrToListOfLists(row[8])] for row in data]
    #SC_illum_bv = [[ (float(elem[1:-1])) for elem in can.recursiveStrToListOfLists(row[4])] for row in data]
    #SC_dark_bv = [[(float(elem[1:-1])) for elem in can.recursiveStrToListOfLists(row[5])] for row in data]
    #PD_illum_bv = [[(float(elem[1:-1])) for elem in can.recursiveStrToListOfLists(row[6])] for row in data]
    #PD_dark_bv = [[(float(elem[1:-1])) for elem in can.recursiveStrToListOfLists(row[7])] for row in data]

    #dial_positions, bias_voltages, SC_illum, SC_dark, PD_illum_bv, PD_dark_bv, SC_illum_bv, SC_dark_bv, PD_illum_bv, PD_dark_bv = can.safeSortOneListByAnother(bias_voltages, [zaber_positions, bias_voltages, SC_illum, SC_dark, PD_illum_bv, PD_dark_bv, SC_illum_bv, SC_dark_bv, PD_illum_bv, PD_dark_bv])
    #return [dial_positions, bias_voltages, SC_illum, SC_dark, PD_illum_bv, PD_dark_bv, SC_illum_bv, SC_dark_bv, PD_illum_bv, PD_dark_bv]
    if recorded_bias_voltages:
        #wavelengths, SC_illum, SC_dark, PD_illum, PD_dark, SC_illum_bv, SC_dark_bv, PD_illum_bv, PD_dark_bv = can.safeSortOneListByAnother(wavelengths, [wavelengths, SC_illum, SC_dark, PD_illum, PD_dark, SC_illum_bv, SC_dark_bv, PD_illum_bv, PD_dark_bv])
        return [LEDs, SC_illum, SC_dark, PD_illum, PD_dark, SC_illum_bv, SC_dark_bv, PD_illum_bv, PD_dark_bv]
    else:
        #wavelengths, SC_illum, SC_dark, PD_illum, PD_dark = can.safeSortOneListByAnother(wavelengths, [wavelengths, SC_illum, SC_dark, PD_illum, PD_dark])
        return [LEDs, SC_illum, SC_dark, PD_illum, PD_dark]

def getSCandPDStats(cal_illum, cal_dark, uncal_illum, uncal_dark):
    cal_means, cal_std = [[np.mean(row) for row in cal_illum ], [np.std(row) for row in cal_illum ] ]
    cal_dark_means, cal_dark_std = [[np.mean(row) for row in cal_dark ], [np.std(row) for row in cal_dark ] ]

    uncal_means, uncal_std = [[np.mean(row) for row in uncal_illum ], [np.std(row) for row in uncal_illum ] ]
    uncal_dark_means, uncal_dark_std = [[np.mean(row) for row in uncal_dark ], [np.std(row) for row in uncal_dark ] ]

    cal_diffs = [cal_means[i] - cal_dark_means[i] for i in range(len(cal_means))]
    cal_diff_stds = [np.sqrt(cal_std[i] ** 2.0 + cal_dark_std[i] ** 2.0) for i in range(len(cal_std))]
    uncal_diffs = [uncal_means[i] - uncal_dark_means[i]  for i in range(len(uncal_means))]
    uncal_diff_stds = [np.sqrt(uncal_std[i] ** 2.0 + uncal_dark_std[i] ** 2.0) for i in range(len(cal_std))]

    return [cal_diffs, cal_diff_stds, uncal_diffs, uncal_diff_stds]


if __name__ == "__main__":
    date_str = '20210505'
    dir_root = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/solarCell/qe/'
    #data_files = [ stem + date_str + '.txt' for stem in ['SC_vs_PD_QE_with_B2987A_cellPosA_inten2_10_', 'SC_vs_PD_QE_with_B2987A_cellPosA_inten4_10_', 'SC_vs_PD_QE_with_B2987A_cellPosB_inten2_10_', 'SC_vs_PD_QE_with_B2987A_cellPosB_inten4_10_'] ]
    stems = ['SC_vs_PD_QE_from_LEDs_with_B2987A_']
    suffixes = ['_A']
    data_files = [ stems[i] + date_str + suffixes[i] + '.txt' for i in range(len(stems)) ]
    file_roots = [f[0:-len('.txt')] for f in data_files]
    save_file_names = [root + '.pdf' for root in file_roots]
    time_sep = 1/60 * 10
    xlims = [330, 1100]
    n_burn_in = 0
    current_scaling_A_to_uA = 10.0 ** 6.0
    recorded_bias_voltages = 0
    leds_to_wavelength_dict = {'M365FP1':365, 'M470F3':470, 'MINTF4':554, 'M625F2':625, 'M660FB1':660, 'M850F2':850, 'M940F3':940, 'M970F3':970, 'M1050F3':1050 }
    wavelength_to_led_dict = {365:'M365FP1', 470:'M470F3', 554:'MINTF4', 625:'M625F2', 660:'M660FB1', 850:'M850F2', 940:'M940F3', 970:'M970F3', 1050:'M1050F3' }

    data_dir = dir_root + date_str + '/'
    PD_QE_data_file = dir_root + 'Hamamatsu_Photodiode_S2281_Spectral_Power_Response.txt'
    save_plot_file = 'SC_QE_measurement_' + date_str + '.pdf'
    PD_dark_plot_suffix, PD_ill_plot_suffix, SC_dark_plot_suffix, SC_ill_plot_suffix = ['_PD_dark_photocurrents.pdf', '_PD_ill_photocurrents.pdf', '_SC_dark_photocurrents.pdf', '_SC_ill_photocurrents.pdf']
    PD_QE_data = can.readInColumnsToList(PD_QE_data_file, delimiter = ' ', n_ignore = 1)
    PD_QE_wavelengths, PD_QE_A_per_W= [ [float(wave) for wave in PD_QE_data[0]], [float(qe) for qe in PD_QE_data[1]] ]
    astro_arch = apa.AstronomicalParameterArchive()
    e_charge = astro_arch.getElectronCharge()  #electron charge in coloumbs
    h_JouleS = astro_arch.getPlancksConstant()  #Planck's constant, in Joule X second
    c = astro_arch.getc() #speed of light in km/s
    PD_QE_e_per_phot = [PD_QE_A_per_W[i] / e_charge * h_JouleS * c / (PD_QE_wavelengths[i] * 10.0 ** -12 ) for i in range(len(PD_QE_wavelengths)) ]
    PD_QE_interp =  scipy.interpolate.interp1d(PD_QE_wavelengths, PD_QE_e_per_phot)

    bad_data_indeces = [[ ] for f in data_files]
    index_break = []
    file_legend_ids = ['']
    PD_color, SC_color, combined_color = ['b','r','purple']
    plot_styles = ['-','--']
    SC_plots = []
    PD_plots = []
    SC_QE_plots = []
    PD_QE_plots = []
    ratio_plots = []

    for i in range(len(data_files)):
        f_QE, axarr_QE = plt.subplots(2,2, figsize = [10, 10] )
        data_file = data_files[i]
        file_root = file_roots[i]

        #z_positions, bias_voltages, SC_illum, SC_dark, PD_illum_bv, PD_dark_bv, SC_illum_bv, SC_dark_bv, PD_illum_bv, PD_dark_bv = readInSCandPDDataFile(data_file, data_dir)
        if recorded_bias_voltages:
            LEDs, SC_illum, SC_dark, PD_illum, PD_dark, SC_illum_bv, SC_dark_bv, PD_illum_bv, PD_dark_bv = readInSCandPDDataFile(data_file, data_dir, recorded_bias_voltages = recorded_bias_voltages )
        else:
            LEDs, SC_illum, SC_dark, PD_illum, PD_dark = readInSCandPDDataFile(data_file, data_dir, recorded_bias_voltages = recorded_bias_voltages )
        wavelengths = [leds_to_wavelength_dict[LED] for LED in LEDs]
        full_SC_illum, full_SC_dark, full_PD_illum, full_PD_dark = [np.array([wave_set for wave_set in SC_illum]) * -10.0 ** 6.0, np.array([wave_set for wave_set in SC_dark]) * -10.0 ** 6.0, np.array([wave_set for wave_set in PD_illum]) * -10.0 ** 6.0, np.array([wave_set for wave_set in PD_dark]) * -10.0 ** 6.0]
        SC_illum, SC_dark, PD_illum, PD_dark = [[wave_set[n_burn_in:] for wave_set in arr] for arr in [full_SC_illum, full_SC_dark, full_PD_illum, full_PD_dark]]
        n_single_shot_cols = int(np.ceil(len(wavelengths) / 4))
        f_PD_ill, axarr_PD_ill = plt.subplots(4, n_single_shot_cols, figsize = (n_single_shot_cols * 2, 4 * 2), sharex = True)
        f_PD_dk, axarr_PD_dk = plt.subplots(4, n_single_shot_cols, figsize = (n_single_shot_cols * 2, 4 * 2), sharex = True)
        f_SC_ill, axarr_SC_ill = plt.subplots(4, n_single_shot_cols, figsize = (n_single_shot_cols * 2, 4 * 2), sharex = True)
        f_SC_dk, axarr_SC_dk = plt.subplots(4, n_single_shot_cols, figsize = (n_single_shot_cols * 2, 4 * 2), sharex = True)
        for j in range(len(wavelengths)):
            wave = wavelengths[j]
            row_num, col_num = [int(j // n_single_shot_cols), int(j % n_single_shot_cols)]
            axarr_SC_ill [row_num, col_num].scatter([time_sep * t for t in range(len(full_SC_illum[j]))], full_SC_illum[j], marker = '.', c = 'k')
            axarr_SC_dk [row_num, col_num].scatter([time_sep * t for t in range(len(full_SC_dark[j]))], full_SC_dark[j], marker = '.', c = 'k')
            axarr_PD_ill [row_num, col_num].scatter([time_sep * t for t in range(len(full_PD_illum[j]))], full_PD_illum[j], marker = '.', c = 'k')
            axarr_PD_dk [row_num, col_num].scatter([time_sep * t for t in range(len(full_PD_dark[j]))], full_PD_dark[j], marker = '.', c = 'k')
            axarr_PD_ill [row_num, col_num].text(0.1, 0.9, r'$\lambda=$' + str(wave) + 'nm', transform=axarr_PD_ill [row_num, col_num].transAxes, color = 'r')
            axarr_PD_dk [row_num, col_num].text(0.1, 0.9, r'$\lambda=$' + str(wave) + 'nm', transform=axarr_PD_dk [row_num, col_num].transAxes, color = 'r')
            axarr_SC_ill [row_num, col_num].text(0.1, 0.9, r'$\lambda=$' + str(wave) + 'nm', transform=axarr_SC_ill [row_num, col_num].transAxes, color = 'r')
            axarr_SC_dk [row_num, col_num].text(0.1, 0.9, r'$\lambda=$' + str(wave) + 'nm', transform=axarr_SC_dk [row_num, col_num].transAxes, color = 'r')
            if col_num == 0:
                axarr_PD_ill [row_num, col_num].set_ylabel(r'Photocurrent ($\mu A$)')
                axarr_PD_dk [row_num, col_num].set_ylabel(r'Photocurrent ($\mu A$)')
                axarr_SC_ill [row_num, col_num].set_ylabel(r'Photocurrent ($\mu A$)')
                axarr_SC_dk [row_num, col_num].set_ylabel(r'Photocurrent ($\mu A$)')
            if row_num == 4-1:
                axarr_PD_ill [row_num, col_num].set_xlabel(r'$\Delta t$ (s)')
                axarr_PD_dk [row_num, col_num].set_xlabel(r'$\Delta t$ (s)')
                axarr_SC_ill [row_num, col_num].set_xlabel(r'$\Delta t$ (s)')
                axarr_SC_dk [row_num, col_num].set_xlabel(r'$\Delta t$ (s)')
        f_PD_ill.suptitle('PD illuminated photocurrent - ' + date_str)
        f_PD_dk.suptitle('PD dark photocurrent - ' + date_str)
        f_SC_ill.suptitle('SC illuminated photocurrent - ' + date_str)
        f_SC_dk.suptitle('SC dark photocurrent - ' + date_str)
        f_PD_ill.tight_layout()
        f_PD_dk.tight_layout()
        f_SC_ill.tight_layout()
        f_SC_dk.tight_layout()
        PD_dark_plot_suffix, PD_ill_plot_suffix, SC_dark_plot_suffix, SC_ill_plot_suffix
        f_PD_ill.savefig(data_dir + file_root + PD_ill_plot_suffix)
        f_PD_dk.savefig(data_dir + file_root + PD_dark_plot_suffix)
        f_SC_ill.savefig(data_dir + file_root + SC_ill_plot_suffix)
        f_SC_dk.savefig(data_dir + file_root + SC_dark_plot_suffix)


        arrs = [SC_illum, SC_dark, PD_illum, PD_dark]
        SC_illum_stats, SC_dark_stats, PD_illum_stats, PD_dark_stats = ([[[np.mean(arr[j]) * current_scaling_A_to_uA, np.std(arr[j]) * current_scaling_A_to_uA] for j in range(len(arr))] for arr in arrs])
        SC_illum_mean = [SC_illum_stats[j][0] for j in range(len(SC_illum_stats)) if not (j in bad_data_indeces[i]) ]
        SC_illum_std = [SC_illum_stats[j][1] for j in range(len(SC_illum_stats)) if not (j in bad_data_indeces[i]) ]
        SC_dark_mean = [SC_dark_stats[j][0] for j in range(len(SC_dark_stats)) if not (j in bad_data_indeces[i]) ]
        SC_dark_std = [SC_dark_stats[j][1] for j in range(len(SC_dark_stats)) if not (j in bad_data_indeces[i]) ]
        PD_illum_mean = [PD_illum_stats[j][0] for j in range(len(PD_illum_stats)) if not (j in bad_data_indeces[i]) ]
        PD_illum_std = [PD_illum_stats[j][1] for j in range(len(PD_illum_stats)) if not (j in bad_data_indeces[i]) ]
        PD_dark_mean = [PD_dark_stats[j][0] for j in range(len(PD_dark_stats)) if not (j in bad_data_indeces[i]) ]
        PD_dark_std = [PD_dark_stats[j][1] for j in range(len(PD_dark_stats)) if not (j in bad_data_indeces[i]) ]
        SC_diff = (np.array(SC_illum_mean) - np.array(SC_dark_mean))
        SC_diff_std = np.sqrt(np.array(SC_illum_std) ** 2.0 + np.array(SC_dark_std) ** 2.0)
        PD_diff = (np.array(PD_illum_mean) - np.array(PD_dark_mean))
        PD_diff_std = np.sqrt(np.array(PD_illum_std) ** 2.0 + np.array(PD_dark_std) ** 2.0)

        ratios = SC_diff / PD_diff
        ratio_errs = np.sqrt((SC_diff_std ** 2.0 * PD_diff ** -2.0) + (SC_diff ** 2.0 * PD_diff_std ** 2.0 ) * PD_diff ** -4.0)
        print ('[len(wavelengths), len(SC_diff)] = ' + str([len(wavelengths), len(SC_diff)] ))
        SC_plots = SC_plots + [axarr_QE[0,0].scatter(wavelengths, SC_diff, c=SC_color, marker = '.')]
        axarr_QE[0,0].plot(wavelengths, SC_diff, c=SC_color, marker = '.')
        axarr_QE[0,0].errorbar(wavelengths, SC_diff, yerr = SC_illum_std, color=SC_color, fmt = 'none')
        PD_plots = PD_plots + [axarr_QE[0,0].scatter(wavelengths, PD_diff, c=PD_color, marker = '.')]
        axarr_QE[0,0].plot(wavelengths, PD_diff, c=PD_color, marker = '.')
        axarr_QE[0,0].errorbar(wavelengths, PD_diff, yerr = PD_illum_std, color=PD_color, fmt = 'none')
        ratio_plots = ratio_plots + [axarr_QE[0,1].scatter(wavelengths, ratios, c = combined_color, marker = '.')]
        axarr_QE[0,1].errorbar(wavelengths, ratios, yerr = ratio_errs, fmt = 'none', color = combined_color)
        SC_QE_plots = SC_QE_plots + [axarr_QE[1,1].scatter(wavelengths, PD_QE_interp(wavelengths) * np.array(ratios), c = SC_color, marker = '.')]
        PD_QE_plots = PD_QE_plots + [axarr_QE[1,0].plot(PD_QE_wavelengths, PD_QE_interp(PD_QE_wavelengths), c = PD_color, marker = '.')[0]]
        axarr_QE[1,1].errorbar(wavelengths, PD_QE_interp(wavelengths) * np.array(ratios), yerr = PD_QE_interp(wavelengths) * np.array(ratio_errs), fmt = 'none', color = SC_color)
        axarr_QE[0,0].set_xticks(wavelengths)
        axarr_QE[0,0].set_xticklabels(LEDs, rotation = 60)
        axarr_QE[1,0].set_xticks(wavelengths)
        axarr_QE[1,0].set_xticklabels(LEDs, rotation = 60)
        axarr_QE[0,1].set_xticks(wavelengths)
        axarr_QE[0,1].set_xticklabels(LEDs, rotation = 60)
        axarr_QE[1,1].set_xticks(wavelengths)
        axarr_QE[1,1].set_xticklabels(LEDs, rotation = 60)
        axarr_QE[0,0].set_xlabel(r'Monochromator Wavelength (nm)', fontsize = 12)
        axarr_QE[1,0].set_xlabel(r'Monochromator Wavelength (nm)', fontsize = 12)
        axarr_QE[0,1].set_xlabel(r'Monochromator Wavelength (nm)', fontsize = 12)
        axarr_QE[1,1].set_xlabel(r'Monochromator Wavelength (nm)', fontsize = 12)
        axarr_QE[0,0].set_ylabel(r'Photocurrent ($\mu A$)', fontsize = 12)
        axarr_QE[0,1].set_ylabel(r'SC Photocurrent / PD photocurrent', fontsize = 12)
        axarr_QE[1,0].set_ylabel(r'Quantum efficiency $(e^-/ \gamma)$', fontsize = 12)
        axarr_QE[1,1].set_ylabel(r'Quantum efficiency $(e^-/ \gamma)$', fontsize = 12)
        axarr_QE[0,0].legend([PD_plots[0], SC_plots[0]], ['Photodiode photocurrent', 'Solar cell photocurrent' ])
        axarr_QE[0,1].legend([ratio_plots[0] ], ['Solar Cell / photodiode photocurrent ratio' ])
    #axarr_QE[2].legend([SC_QE_plots[0], PD_QE_plots[0], lp_in_third], ['Inferred Solar Cell QE', 'NIST Photodiode QE', 'Longpass Inserted'])
        axarr_QE[1,0].legend([PD_QE_plots[0]], ['NIST Photodiode QE (NOT OUR DATA)'])
        axarr_QE[1,1].legend([SC_QE_plots[0] ], ['Inferred Solar Cell QE' ])
        [axarr_QE[i % 2 ,i // 2].set_xlim(xlims) for i in range(4)]
        f_QE.suptitle('Solar Cell QE measurements - ' + date_str, fontsize = 16)
        plt.tight_layout()
        f_QE.savefig(data_dir + save_file_names[i])
