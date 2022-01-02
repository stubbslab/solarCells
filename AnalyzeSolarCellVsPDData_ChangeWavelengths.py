import numpy as np
import matplotlib.pyplot as plt
import cantrips as can
import sys
import AstronomicalParameterArchive as apa
import scipy

def readInSCandPDDataFile(data_file, data_dir, SC_neg_current, PD_neg_current, recorded_bias_voltages = 1, ignore_str = 'DELETED'):
    with open (data_dir + data_file, 'r') as open_file:
        data = open_file.readlines()
    data = [row.split('|') for row in data]
    wavelengths = [(row[0])  for row in data]
    PD_illum = [[ ( float(elem)) for elem in can.recursiveStrToListOfLists(row[1]) if elem != ignore_str] for row in data]
    PD_dark = [[(float(elem)) for elem in can.recursiveStrToListOfLists(row[2]) if elem != ignore_str] for row in data]
    SC_illum = [[(float(elem)) for elem in can.recursiveStrToListOfLists(row[3]) if elem != ignore_str] for row in data]
    SC_dark = [[(float(elem)) for elem in can.recursiveStrToListOfLists(row[4]) if elem != ignore_str] for row in data]
    if recorded_bias_voltages:
        PD_illum_bv = [[ (float(elem)) for elem in can.recursiveStrToListOfLists(row[5]) if elem != ignore_str] for row in data]
        PD_dark_bv = [[(float(elem)) for elem in can.recursiveStrToListOfLists(row[6]) if elem != ignore_str] for row in data]
        SC_illum_bv = [[(float(elem)) for elem in can.recursiveStrToListOfLists(row[7]) if elem != ignore_str] for row in data]
        SC_dark_bv = [[(float(elem)) for elem in can.recursiveStrToListOfLists(row[8]) if elem != ignore_str] for row in data]
        SC_temps = [float(row[9]) for row in data]
    else:
        SC_temps = [float(row[5]) for row in data]
    #SC_illum_bv = [[ (float(elem[1:-1])) for elem in can.recursiveStrToListOfLists(row[4])] for row in data]
    #SC_dark_bv = [[(float(elem[1:-1])) for elem in can.recursiveStrToListOfLists(row[5])] for row in data]
    #PD_illum_bv = [[(float(elem[1:-1])) for elem in can.recursiveStrToListOfLists(row[6])] for row in data]
    #PD_dark_bv = [[(float(elem[1:-1])) for elem in can.recursiveStrToListOfLists(row[7])] for row in data]

    #dial_positions, bias_voltages, SC_illum, SC_dark, PD_illum_bv, PD_dark_bv, SC_illum_bv, SC_dark_bv, PD_illum_bv, PD_dark_bv = can.safeSortOneListByAnother(bias_voltages, [zaber_positions, bias_voltages, SC_illum, SC_dark, PD_illum_bv, PD_dark_bv, SC_illum_bv, SC_dark_bv, PD_illum_bv, PD_dark_bv])
    #return [dial_positions, bias_voltages, SC_illum, SC_dark, PD_illum_bv, PD_dark_bv, SC_illum_bv, SC_dark_bv, PD_illum_bv, PD_dark_bv]
    if recorded_bias_voltages:
        wavelengths, SC_illum, SC_dark, PD_illum, PD_dark, SC_illum_bv, SC_dark_bv, PD_illum_bv, PD_dark_bv, SC_temps = can.safeSortOneListByAnother(wavelengths, [wavelengths, SC_illum, SC_dark, PD_illum, PD_dark, SC_illum_bv, SC_dark_bv, PD_illum_bv, PD_dark_bv, SC_temps])
        return [wavelengths, SC_illum, SC_dark, PD_illum, PD_dark, SC_illum_bv, SC_dark_bv, PD_illum_bv, PD_dark_bv, SC_temps]
    else:
        wavelengths, SC_illum, SC_dark, PD_illum, PD_dark = can.safeSortOneListByAnother(wavelengths, [wavelengths, SC_illum, SC_dark, PD_illum, PD_dark, SC_temps])
        return [wavelengths, SC_illum, SC_dark, PD_illum, PD_dark, SC_temps]

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
    dir_root = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/solarCell/qe/'
    #data_files = [ stem + date_str + '.txt' for stem in ['SC_vs_PD_QE_with_B2987A_cellPosA_inten2_10_', 'SC_vs_PD_QE_with_B2987A_cellPosA_inten4_10_', 'SC_vs_PD_QE_with_B2987A_cellPosB_inten2_10_', 'SC_vs_PD_QE_with_B2987A_cellPosB_inten4_10_'] ]

    stems = ['SC_QE_from_mono_SC_ED_' ]
    suffixes = ['_1']
    date_strs = ['20210531']
    extra_title_str = 'Cell ID ED, SS'
    cell_ids = ['Cell ED - ' + date_strs[0] + 'SS']
    shunt_resistances = [610]

    zoomed_in_QE_lims = [0.97, 1.01]
    lp_filter_wavelength = 655
    SC_neg_current = 0
    PD_neg_current = 0
    deleted_element_id = 'DELETED'
    data_files = [ stems[i] + date_strs[i] + suffixes[i] + '.txt' for i in range(len(stems)) ]
    file_roots = [f[0:-len('.txt')] for f in data_files]
    save_file_names = [root + '.pdf' for root in file_roots]
    xlims = [340, 1060]
    recorded_bias_voltages = 1
    title_font_size = 24
    labelsize = 18
    ticksize = 16
    plt.rc('xtick',labelsize=ticksize)
    plt.rc('ytick',labelsize=ticksize)
    time_sep = 0.1

    data_dirs = [dir_root + date_str + '/' for date_str in date_strs]
    save_plot_files = ['SC_QE_measurement_' + date_str + '.pdf' for date_str in date_strs]
    PD_dark_plot_suffix, PD_ill_plot_suffix, SC_dark_plot_suffix, SC_ill_plot_suffix = ['_PD_dark_photocurrents.pdf', '_PD_ill_photocurrents.pdf', '_SC_dark_photocurrents.pdf', '_SC_ill_photocurrents.pdf']
    if recorded_bias_voltages:
        PD_dark_bv_plot_suffix, PD_ill_bv_plot_suffix, SC_dark_bv_plot_suffix, SC_ill_bv_plot_suffix = ['_PD_dark_voltages.pdf', '_PD_ill_voltages.pdf', '_SC_dark_voltages.pdf', '_SC_ill_voltages.pdf']
    PD_QE_data_file = dir_root + 'Hamamatsu_Photodiode_S2281_Spectral_Power_Response.txt'
    PD_QE_data = can.readInColumnsToList(PD_QE_data_file, delimiter = ' ', n_ignore = 1)
    PD_QE_wavelengths, PD_QE_A_per_W = [ [float( wave) for wave in PD_QE_data[0]], [float(qe) for qe in PD_QE_data[1]] ]
    astro_arch = apa.AstronomicalParameterArchive()
    e_charge = astro_arch.getElectronCharge()  #electron charge in coloumbs
    h_JouleS = astro_arch.getPlancksConstant()  #Planck's constant, in Joule X second
    c = astro_arch.getc() #speed of light in km/s
    PD_QE_e_per_phot = [PD_QE_A_per_W[i] / e_charge * h_JouleS * c / (PD_QE_wavelengths[i] * 10.0 ** -12 ) for i in range(len(PD_QE_wavelengths)) ]
    PD_QE_interp =  scipy.interpolate.interp1d(PD_QE_wavelengths, PD_QE_e_per_phot)


    bad_data_indeces = [[ ] for f in data_files]
    current_scaling_v_to_uA = 1
    index_break = []
    file_legend_ids = ['']
    PD_colors, SC_colors, combined_colors, indicator_color = [['b', 'yellow', 'orange', 'pink'],['r', 'cyan', 'brown', 'magenta'], ['purple', 'green', 'grey', 'gold'], 'grey']
    plot_styles = ['-','--']
    SC_plots = []
    PD_plots = []
    PD_QE_plots = []
    full_SC_QE_plots = []
    zoomed_SC_QE_plots = []
    ratio_plots = []
    burn_in_PD = 11
    burn_in_SC = 41

    f_QE, axarr_QE = plt.subplots(2,2, figsize = [15, 15] )
    for i in range(len(data_files)):
        print ('Working on ' + str(i+1) + 'th data file of ' + str(len(data_files)))
        Rsh = shunt_resistances[i]

        data_file = data_files[i]
        file_root = file_roots[i]
        PD_color = PD_colors[i]
        SC_color = SC_colors[i]
        combined_color = combined_colors[i]
        data_dir = data_dirs[i]
        date_str = date_strs[i]

        #z_positions, bias_voltages, SC_illum, SC_dark, PD_illum_bv, PD_dark_bv, SC_illum_bv, SC_dark_bv, PD_illum_bv, PD_dark_bv = readInSCandPDDataFile(data_file, data_dir)
        if recorded_bias_voltages:
            wavelengths, SC_illum, SC_dark, PD_illum, PD_dark, SC_illum_bv, SC_dark_bv, PD_illum_bv, PD_dark_bv, SC_temps = readInSCandPDDataFile(data_file, data_dir, SC_neg_current, PD_neg_current, recorded_bias_voltages = recorded_bias_voltages , ignore_str = deleted_element_id)
            print ('SC_temps = ' + str(SC_temps))
        else:
            wavelengths, SC_illum, SC_dark, PD_illum, PD_dark = readInSCandPDDataFile(data_file, data_dir, SC_neg_current, PD_neg_current, recorded_bias_voltages = recorded_bias_voltages , ignore_str = deleted_element_id )

        full_SC_illum, full_SC_dark = [ [np.array(wave_set) * (-1 if SC_neg_current else 1) * 10.0 ** 6.0 for wave_set in SC_illum] , [np.array(wave_set) * (-1 if SC_neg_current else 1) * 10.0 ** 6.0 for wave_set in SC_dark] ]
        #full_PD_illum, full_PD_dark = [np.array([wave_set for wave_set in PD_illum]) * (-1 if PD_neg_current else 1) * 10.0 ** 6.0, np.array([wave_set for wave_set in PD_dark]) * (-1 if PD_neg_current else 1) * 10.0 ** 6.0]
        full_PD_illum, full_PD_dark = [ [np.array(wave_set) * (-1 if PD_neg_current else 1) * 10.0 ** 6.0 for wave_set in PD_illum], [np.array(wave_set) * (-1 if PD_neg_current else 1) * 10.0 ** 6.0 for wave_set in PD_dark] ]

        if recorded_bias_voltages:
            full_SC_illum_bv, full_SC_dark_bv = [[np.array(wave_set) * (-1 if SC_neg_current else 1) * 10.0 ** 6.0 for wave_set in SC_illum_bv], [np.array(wave_set) * (-1 if SC_neg_current else 1) * 10.0 ** 6.0 for wave_set in SC_dark_bv] ]
            full_PD_illum_bv, full_PD_dark_bv = [[np.array(wave_set) * (-1 if PD_neg_current else 1) * 10.0 ** 6.0 for wave_set in PD_illum_bv], [np.array(wave_set) * (-1 if PD_neg_current else 1) * 10.0 ** 6.0 for wave_set in PD_dark_bv] ]
            SC_ill_ylims, SC_dk_ylims, PD_ill_ylims, PD_dk_ylims, SC_ill_bv_ylims, SC_dk_bv_ylims, PD_ill_bv_ylims, PD_dk_bv_ylims = [[np.median(arr) - 3.0 * np.std(arr), np.median(arr) + 3.0 * np.std(arr)] for arr in [full_SC_illum, full_SC_dark, full_PD_illum, full_PD_dark, full_SC_illum_bv, full_SC_dark_bv, full_PD_illum_bv, full_PD_dark_bv]]
        else:
            SC_ill_ylims, SC_dk_ylims, PD_ill_ylims, PD_dk_ylims = [[np.median(arr) - 3.0 * np.std(arr), np.median(arr) + 3.0 * np.std(arr)] for arr in [full_SC_illum, full_SC_dark, full_PD_illum, full_PD_dark]]

        #We should correct each data point in the measured currents with the measured voltages, since the measured voltages drift a little
        #full_SC_illum = np.array(full_SC_illum_bv) / Rsh + np.array(full_SC_illum)
        #full_SC_dark = np.array(full_SC_dark_bv) / Rsh + np.array(full_SC_dark)

        SC_illum, SC_dark = [[wave_set[burn_in_SC:] for wave_set in arr] for arr in [full_SC_illum, full_SC_dark]]
        PD_illum, PD_dark = [[wave_set[burn_in_PD:] for wave_set in arr] for arr in [full_PD_illum, full_PD_dark]]
        SC_illum_bv, SC_dark_bv = [[wave_set[burn_in_SC:] for wave_set in arr] for arr in [full_SC_illum_bv, full_SC_dark_bv]]
        PD_illum_bv, PD_dark_bv = [[wave_set[burn_in_PD:] for wave_set in arr] for arr in [full_PD_illum_bv, full_PD_dark_bv]]
        n_rows = min(8, len(wavelengths))
        n_single_shot_cols = int(np.ceil(len(wavelengths) / n_rows))
        f_PD_ill, axarr_PD_ill = plt.subplots(n_rows, n_single_shot_cols, figsize = (n_single_shot_cols * 3, n_rows * 3), sharex = True, squeeze = False)
        f_PD_dk, axarr_PD_dk = plt.subplots(n_rows, n_single_shot_cols, figsize = (n_single_shot_cols * 3, n_rows * 3), sharex = True, squeeze = False)
        f_SC_ill, axarr_SC_ill = plt.subplots(n_rows, n_single_shot_cols, figsize = (n_single_shot_cols * 3, n_rows * 3), sharex = True, squeeze = False)
        f_SC_dk, axarr_SC_dk = plt.subplots(n_rows, n_single_shot_cols, figsize = (n_single_shot_cols * 3, n_rows * 3), sharex = True, squeeze = False)
        if recorded_bias_voltages:
            f_PD_bv_ill, axarr_PD_bv_ill = plt.subplots(n_rows, n_single_shot_cols, figsize = (n_single_shot_cols * 3, n_rows * 3), sharex = True, squeeze = False)
            f_PD_bv_dk, axarr_PD_bv_dk = plt.subplots(n_rows, n_single_shot_cols, figsize = (n_single_shot_cols * 3, n_rows * 3), sharex = True, squeeze = False)
            f_SC_bv_ill, axarr_SC_bv_ill = plt.subplots(n_rows, n_single_shot_cols, figsize = (n_single_shot_cols * 3, n_rows * 3), sharex = True, squeeze = False)
            f_SC_bv_dk, axarr_SC_bv_dk = plt.subplots(n_rows, n_single_shot_cols, figsize = (n_single_shot_cols * 3, n_rows * 3), sharex = True, squeeze = False)
        print ('Making plots for SC and PD traces...')
        try:
            wavelength_strs = wavelengths[:]
            true_wavelengths = [float(wave) for wave in wavelengths]
            wavelengths = [float(wave) for wave in wavelengths]
            wavelength_strs = [r'$\lambda=$' + str(wave) + 'nm' for wave in wavelengths]
            xlabel = r'Monochromator Wavelength (nm)'
        except  ValueError:
            wavelength_strs = wavelengths[:]
            true_wavelengths = [625 for i in range(len(wavelengths))]
            wavelengths = [i for i in range(len(wavelengths))]
            wavelength_strs = ['Obs Number ' + str(wave) for wave in wavelengths]
            xlabel = 'Sequence number'
            PD_QE_wavelengths = true_wavelengths
        wavelength_lims = [min(wavelengths), max(wavelengths)]
        true_wavelength_lims = [min(true_wavelengths), max(true_wavelengths)]

        print ('wavelengths = ' + str(wavelengths))
        print ('wavelength_strs = ' + str(wavelength_strs))
        print ('n_single_shot_cols = ' + str(n_single_shot_cols))

        if recorded_bias_voltages:
            all_data_axes = [axarr_PD_ill, axarr_PD_dk, axarr_SC_ill, axarr_SC_dk, axarr_PD_bv_ill, axarr_PD_bv_dk, axarr_SC_bv_ill, axarr_SC_bv_dk, ]
            all_data_figures = [f_PD_ill, f_PD_dk, f_SC_ill, f_SC_dk, f_PD_bv_ill, f_PD_bv_dk, f_SC_bv_ill, f_SC_bv_dk, ]
            all_data_sets = [full_PD_illum, full_PD_dark, full_SC_illum, full_SC_dark, full_PD_illum_bv, full_PD_dark_bv, full_SC_illum_bv, full_SC_dark_bv, ]
            ylims_set = [PD_ill_ylims, PD_dk_ylims, SC_ill_ylims, SC_dk_ylims, PD_ill_bv_ylims, PD_dk_bv_ylims, SC_ill_bv_ylims, SC_dk_bv_ylims]
            full_data_ylabels = [r'Photocurrent ($\mu A$)' for i in range(4)] + [r'Output Voltage ($\mu V$)' for i in range(4)]
            full_data_xlabels = [r'$\Delta t$ (s)' for i in range(8)]
            suptitles = ['PD illuminated photocurrent - ' + date_str, 'PD illuminated photocurrent - ' + date_str, 'SC illuminated photocurrent - ' + date_str, 'SC dark photocurrent - ' + date_str,
                         'PD illuminated output voltage - ' + date_str, 'PD illuminated output voltage - ' + date_str, 'SC illuminated output voltage - ' + date_str, 'SC dark output voltage - ' + date_str]
            plot_suffixes = [ PD_ill_plot_suffix, PD_dark_plot_suffix, SC_ill_plot_suffix, SC_dark_plot_suffix, PD_ill_bv_plot_suffix, PD_dark_bv_plot_suffix, SC_ill_bv_plot_suffix, SC_dark_bv_plot_suffix, ]
        else:
            all_data_axes = [axarr_PD_ill, axarr_PD_dk, axarr_SC_ill, axarr_SC_dk, ]
            all_data_figures = [f_PD_ill, f_PD_dk, f_SC_ill, f_SC_dk, ]
            all_data_sets = [full_PD_illum, full_PD_dark, full_SC_illum, full_SC_dark]
            ylims_set = [PD_ill_ylims, PD_dk_ylims, SC_ill_ylims, SC_dk_ylims]
            full_data_ylabels = [r'Photocurrent ($\mu A$)' for i in range(4)]
            full_data_xlabels = [r'$\Delta t$ (s)' for i in range(4)]
            suptitles = ['PD illuminated photocurrent - ' + date_str, 'PD illuminated photocurrent - ' + date_str, 'SC illuminated photocurrent - ' + date_str, 'SC dark photocurrent - ' + date_str]
            plot_suffixes = [ PD_ill_plot_suffix, PD_dark_plot_suffix, SC_ill_plot_suffix, SC_dark_plot_suffix, ]
        for j in range(len(wavelengths)):
            wave = wavelengths[j]
            wavelength_str = wavelength_strs[j]
            row_num, col_num = [int(j // n_single_shot_cols), int(j % n_single_shot_cols)]
            for k in range(len(all_data_axes)):
                data_axis_arr = all_data_axes [k]
                all_data_set = all_data_sets[k]
                ylims = ylims_set [k]
                full_data_ylabel = full_data_ylabels[k]
                full_data_xlabel = full_data_xlabels[k]
                data_axis_arr [row_num, col_num].plot([time_sep * t for t in range(len(all_data_set[j]))], all_data_set[j], color = 'k')
                data_axis_arr[row_num, col_num].axvline(time_sep * (burn_in_SC - 0.5), linestyle = '--', c = 'grey', alpha = 0.5)
                data_axis_arr[row_num, col_num].set_ylim(ylims)
                data_axis_arr [row_num, col_num].text(0.1, 0.9, wavelength_str, transform=data_axis_arr [row_num, col_num].transAxes, color = 'r')
                if col_num == 0:
                    data_axis_arr [row_num, col_num].set_ylabel(full_data_ylabel, fontsize = labelsize)
                if row_num == n_rows-1:
                    data_axis_arr [row_num, col_num].set_xlabel(full_data_xlabel, fontsize = labelsize)
        print ('Done making plots for SC and PD traces. ')

        for k in range(len(all_data_axes)):
            single_data_type_axes = all_data_axes[k]
            single_data_figure =all_data_figures [k]
            single_data_figure.suptitle(suptitles[k])
            single_data_figure.tight_layout()
            single_data_figure.savefig(data_dir + file_root + plot_suffixes[k])

        arrs = [SC_illum, SC_dark, PD_illum, PD_dark]
        SC_illum_stats, SC_dark_stats, PD_illum_stats, PD_dark_stats = ([[[np.mean(arr[j]) * current_scaling_v_to_uA, np.std(arr[j]) / np.sqrt(len(arr[j])) * current_scaling_v_to_uA] for j in range(len(arr))] for arr in arrs])
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

        """
        print ('Working on QE plots...')
        print ('[SC_illum_mean[4], SC_illum_mean[5]] = ' + str([SC_illum_mean[4], SC_illum_mean[5]]))
        print ('[SC_illum_std[4], SC_illum_std[5]] = ' + str([SC_illum_std[4], SC_illum_std[5]]))
        print ('(SC_illum_mean[4] - SC_illum_mean[5]) / np.sqrt(SC_illum_std[4] ** 2.0 + SC_illum_std[5] ** 2.0) = ' + str((SC_illum_mean[4] - SC_illum_mean[5]) / np.sqrt(SC_illum_std[4] ** 2.0 + SC_illum_std[5] ** 2.0)))
        print ('[PD_illum_mean[4], PD_illum_mean[5]] = ' + str([PD_illum_mean[4], PD_illum_mean[5]]))
        print ('[PD_illum_std[4], PD_illum_std[5]] = ' + str([PD_illum_std[4], PD_illum_std[5]]))
        print ('(PD_illum_mean[4] - PD_illum_mean[5]) / np.sqrt(PD_illum_std[4] ** 2.0 + PD_illum_std[5] ** 2.0) = ' + str((PD_illum_mean[4] - PD_illum_mean[5]) / np.sqrt(PD_illum_std[4] ** 2.0 + PD_illum_std[5] ** 2.0)))
        print ('[SC_dark_mean[4], SC_dark_mean[5]] = ' + str([SC_dark_mean[4], SC_dark_mean[5]]))
        print ('[SC_dark_std[4], SC_dark_std[5]] = ' + str([SC_dark_std[4], SC_dark_std[5]]))
        print ('(SC_dark_mean[4] - SC_dark_mean[5]) / np.sqrt(SC_dark_std[4] ** 2.0 + SC_dark_std[5] ** 2.0) = ' + str((SC_dark_mean[4] - SC_dark_mean[5]) / np.sqrt(SC_dark_std[4] ** 2.0 + SC_dark_std[5] ** 2.0)))
        print ('[PD_dark_mean[4], PD_dark_mean[5]] = ' + str([PD_dark_mean[4], PD_dark_mean[5]]))
        print ('[PD_dark_std[4], PD_dark_std[5]] = ' + str([PD_dark_std[4], PD_dark_std[5]]))
        print ('(PD_dark_mean[4] - PD_dark_mean[5]) / np.sqrt(PD_dark_std[4] ** 2.0 + PD_dark_std[5] ** 2.0) = ' + str((PD_dark_mean[4] - PD_dark_mean[5]) / np.sqrt(PD_dark_std[4] ** 2.0 + PD_dark_std[5] ** 2.0)))
        print ('[SC_diff[4], SC_diff[5]] = ' + str([SC_diff[4], SC_diff[5]]))
        print ('[PD_diff[4], PD_diff[5]] = ' + str([PD_diff[4], PD_diff[5]]))
        print ('[ratios[4], ratios[5]] = ' + str([ratios[4], ratios[5]] ))
        print ('[len(wavelengths), len(SC_diff)] = ' + str([len(wavelengths), len(SC_diff)] ))
        """

        SC_plots = SC_plots + [axarr_QE[0,0].scatter(wavelengths, SC_diff, c=SC_color, marker = '.')]
        axarr_QE[0,0].plot(wavelengths, SC_diff, c=SC_color, marker = '.')
        axarr_QE[0,0].errorbar(wavelengths, SC_diff, yerr = SC_illum_std, color=SC_color, fmt = 'none')
        PD_plots = PD_plots + [axarr_QE[0,0].scatter(wavelengths, PD_diff, c=PD_color, marker = '.')]
        axarr_QE[0,0].plot(wavelengths, PD_diff, c=PD_color, marker = '.')
        axarr_QE[0,0].errorbar(wavelengths, PD_diff, yerr = PD_illum_std, color=PD_color, fmt = 'none')
        lp_in_first = axarr_QE[0,0].axvline(lp_filter_wavelength, color = indicator_color, alpha = 0.5, linestyle = 'dashed')
        ratio_plots = ratio_plots + [axarr_QE[0,1].scatter(wavelengths, ratios, c = combined_color, marker = '.')]
        axarr_QE[0,1].errorbar(wavelengths, ratios, yerr = ratio_errs, fmt = 'none', color = combined_color)
        lp_in_second = axarr_QE[0,1].axvline(lp_filter_wavelength, color = indicator_color, alpha = 0.5, linestyle = 'dashed')
        PD_QE_plots = PD_QE_plots + [axarr_QE[1,0].plot(wavelengths, PD_QE_interp(wavelengths), c = PD_color, marker = '.')[0]]
        full_SC_QE_plots = full_SC_QE_plots + [axarr_QE[1,0].scatter(wavelengths, PD_QE_interp(true_wavelengths) * np.array(ratios), c = SC_color, marker = '.')]
        axarr_QE[1,1].axhline(1.0, 0, 1, color = 'k', linestyle = 'dashed')
        zoomed_SC_QE_plots = zoomed_SC_QE_plots + [axarr_QE[1,1].scatter(wavelengths, PD_QE_interp(true_wavelengths) * np.array(ratios), c = SC_color, marker = '.')]
        axarr_QE[1,1].errorbar(wavelengths, PD_QE_interp(true_wavelengths) * np.array(ratios), yerr = PD_QE_interp(true_wavelengths) * np.array(ratio_errs), fmt = 'none', color = SC_color)

        axarr_QE[0,0].set_xlim(wavelength_lims)
        axarr_QE[0,1].set_xlim(wavelength_lims)
        axarr_QE[1,0].set_xlim(wavelength_lims)
        axarr_QE[1,1].set_xlim(wavelength_lims)

        lp_in_forth = axarr_QE[1,1].axvline(lp_filter_wavelength, 0, 1, color = indicator_color, alpha = 0.5, linestyle = 'dashed')
        axarr_QE[0,0].set_xlabel(xlabel, fontsize = labelsize)
        axarr_QE[1,0].set_xlabel(xlabel, fontsize = labelsize)
        axarr_QE[0,1].set_xlabel(xlabel, fontsize = labelsize)
        axarr_QE[1,1].set_xlabel(xlabel, fontsize = labelsize)
        axarr_QE[0,0].set_ylabel(r'Photocurrent ($\mu A$)', fontsize = labelsize)
        axarr_QE[0,1].set_ylabel(r'SC Photocurrent / PD photocurrent', fontsize = labelsize)
        axarr_QE[1,0].set_ylabel(r'Quantum efficiency $(e^-/ \gamma)$', fontsize = labelsize)
        axarr_QE[1,1].set_ylabel(r'Quantum efficiency $(e^-/ \gamma)$', fontsize = labelsize)

        axarr_QE[1,1].set_ylim(zoomed_in_QE_lims)
        #[axarr_QE[i % 2 ,i // 2].set_xlim(xlims) for i in range(4)]
    axarr_QE[0,0].legend(PD_plots + SC_plots + [lp_in_first], ['Photodiode photocurrent - ' + cell_id for cell_id in cell_ids] + ['Solar cell photocurrent - ' + cell_id for cell_id in cell_ids] + ['>600nm Longpass Inserted'] )
    axarr_QE[0,1].legend(ratio_plots + [lp_in_second], ['Solar Cell / photodiode photocurrent ratio - ' + cell_id for cell_id in cell_ids] + ['>600nm Longpass Inserted'])
#axarr_QE[2].legend([SC_QE_plots[0], PD_QE_plots[0], lp_in_third], ['Inferred Solar Cell QE', 'NIST Photodiode QE', 'Longpass Inserted'])
    axarr_QE[1,0].legend([PD_QE_plots[-1]] + full_SC_QE_plots, ['NIST Photodiode QE (NOT our data)'] + ['Inferred Solar Cell QE - ' + cell_id for cell_id in cell_ids])
    axarr_QE[1,1].legend(zoomed_SC_QE_plots + [lp_in_forth], ['Inferred Solar Cell QE - ' + cell_id for cell_id in cell_ids] + ['>600nm Longpass Inserted'])
    f_QE.suptitle('Solar Cell QE measurements - ' + extra_title_str , fontsize = title_font_size)
    plt.tight_layout()
    print ('Saving QE plots...')
    f_QE.savefig(data_dir + save_file_names[-1])
    print ('Done.')
