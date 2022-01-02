import numpy as np
import matplotlib.pyplot as plt
import cantrips as can
import sys
import AstronomicalParameterArchive as apa
import scipy

def readInSCandPDDataFile(data_file, data_dir, SC_neg_current, PD_neg_current, ignore_str = 'DELETED'):
    with open (data_dir + data_file, 'r') as open_file:
        data = open_file.readlines()
    data = [row.split('|') for row in data]
    stage_pos = [float(row[0])  for row in data]
    PD_ill_times = [[ ( float(elem)) for elem in can.recursiveStrToListOfLists(row[1]) if elem != ignore_str] for row in data]
    PD_illum = [[ ( float(elem)) for elem in can.recursiveStrToListOfLists(row[2]) if elem != ignore_str] for row in data]
    PD_illum_bv = [[(float(elem)) for elem in can.recursiveStrToListOfLists(row[3]) if elem != ignore_str] for row in data]
    PD_dk_times = [[ ( float(elem)) for elem in can.recursiveStrToListOfLists(row[4]) if elem != ignore_str] for row in data]
    PD_dark = [[(float(elem)) for elem in can.recursiveStrToListOfLists(row[5]) if elem != ignore_str] for row in data]
    PD_dark_bv = [[(float(elem)) for elem in can.recursiveStrToListOfLists(row[6]) if elem != ignore_str] for row in data]
    SC_ill_times = [[ (float(elem)) for elem in can.recursiveStrToListOfLists(row[7]) if elem != ignore_str] for row in data]
    SC_illum = [[(float(elem)) for elem in can.recursiveStrToListOfLists(row[8]) if elem != ignore_str] for row in data]
    SC_illum_bv = [[(float(elem)) for elem in can.recursiveStrToListOfLists(row[9]) if elem != ignore_str] for row in data]
    SC_dk_times = [[ (float(elem)) for elem in can.recursiveStrToListOfLists(row[10]) if elem != ignore_str] for row in data]
    SC_dark = [[(float(elem)) for elem in can.recursiveStrToListOfLists(row[11]) if elem != ignore_str] for row in data]
    SC_dark_bv = [[(float(elem)) for elem in can.recursiveStrToListOfLists(row[12]) if elem != ignore_str] for row in data]
    #SC_dark_bv = [[(float(elem)) for elem in can.recursiveStrToListOfLists(row[7]) if elem != ignore_str] for row in data]
    SC_temps = [float(row[13]) for row in data]
    #SC_illum_bv = [[ (float(elem[1:-1])) for elem in can.recursiveStrToListOfLists(row[4])] for row in data]
    #SC_dark_bv = [[(float(elem[1:-1])) for elem in can.recursiveStrToListOfLists(row[5])] for row in data]
    #PD_illum_bv = [[(float(elem[1:-1])) for elem in can.recursiveStrToListOfLists(row[6])] for row in data]
    #PD_dark_bv = [[(float(elem[1:-1])) for elem in can.recursiveStrToListOfLists(row[7])] for row in data]

    #dial_positions, bias_voltages, SC_illum, SC_dark, PD_illum_bv, PD_dark_bv, SC_illum_bv, SC_dark_bv, PD_illum_bv, PD_dark_bv = can.safeSortOneListByAnother(bias_voltages, [zaber_positions, bias_voltages, SC_illum, SC_dark, PD_illum_bv, PD_dark_bv, SC_illum_bv, SC_dark_bv, PD_illum_bv, PD_dark_bv])
    #return [dial_positions, bias_voltages, SC_illum, SC_dark, PD_illum_bv, PD_dark_bv, SC_illum_bv, SC_dark_bv, PD_illum_bv, PD_dark_bv]
    #wavelengths, stage_pos, SC_ill_times, SC_illum, SC_illum_bv, SC_dk_times, SC_dark, SC_dark_bv, PD_ill_times, PD_illum, PD_illum_bv, PD_dk_times, PD_dark, PD_dark_bv, SC_temps = can.safeSortOneListByAnother(wavelengths, [wavelengths, stage_pos, SC_ill_times, SC_illum, SC_illum_bv, SC_dk_times, SC_dark, SC_dark_bv, PD_ill_times, PD_illum, PD_illum_bv, PD_dk_times, PD_dark, PD_dark_bv, SC_temps ])
    return [stage_pos, SC_ill_times, SC_illum, SC_illum_bv, SC_dk_times, SC_dark, SC_dark_bv, PD_ill_times, PD_illum, PD_illum_bv, PD_dk_times, PD_dark, PD_dark_bv, SC_temps]


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

    stems = ['SC_QE_moved_stage_SC_ED_']
    #stems = ['SC_QE_from_mono_SC_ED_']
    suffixes = ['_singleStream1_straightSource']
    #suffixes = ['_singleStream1_straightSource']
    date_strs = ['20210611']
    #date_strs = ['20210607']
    extra_title_str = 'Cell ID ED, Move SC '
    cell_ids = ['Cell ED - ' + date_strs[0] + ' 1']
    #cell_ids = ['Cell ED - ' + date_strs[0] + ' AM']
    shunt_resistances = [610 ]
    #shunt_resistances = [610]
    plot_wavelengths = 0
    #if we don't plot wavelengths, we plot in times

    zoomed_in_QE_lims = [0.96, 1.005]
    lp_filter_wavelength = 655
    insert_longpass = 0
    SC_neg_current = 0
    PD_neg_current = 0
    deleted_element_id = 'DELETED'
    data_files = [ stems[i] + date_strs[i] + suffixes[i] + '.txt' for i in range(len(stems)) ]
    file_roots = [f[0:-len('.txt')] for f in data_files]
    save_file_names = [root + '.pdf' for root in file_roots]
    xlims = [340, 1060]
    title_font_size = 24
    labelsize = 18
    ticksize = 16
    plt.rc('xtick',labelsize=ticksize)
    plt.rc('ytick',labelsize=ticksize)
    time_sep = 0.1

    data_dirs = [dir_root + date_str + '/' for date_str in date_strs]
    save_plot_files = ['SC_QE_measurement_' + date_str + '.pdf' for date_str in date_strs]
    PD_dark_plot_suffix, PD_ill_plot_suffix, SC_dark_plot_suffix, SC_ill_plot_suffix = ['_PD_dark_photocurrents.pdf', '_PD_ill_photocurrents.pdf', '_SC_dark_photocurrents.pdf', '_SC_ill_photocurrents.pdf']
    #PD_ill_plot_suffix, SC_ill_plot_suffix = ['_PD_ill_photocurrents.pdf', '_SC_ill_photocurrents.pdf']
    PD_dark_bv_plot_suffix, PD_ill_bv_plot_suffix, SC_dark_bv_plot_suffix, SC_ill_bv_plot_suffix = ['_PD_dark_voltages.pdf', '_PD_ill_voltages.pdf', '_SC_dark_voltages.pdf', '_SC_ill_voltages.pdf']
    #PD_ill_bv_plot_suffix, SC_ill_bv_plot_suffix = ['_PD_ill_voltages.pdf', '_SC_ill_voltages.pdf']
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
    PD_colors, SC_colors, combined_colors, indicator_color = [['b', 'yellow', 'orange', 'pink', 'purple', 'grey'], ['r', 'cyan', 'brown', 'magenta', 'green', 'gold'], ['purple', 'green', 'grey', 'gold', 'r', 'cyan'], 'grey']
    plot_styles = ['-','--']
    SC_plots = []
    PD_plots = []
    PD_QE_plots = []
    full_SC_QE_plots = []
    zoomed_SC_QE_plots = []
    ratio_plots = []
    burn_in_PD = 11
    burn_in_SC = 51
    x_lims = [np.inf, 0.0] #Initialized here - this will be filled in by the code.

    f_QE, axarr_QE = plt.subplots(2,2, figsize = [15, 15],  )
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

        stage_positions, SC_illum_times, SC_illum, SC_illum_bv, SC_dark_times, SC_dark, SC_dark_bv, PD_illum_times, PD_illum, PD_illum_bv, PD_dark_times, PD_dark, PD_dark_bv, SC_temps = readInSCandPDDataFile(data_file, data_dir, SC_neg_current, PD_neg_current, ignore_str = deleted_element_id)
        print ('SC_temps = ' + str(SC_temps))
        print ('len(SC_temps) = ' + str(len(SC_temps)))

        full_SC_illum_times, full_PD_illum_times, full_SC_dark_times, full_PD_dark_times = [[np.array(wave_set) for wave_set in arr] for arr in [SC_illum_times, PD_illum_times, SC_dark_times, PD_dark_times ]]
        full_SC_illum, full_SC_dark =  [[np.array(wave_set) * (-1 if SC_neg_current else 1) * 10.0 ** 6.0 for wave_set in arr] for arr in [SC_illum, SC_dark]]
        #full_PD_illum, full_PD_dark = [np.array([wave_set for wave_set in PD_illum]) * (-1 if PD_neg_current else 1) * 10.0 ** 6.0, np.array([wave_set for wave_set in PD_dark]) * (-1 if PD_neg_current else 1) * 10.0 ** 6.0]
        full_PD_illum, full_PD_dark = [[np.array(wave_set) * (-1 if PD_neg_current else 1) * 10.0 ** 6.0 for wave_set in arr] for arr in [PD_illum, PD_dark]]
        full_SC_illum_bv, full_SC_dark_bv = [[np.array(wave_set) * (-1 if SC_neg_current else 1) * 10.0 ** 6.0 for wave_set in arr] for arr in [SC_illum_bv, SC_dark_bv]]
        full_PD_illum_bv, full_PD_dark_bv = [[np.array(wave_set) * (-1 if PD_neg_current else 1) * 10.0 ** 6.0 for wave_set in arr] for arr in [PD_illum_bv, PD_dark_bv]]
        arrs_to_display = [full_SC_illum, full_PD_illum, full_SC_dark, full_PD_dark, full_SC_illum_bv, full_PD_illum_bv, full_SC_dark_bv, full_PD_dark_bv ]
        clipped_sigs_of_data_type = [0.0 for arr in arrs_to_display ]
        for j in range(len(arrs_to_display)):
            arr = arrs_to_display[j]
            clipped_sigs_of_data_type[j] = np.std(arr)
        print ('clipped_sigs_of_data_type = ' + str(clipped_sigs_of_data_type))

        SC_ill_ylims, PD_ill_ylims, SC_dk_ylims, PD_dk_ylims, SC_ill_bv_ylims, PD_ill_bv_ylims, SC_dk_bv_ylims, PD_dk_bv_ylims = [[max(np.median(arrs_to_display[j]) - clipped_sigs_of_data_type[j] * 3.0, np.min(arrs_to_display[j])), min(np.median(arrs_to_display[j]) + clipped_sigs_of_data_type[j] * 3.0, np.max(arrs_to_display[j]))] for j in range(len(arrs_to_display)) ]
        print ('[SC_ill_ylims, PD_ill_ylims, SC_dk_ylims, PD_dk_ylims, SC_ill_bv_ylims, PD_ill_bv_ylims, SC_dk_bv_ylims, PD_dk_bv_ylims] = ' + str([SC_ill_ylims, PD_ill_ylims, SC_dk_ylims, PD_dk_ylims, SC_ill_bv_ylims, PD_ill_bv_ylims, SC_dk_bv_ylims, PD_dk_bv_ylims] ))

        #We should correct each data point in the measured currents with the measured voltages, since the measured voltages drift a little
        #full_SC_illum = np.array(full_SC_illum_bv) / Rsh + np.array(full_SC_illum)
        #full_SC_dark = np.array(full_SC_dark_bv) / Rsh + np.array(full_SC_dark)

        print ('Here 0')
        SC_illum_times, SC_illum, SC_illum_bv, SC_dark_times, SC_dark, SC_dark_bv = [[wave_set[burn_in_SC:] for wave_set in arr] for arr in [full_SC_illum_times, full_SC_illum, full_SC_illum_bv, full_SC_dark_times, full_SC_dark, full_SC_dark_bv]]
        print ('Here 01')
        PD_illum_times, PD_illum, PD_illum_bv, PD_dark_times, PD_dark, PD_dark_bv = [[wave_set[burn_in_PD:] for wave_set in arr] for arr in [full_PD_illum_times, full_PD_illum, full_PD_illum_bv, full_PD_dark_times, full_PD_dark, full_PD_dark_bv]]

        print ('Here 1')
        n_rows = min(16, len(stage_positions))
        n_single_shot_cols = int(np.ceil(len(stage_positions) / n_rows))
        print ('Here 2')
        f_PD_ill, axarr_PD_ill = plt.subplots(n_rows, n_single_shot_cols, figsize = (n_single_shot_cols * 3, n_rows * 3), sharex = True, squeeze = False)
        print ('Here 3')
        f_PD_dk, axarr_PD_dk = plt.subplots(n_rows, n_single_shot_cols, figsize = (n_single_shot_cols * 3, n_rows * 3), sharex = True, squeeze = False)
        print ('Here 4')
        f_SC_ill, axarr_SC_ill = plt.subplots(n_rows, n_single_shot_cols, figsize = (n_single_shot_cols * 3, n_rows * 3), sharex = True, squeeze = False)
        print ('Here 5')
        f_SC_dk, axarr_SC_dk = plt.subplots(n_rows, n_single_shot_cols, figsize = (n_single_shot_cols * 3, n_rows * 3), sharex = True, squeeze = False)
        print ('Here 6')
        f_PD_bv_ill, axarr_PD_bv_ill = plt.subplots(n_rows, n_single_shot_cols, figsize = (n_single_shot_cols * 3, n_rows * 3), sharex = True, squeeze = False)
        print ('Here 7')
        f_PD_bv_dk, axarr_PD_bv_dk = plt.subplots(n_rows, n_single_shot_cols, figsize = (n_single_shot_cols * 3, n_rows * 3), sharex = True, squeeze = False)
        print ('Here 8')
        f_SC_bv_ill, axarr_SC_bv_ill = plt.subplots(n_rows, n_single_shot_cols, figsize = (n_single_shot_cols * 3, n_rows * 3), sharex = True, squeeze = False)
        print ('Here 9')
        f_SC_bv_dk, axarr_SC_bv_dk = plt.subplots(n_rows, n_single_shot_cols, figsize = (n_single_shot_cols * 3, n_rows * 3), sharex = True, squeeze = False)
        print ('Making plots for SC and PD traces...')

        true_wavelengths = [650 for i in range(len(stage_positions))]
        data_pos_strs = [str(can.round_to_n(pos, 4)) + 'mm' for pos in stage_positions]
        SC_x_vals = stage_positions[:]
        PD_x_vals = stage_positions[:]
        plot_strs = ['Obs Number ' + str(i) for i in range(len(stage_positions))]
        xlabel = r'SC Stage position (mm)'
        local_x_min = min(SC_x_vals + PD_x_vals) - (PD_x_vals[1] - PD_x_vals[0]) / 2.0
        local_x_max = max(SC_x_vals + PD_x_vals) + (PD_x_vals[1] - PD_x_vals[0]) / 2.0
        x_lims = [min(local_x_min, x_lims[0]), max(local_x_max, x_lims[1])]


        all_data_axes = [axarr_PD_ill, axarr_SC_ill, axarr_PD_dk, axarr_SC_dk, axarr_PD_bv_ill, axarr_SC_bv_ill, axarr_PD_bv_dk, axarr_SC_bv_dk ]
        all_data_figures = [f_PD_ill, f_SC_ill, f_PD_dk, f_SC_dk, f_PD_bv_ill, f_SC_bv_ill, f_PD_bv_dk, f_SC_bv_dk]
        all_data_x_sets = [[np.array(time_set) - time_set[0] for time_set in full_PD_illum_times], [np.array(time_set) - time_set[0] for time_set in full_SC_illum_times],
                           [np.array(time_set) - time_set[0] for time_set in full_PD_dark_times], [np.array(time_set) - time_set[0] for time_set in full_SC_dark_times],
                            [np.array(time_set) - time_set[0] for time_set in full_PD_illum_times], [np.array(time_set) - time_set[0] for time_set in full_SC_illum_times ],
                             [np.array(time_set) - time_set[0] for time_set in full_PD_dark_times], [np.array(time_set) - time_set[0] for time_set in full_SC_dark_times ] ]
        print ('full_PD_illum_times[0] = ' + str(full_PD_illum_times[0] ))
        all_data_y_sets = [full_PD_illum, full_SC_illum, full_PD_dark, full_SC_dark, full_PD_illum_bv, full_SC_illum_bv, full_PD_dark_bv, full_SC_dark_bv ]
        ylims_set = [PD_ill_ylims, SC_ill_ylims,  PD_dk_ylims, SC_dk_ylims,  PD_ill_bv_ylims,  SC_ill_bv_ylims,  PD_dk_bv_ylims,  SC_dk_bv_ylims ]
        full_data_ylabels = [r'Photocurrent ($\mu A$)' for i in range(4)] + [r'Output Voltage ($\mu V$)' for i in range(4)]
        full_data_xlabels = [r'SC stage position (mm)' for i in range(8)]
        suptitles = ['PD illuminated photocurrent - ' + date_str, 'SC illuminated photocurrent - ' + date_str, 'PD dark photocurrent - ' + date_str, 'SC dark photocurrent - ' + date_str,
                     'PD illuminated output voltage - ' + date_str, 'SC illuminated output voltage - ' + date_str, 'PD dark output voltage - ' + date_str, 'SC dark output voltage - ' + date_str]
        #suptitles = ['PD illuminated photocurrent - ' + date_str, 'SC illuminated photocurrent - ' + date_str,
        #             'PD illuminated output voltage - ' + date_str, 'SC illuminated output voltage - ' + date_str ]
        plot_suffixes = [ PD_ill_plot_suffix, SC_ill_plot_suffix, PD_dark_plot_suffix, SC_dark_plot_suffix, PD_ill_bv_plot_suffix, SC_ill_bv_plot_suffix, PD_dark_bv_plot_suffix, SC_dark_bv_plot_suffix]
        v_line_indeces = [burn_in_PD, burn_in_SC, burn_in_PD, burn_in_SC, burn_in_PD, burn_in_SC, burn_in_PD, burn_in_SC]

        for j in range(len(stage_positions)):
            print ('Working on j = ' + str(j) + ' of ' + str(len(stage_positions)) )
            row_num, col_num = [int(j // n_single_shot_cols), int(j % n_single_shot_cols)]
            for k in range(len(all_data_axes)):
                data_pos_str = data_pos_strs[k]
                data_axis_arr = all_data_axes [k]
                all_data_x_set = all_data_x_sets[k]
                all_data_y_set = all_data_y_sets[k]
                ylims = ylims_set [k]
                full_data_ylabel = full_data_ylabels[k]
                full_data_xlabel = full_data_xlabels[k]
                data_axis_arr [row_num, col_num].plot(all_data_x_set[j], all_data_y_set[j], color = 'k')
                data_axis_arr[row_num, col_num].axvline(all_data_x_set[j][v_line_indeces[k]], linestyle = '--', c = 'grey', alpha = 0.5)
                data_axis_arr[row_num, col_num].set_ylim(ylims)
                data_axis_arr [row_num, col_num].text(0.1, 0.9, data_pos_str, transform=data_axis_arr [row_num, col_num].transAxes, color = 'r')
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
        #SC_illum_stats, PD_illum_stats,  = ([[[np.mean(arr[j]) * current_scaling_v_to_uA, np.std(arr[j]) / np.sqrt(len(arr[j])) * current_scaling_v_to_uA] for j in range(len(arr))] for arr in arrs])
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
        #SC_diff_std = np.array(SC_illum_std)
        PD_diff = (np.array(PD_illum_mean) - np.array(PD_dark_mean))
        PD_diff_std = np.sqrt(np.array(PD_illum_std) ** 2.0 + np.array(PD_dark_std) ** 2.0)
        #PD_diff_std = np.array(PD_illum_std)
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

        SC_plots = SC_plots + [axarr_QE[0,0].scatter(SC_x_vals, SC_diff, c=SC_color, marker = '.')]
        #axarr_QE[0,0].plot(SC_x_vals, SC_diff, c=SC_color, marker = '.')
        axarr_QE[0,0].errorbar(SC_x_vals, SC_diff, yerr = SC_illum_std, color=SC_color, fmt = 'none')
        PD_plots = PD_plots + [axarr_QE[0,0].scatter(PD_x_vals, PD_diff, c=PD_color, marker = '.')]
        #axarr_QE[0,0].plot(PD_x_vals, PD_diff, c=PD_color, marker = '.')
        axarr_QE[0,0].errorbar(PD_x_vals, PD_diff, yerr = PD_illum_std, color=PD_color, fmt = 'none')
        #lp_in_first = axarr_QE[0,0].axvline(lp_filter_wavelength, color = indicator_color, alpha = 0.5, linestyle = 'dashed')
        ratio_plots = ratio_plots + [axarr_QE[0,1].scatter(SC_x_vals, ratios, c = combined_color, marker = '.')]
        axarr_QE[0,1].errorbar(SC_x_vals, ratios, yerr = ratio_errs, fmt = 'none', color = combined_color)
        #lp_in_second = axarr_QE[0,1].axvline(lp_filter_wavelength, color = indicator_color, alpha = 0.5, linestyle = 'dashed')
        PD_QE_plots = PD_QE_plots + [axarr_QE[1,0].plot(SC_x_vals, PD_QE_interp(true_wavelengths), c = PD_color, marker = '.')[0]]
        full_SC_QE_plots = full_SC_QE_plots + [axarr_QE[1,0].scatter(SC_x_vals, PD_QE_interp(true_wavelengths) * np.array(ratios), c = SC_color, marker = '.')]
        axarr_QE[1,1].axhline(1.0, 0, 1, color = 'k', linestyle = 'dashed')
        zoomed_SC_QE_plots = zoomed_SC_QE_plots + [axarr_QE[1,1].scatter(SC_x_vals, PD_QE_interp(true_wavelengths) * np.array(ratios), c = SC_color, marker = '.')]
        axarr_QE[1,1].errorbar(SC_x_vals, PD_QE_interp(true_wavelengths) * np.array(ratios), yerr = PD_QE_interp(true_wavelengths) * np.array(ratio_errs), fmt = 'none', color = SC_color)

        axarr_QE[0,0].set_xlim(x_lims)
        axarr_QE[0,1].set_xlim(x_lims)
        axarr_QE[1,0].set_xlim(x_lims)
        axarr_QE[1,1].set_xlim(x_lims)

        #lp_in_forth = axarr_QE[1,1].axvline(lp_filter_wavelength, 0, 1, color = indicator_color, alpha = 0.5, linestyle = 'dashed')
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
    if insert_longpass:
        axarr_QE[0,0].legend(PD_plots + SC_plots + [lp_in_first], ['Photodiode photocurrent - ' + cell_id for cell_id in cell_ids] + ['Solar cell photocurrent - ' + cell_id for cell_id in cell_ids] + ['>600nm Longpass Inserted'] )
        axarr_QE[0,1].legend(ratio_plots + [lp_in_second], ['Solar Cell / photodiode photocurrent ratio - ' + cell_id for cell_id in cell_ids] + ['>600nm Longpass Inserted'])
#axarr_QE[2].legend([SC_QE_plots[0], PD_QE_plots[0], lp_in_third], ['Inferred Solar Cell QE', 'NIST Photodiode QE', 'Longpass Inserted'])
        axarr_QE[1,0].legend([PD_QE_plots[-1]] + full_SC_QE_plots, ['NIST Photodiode QE (NOT our data)'] + ['Inferred Solar Cell QE - ' + cell_id for cell_id in cell_ids])
        axarr_QE[1,1].legend(zoomed_SC_QE_plots + [lp_in_forth], ['Inferred Solar Cell QE - ' + cell_id for cell_id in cell_ids] + ['>600nm Longpass Inserted'])
    else:
        axarr_QE[0,0].legend(SC_plots + PD_plots, ['Solar cell photocurrent - ' + cell_id for cell_id in cell_ids] + ['Photodiode photocurrent - ' + cell_id for cell_id in cell_ids]  )
        axarr_QE[0,1].legend(ratio_plots  , ['Solar Cell / photodiode photocurrent ratio - ' + cell_id for cell_id in cell_ids] )
#axarr_QE[2].legend([SC_QE_plots[0], PD_QE_plots[0], lp_in_third], ['Inferred Solar Cell QE', 'NIST Photodiode QE', 'Longpass Inserted'])
        axarr_QE[1,0].legend(full_SC_QE_plots + [PD_QE_plots[-1]], ['Inferred Solar Cell QE - ' + cell_id for cell_id in cell_ids] + ['NIST Photodiode QE (NOT our data)'])
        axarr_QE[1,1].legend(zoomed_SC_QE_plots , ['Inferred Solar Cell QE - ' + cell_id for cell_id in cell_ids])
    f_QE.suptitle('Solar Cell QE measurements - ' + extra_title_str , fontsize = title_font_size)
    f_QE.subplots_adjust(wspace=0.4)
    plt.tight_layout()
    print ('Saving QE plots...')
    f_QE.savefig(data_dir + save_file_names[-1])
    print ('Done.')
