import numpy as np
import matplotlib.pyplot as plt
import cantrips as can
import sys
import AstronomicalParameterArchive as apa
import scipy

def readInSCandPDDataFile(data_file, data_dir, SC_neg_current, PD_neg_current, ignore_str = 'DELETED', indeces_to_include = 'all'):
    with open (data_dir + data_file, 'r') as open_file:
        data = open_file.readlines()
    data = [row.split('|') for row in data]
    print ('data[0] = ' + str(data[0]))
    print ('len(data) = ' + str(len(data)))
    print ('len(data[0]) = ' + str(len(data[0])))
    wavelengths = [(row[0])  for row in data]
    for i in range(len(data)):
        row = data[i]
        print ('i = ' + str(i))
        print ('row[12] = ' + str(row[12]))
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
    #SC_temps = [float(row[13]) for row in data]
    SC_temps = [-1 for row in data]

    arrs = [wavelengths, SC_ill_times, SC_illum, SC_illum_bv, SC_dk_times, SC_dark, SC_dark_bv, PD_ill_times, PD_illum, PD_illum_bv, PD_dk_times, PD_dark, PD_dark_bv, SC_temps]
    if indeces_to_include != 'all':
        arrs = [arr[indeces_to_include[0]:indeces_to_include[1]] for arr in arrs]
    #SC_illum_bv = [[ (float(elem[1:-1])) for elem in can.recursiveStrToListOfLists(row[4])] for row in data]
    #SC_dark_bv = [[(float(elem[1:-1])) for elem in can.recursiveStrToListOfLists(row[5])] for row in data]
    #PD_illum_bv = [[(float(elem[1:-1])) for elem in can.recursiveStrToListOfLists(row[6])] for row in data]
    #PD_dark_bv = [[(float(elem[1:-1])) for elem in can.recursiveStrToListOfLists(row[7])] for row in data]

    #dial_positions, bias_voltages, SC_illum, SC_dark, PD_illum_bv, PD_dark_bv, SC_illum_bv, SC_dark_bv, PD_illum_bv, PD_dark_bv = can.safeSortOneListByAnother(bias_voltages, [zaber_positions, bias_voltages, SC_illum, SC_dark, PD_illum_bv, PD_dark_bv, SC_illum_bv, SC_dark_bv, PD_illum_bv, PD_dark_bv])
    #return [dial_positions, bias_voltages, SC_illum, SC_dark, PD_illum_bv, PD_dark_bv, SC_illum_bv, SC_dark_bv, PD_illum_bv, PD_dark_bv]
    #wavelengths, SC_ill_times, SC_illum, SC_illum_bv, SC_dk_times, SC_dark, SC_dark_bv, PD_ill_times, PD_illum, PD_illum_bv, PD_dk_times, PD_dark, PD_dark_bv, SC_temps = can.safeSortOneListByAnother(wavelengths, [wavelengths, SC_ill_times, SC_illum, SC_illum_bv, SC_dk_times, SC_dark, SC_dark_bv, PD_ill_times, PD_illum, PD_illum_bv, PD_dk_times, PD_dark, PD_dark_bv, SC_temps ])
    return arrs


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

    stems = ['SC_QE_NewCell_']
    #suffixes = ['_1', '_1', '_MultiDay']
    suffixes = ['']
    #extra_save_suffixes = ['', '', '']
    extra_save_suffixes = ['_LongQECurve']
    date_strs = ['20220401']
    extra_title_str = 'Cell ID ??'
    #cell_ids = [ 'Cell ED - ' + date_strs[0],  'Cell ED - ' + date_strs[1], 'Cell ED - ' + date_strs[2]]
    #cell_ids = ['Cell ED - ' + date_strs[0]]
    cell_ids = ['Cell ??' ]
    shunt_resistances = [340, 340, 340]
    shunt_resistances = [2000]
    plot_wavelengths = 1
    #if we don't plot wavelengths, we plot in times
    data_file_indeces_to_include = 'all'
    adj_vs_overall_std_ratio_for_inclusion = 0.0 #0.8

    zoomed_in_QE_lims = [0.89, 1.01]
    lp_filter_wavelength = 550
    wavelength_correction = -5
    n_qe_spline_points = 21
    insert_longpass = 0
    SC_neg_current = 0
    PD_neg_current = 0
    deleted_element_id = 'DELETED'
    data_files = [ stems[i] + date_strs[i] + suffixes[i] + '.txt' for i in range(len(stems)) ]
    file_roots = [f[0:-len('.txt')] for f in data_files]
    save_file_names = [file_roots[i] + extra_save_suffixes[i] + '.pdf' for i in range(len(file_roots))]
    save_data_file_names = [file_roots[i] + extra_save_suffixes[i] + '.txt' for i in range(len(file_roots))]
    xlims = [340, 1060]
    title_font_size = 24
    labelsize = 17
    ticksize = 15
    plt.rc('xtick',labelsize=ticksize)
    plt.rc('ytick',labelsize=ticksize)

    data_dirs = [dir_root + date_str + '/' for date_str in date_strs]
    save_plot_files = ['SC_QE_measurement_' + date_strs[i] + extra_save_suffixes[i] + '.pdf' for i in range(len(date_strs))]
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
    PD_colors, SC_colors, combined_colors, indicator_color = [['b', 'yellow', 'orange', 'pink', 'purple'], ['r', 'cyan', 'brown', 'magenta', 'green'], ['purple', 'green', 'grey', 'gold', 'r'], 'grey']
    plot_styles = ['-','--']
    SC_plots = []
    PD_plots = []
    PD_QE_plots = []
    full_SC_QE_plots = []
    zoomed_SC_QE_plots = []
    ratio_plots = []
    burn_in_PD = 5
    burn_in_SC = 5
    x_lims = [0.0, 0.0] #Initialized here - this will be filled in by the code.
    #An indicator of a bad observation can be a non-flat photocurrent.
    #  We can, roughly, identify this by looking at the relative std of
    #   all data points in photocurrent (PD and SC) data vector vs the
    #   std of adjacent data points.
    #  A very flat data vector should have this ratio be ~1.  A time
    #   varying data vector has this ratio < 1.  This threshold
    #   determines how willing we are to accept potentially bad data
    #   vs excluding possibly good data.

    stacked_wavelengths = []
    stacked_unique_SC_x_vals = []
    stacked_ratios = []

    f_QE, axarr_QE = plt.subplots(2,2, figsize = [15, 15],  )
    for i in range(len(data_files)):
        print ('Working on ' + str(i+1) + 'th data file of ' + str(len(data_files)))
        Rsh = shunt_resistances[i]

        data_file = data_files[i]
        file_root = file_roots[i]
        extra_save_suffix = extra_save_suffixes[i]
        PD_color = PD_colors[i]
        SC_color = SC_colors[i]
        combined_color = combined_colors[i]
        data_dir = data_dirs[i]
        date_str = date_strs[i]

        #z_positions, bias_voltages, SC_illum, SC_dark, PD_illum_bv, PD_dark_bv, SC_illum_bv, SC_dark_bv, PD_illum_bv, PD_dark_bv = readInSCandPDDataFile(data_file, data_dir)

        wavelengths, SC_illum_times, SC_illum, SC_illum_bv, SC_dark_times, SC_dark, SC_dark_bv, PD_illum_times, PD_illum, PD_illum_bv, PD_dark_times, PD_dark, PD_dark_bv, SC_temps = readInSCandPDDataFile(data_file, data_dir, SC_neg_current, PD_neg_current, ignore_str = deleted_element_id, indeces_to_include = data_file_indeces_to_include)
        wavelengths = [float(wave) + wavelength_correction for wave in wavelengths]
        print ('SC_temps = ' + str(SC_temps))
        print ('wavelengths = ' + str(wavelengths))
        bad_measurements = [0 for wave in wavelengths]

        full_SC_illum_times, full_PD_illum_times, full_SC_dark_times, full_PD_dark_times = [[np.array(wave_set) for wave_set in arr] for arr in [SC_illum_times, PD_illum_times, SC_dark_times, PD_dark_times ]]
        print ('Here 1')
        full_SC_illum, full_SC_dark =  [[np.array(wave_set) * (-1 if SC_neg_current else 1) * 10.0 ** 6.0 for wave_set in arr] for arr in [SC_illum, SC_dark]]
        #full_PD_illum, full_PD_dark = [np.array([wave_set for wave_set in PD_illum]) * (-1 if PD_neg_current else 1) * 10.0 ** 6.0, np.array([wave_set for wave_set in PD_dark]) * (-1 if PD_neg_current else 1) * 10.0 ** 6.0]
        full_PD_illum, full_PD_dark = [[np.array(wave_set) * (-1 if PD_neg_current else 1) * 10.0 ** 6.0 for wave_set in arr] for arr in [PD_illum, PD_dark]]
        print ('Here 2')
        full_SC_illum_bv, full_SC_dark_bv = [[np.array(wave_set) * (-1 if SC_neg_current else 1) * 10.0 ** 6.0 for wave_set in arr] for arr in [SC_illum_bv, SC_dark_bv]]
        print ('Here 3')
        full_PD_illum_bv, full_PD_dark_bv = [[np.array(wave_set) * (-1 if PD_neg_current else 1) * 10.0 ** 6.0 for wave_set in arr] for arr in [PD_illum_bv, PD_dark_bv]]
        print ('Here 4')
        arrs_to_display = [full_SC_illum, full_PD_illum, full_SC_dark, full_PD_dark, full_SC_illum_bv, full_PD_illum_bv, full_SC_dark_bv, full_PD_dark_bv ]
        print ('Here 5')
        #clipped_sigs_of_data_type = [can.sigClipStd(can.flattenListOfLists(arr), sig_clip = 3.0) for arr in arrs_to_display ]
        clipped_sigs_of_data_type = [np.std(arr) for arr in arrs_to_display ]
        print ('clipped_sigs_of_data_type = ' + str(clipped_sigs_of_data_type))
        SC_ill_ylims, PD_ill_ylims, SC_dk_ylims, PD_dk_ylims, SC_ill_bv_ylims, PD_ill_bv_ylims, SC_dk_bv_ylims, PD_dk_bv_ylims = [[np.median(arrs_to_display[i]) - clipped_sigs_of_data_type[i] * 3.0, np.median(arrs_to_display[i]) + clipped_sigs_of_data_type[i] * 3.0] for i in range(len(arrs_to_display)) ]
        print ('[SC_ill_ylims, PD_ill_ylims, SC_dk_ylims, PD_dk_ylims, SC_ill_bv_ylims, PD_ill_bv_ylims, SC_dk_bv_ylims, PD_dk_bv_ylims] = ' + str([SC_ill_ylims, PD_ill_ylims, SC_dk_ylims, PD_dk_ylims, SC_ill_bv_ylims, PD_ill_bv_ylims, SC_dk_bv_ylims, PD_dk_bv_ylims] ))

        #We should correct each data point in the measured currents with the measured voltages, since the measured voltages drift a little
        #full_SC_illum = np.array(full_SC_illum_bv) / Rsh + np.array(full_SC_illum)
        #full_SC_dark = np.array(full_SC_dark_bv) / Rsh + np.array(full_SC_dark)

        SC_illum_times, SC_illum, SC_illum_bv, SC_dark_times, SC_dark, SC_dark_bv = [[wave_set[burn_in_SC:] for wave_set in arr] for arr in [full_SC_illum_times, full_SC_illum, full_SC_illum_bv, full_SC_dark_times, full_SC_dark, full_SC_dark_bv]]
        PD_illum_times, PD_illum, PD_illum_bv, PD_dark_times, PD_dark, PD_dark_bv = [[wave_set[burn_in_PD:] for wave_set in arr] for arr in [full_PD_illum_times, full_PD_illum, full_PD_illum_bv, full_PD_dark_times, full_PD_dark, full_PD_dark_bv]]

        print ('Here 6')
        n_rows = min(8, len(wavelengths))
        n_single_shot_cols = int(np.ceil(len(wavelengths) / n_rows))
        """
        f_PD_ill, axarr_PD_ill = plt.subplots(n_rows, n_single_shot_cols, figsize = (n_single_shot_cols * 3, n_rows * 3), sharex = True, squeeze = False)
        f_PD_dk, axarr_PD_dk = plt.subplots(n_rows, n_single_shot_cols, figsize = (n_single_shot_cols * 3, n_rows * 3), sharex = True, squeeze = False)
        f_SC_ill, axarr_SC_ill = plt.subplots(n_rows, n_single_shot_cols, figsize = (n_single_shot_cols * 3, n_rows * 3), sharex = True, squeeze = False)
        f_SC_dk, axarr_SC_dk = plt.subplots(n_rows, n_single_shot_cols, figsize = (n_single_shot_cols * 3, n_rows * 3), sharex = True, squeeze = False)
        f_PD_bv_ill, axarr_PD_bv_ill = plt.subplots(n_rows, n_single_shot_cols, figsize = (n_single_shot_cols * 3, n_rows * 3), sharex = True, squeeze = False)
        f_PD_bv_dk, axarr_PD_bv_dk = plt.subplots(n_rows, n_single_shot_cols, figsize = (n_single_shot_cols * 3, n_rows * 3), sharex = True, squeeze = False)
        f_SC_bv_ill, axarr_SC_bv_ill = plt.subplots(n_rows, n_single_shot_cols, figsize = (n_single_shot_cols * 3, n_rows * 3), sharex = True, squeeze = False)
        f_SC_bv_dk, axarr_SC_bv_dk = plt.subplots(n_rows, n_single_shot_cols, figsize = (n_single_shot_cols * 3, n_rows * 3), sharex = True, squeeze = False)
        """
        if plot_wavelengths:
            data_set_strs = wavelengths[:]
            true_wavelengths = [float(wave) for wave in wavelengths]
            wavelengths = [float(wave) for wave in wavelengths]
            wavelength_strs = [str(wave) + 'nm' for wave in wavelengths]
            SC_x_vals = PD_x_vals = wavelengths[:]
            plot_strs = [r'$\lambda=$' + str(wave) + 'nm' for wave in wavelengths]
            xlabel = r'Monochromator Wavelength (nm)'
            min_wavelength, max_wavelength = [min(wavelengths), max(wavelengths)]
            delta_wavelength = max_wavelength - min_wavelength
            x_lims = [min_wavelength - delta_wavelength * 0.01, max_wavelength + delta_wavelength * 0.01]
        else:
            wavelength_strs = wavelengths[:]
            true_wavelengths = [650 for i in range(len(wavelengths))]
            mean_SC_exp_times = np.mean(SC_illum_times, axis = 1)
            mean_PD_exp_times = np.mean(PD_illum_times, axis = 1)
            earliest = np.min(mean_SC_exp_times.tolist() + mean_PD_exp_times.tolist())
            SC_x_vals = (mean_SC_exp_times- earliest).tolist()
            PD_x_vals = (mean_PD_exp_times - earliest).tolist()
            plot_strs = ['Obs Number ' + str(i) for i in range(len(wavelengths))]
            xlabel = r'$\Delta t$ (s)'
            x_lims = [0.0, max(SC_x_vals + PD_x_vals + [x_lims[1]])]
        unique_SC_x_vals = np.unique(SC_x_vals)


        #all_data_axes = [axarr_PD_ill, axarr_SC_ill, axarr_PD_dk, axarr_SC_dk, axarr_PD_bv_ill, axarr_SC_bv_ill, axarr_PD_bv_dk, axarr_SC_bv_dk ]
        #all_data_figures = [f_PD_ill, f_SC_ill, f_PD_dk, f_SC_dk, f_PD_bv_ill, f_SC_bv_ill, f_PD_bv_dk, f_SC_bv_dk]
        all_data_x_sets = [[np.array(time_set) - time_set[0] for time_set in full_PD_illum_times], [np.array(time_set) - time_set[0] for time_set in full_SC_illum_times],
                           [np.array(time_set) - time_set[0] for time_set in full_PD_dark_times], [np.array(time_set) - time_set[0] for time_set in full_SC_dark_times],
                            [np.array(time_set) - time_set[0] for time_set in full_PD_illum_times], [np.array(time_set) - time_set[0] for time_set in full_SC_illum_times ],
                             [np.array(time_set) - time_set[0] for time_set in full_PD_dark_times], [np.array(time_set) - time_set[0] for time_set in full_SC_dark_times ] ]
        print ('full_PD_illum_times[0] = ' + str(full_PD_illum_times[0] ))
        all_data_y_sets = [full_PD_illum, full_SC_illum, full_PD_dark, full_SC_dark, full_PD_illum_bv, full_SC_illum_bv, full_PD_dark_bv, full_SC_dark_bv ]
        ylims_set = [PD_ill_ylims, SC_ill_ylims,  PD_dk_ylims, SC_dk_ylims,  PD_ill_bv_ylims,  SC_ill_bv_ylims,  PD_dk_bv_ylims,  SC_dk_bv_ylims ]
        full_data_ylabels = [r'Photocurrent ($\mu A$)' for i in range(4)] + [r'Output Voltage ($\mu V$)' for i in range(4)]
        full_data_xlabels = [r'$\Delta t$ (s)' for i in range(8)]
        suptitles = ['PD illuminated photocurrent - ' + date_str, 'SC illuminated photocurrent - ' + date_str, 'PD dark photocurrent - ' + date_str, 'SC dark photocurrent - ' + date_str,
                     'PD illuminated output voltage - ' + date_str, 'SC illuminated output voltage - ' + date_str, 'PD dark output voltage - ' + date_str, 'SC dark output voltage - ' + date_str]
        #suptitles = ['PD illuminated photocurrent - ' + date_str, 'SC illuminated photocurrent - ' + date_str,
        #             'PD illuminated output voltage - ' + date_str, 'SC illuminated output voltage - ' + date_str ]
        plot_suffixes = [ PD_ill_plot_suffix, SC_ill_plot_suffix, PD_dark_plot_suffix, SC_dark_plot_suffix, PD_ill_bv_plot_suffix, SC_ill_bv_plot_suffix, PD_dark_bv_plot_suffix, SC_dark_bv_plot_suffix]
        v_line_indeces = [burn_in_PD, burn_in_SC, burn_in_PD, burn_in_SC, burn_in_PD, burn_in_SC, burn_in_PD, burn_in_SC]

        print ('Making plots for SC and PD traces...')
        for k in range(len([full_PD_illum, full_SC_illum])):
            #Cut out measurements with large drift term
            print ('bad_measurements = ' + str(bad_measurements))
            all_data_y_set = [full_PD_illum, full_SC_illum][k]
            v_line_index = [burn_in_PD, burn_in_SC][k]
            for j in range(len(wavelengths)):
                y_set = all_data_y_set[j]
                full_y_std = np.std(y_set[v_line_index:])
                adj_y_std = np.std([y_set[i] - y_set[i-1] for i in range(v_line_index, len(y_set))])
                y_adj_to_full_std_ratio = adj_y_std / full_y_std
                #print ('wavelengths: ' + str(wavelengths[j]) + '-> y_adj_to_full_std_ratio: ' + str(y_adj_to_full_std_ratio))
                bad_measurements[j] = (bad_measurements[j] or (y_adj_to_full_std_ratio < adj_vs_overall_std_ratio_for_inclusion))

        bad_data_indeces[i] = bad_data_indeces[i] + [index for index in range(len(bad_measurements)) if bad_measurements[index]]
        print ('bad_data_indeces[i] = ' + str(bad_data_indeces[i]))
        for k in range(len(all_data_x_sets)):
            print ('k = ' + str(k))
            single_data_figure, data_axis_arr = plt.subplots(n_rows, n_single_shot_cols, figsize = (n_single_shot_cols * 1, n_rows * 1), sharex = True, squeeze = False)
            #data_axis_arr = all_data_axes [k]
            all_data_x_set = all_data_x_sets[k]
            all_data_y_set = all_data_y_sets[k]
            full_data_ylabel = full_data_ylabels[k]
            full_data_xlabel = full_data_xlabels[k]
            for j in range(len(wavelengths)):
                print ('Working on j = ' + str(j) + ' of ' + str(len(wavelengths)) )
                bad_measurement = bad_measurements[j]
                wave = wavelengths[j]
                x_set = all_data_x_set[j]
                y_set = all_data_y_set[j]
                full_y_std = np.std(y_set[v_line_indeces[k]:])
                adj_y_std = np.std([y_set[i] - y_set[i-1] for i in range(v_line_indeces[k], len(y_set))])
                y_adj_to_full_std_ratio = adj_y_std / full_y_std
                wavelength_str = wavelength_strs[j] + r' - $\sigma_A/\sigma_F$=' + str(can.round_to_n(y_adj_to_full_std_ratio, 3))
                row_num, col_num = [int(j // n_single_shot_cols), int(j % n_single_shot_cols)]
                data_axis_arr [row_num, col_num].plot(x_set, y_set, color = 'k')
                data_axis_arr[row_num, col_num].axvline(x_set[v_line_indeces[k]], linestyle = '--', c = 'grey', alpha = 0.5)
                #data_axis_arr[row_num, col_num].set_ylim(ylims)
                data_axis_arr [row_num, col_num].text(0.1, 0.9, wavelength_str, transform=data_axis_arr [row_num, col_num].transAxes, color = ('r' if bad_measurement else 'g') )
                if col_num == 0:
                    data_axis_arr [row_num, col_num].set_ylabel(full_data_ylabel, fontsize = labelsize)
                if row_num == n_rows-1:
                    data_axis_arr [row_num, col_num].set_xlabel(full_data_xlabel, fontsize = labelsize)

            #single_data_type_axes = all_data_axes[k]
            #single_data_figure = all_data_figures [k]
            single_data_figure.suptitle(suptitles[k])
            single_data_figure.tight_layout()
            single_data_figure.savefig(data_dir + file_root + extra_save_suffix + plot_suffixes[k])
            print ('Done making ' + str(k) + 'th plot for SC and PD traces of ' + str(len(all_data_x_sets) ))

        arrs = [SC_illum, SC_dark, PD_illum, PD_dark]
        SC_illum_stats, SC_dark_stats, PD_illum_stats, PD_dark_stats = ([[[np.mean(arr[j]) * current_scaling_v_to_uA, np.std(arr[j]) / np.sqrt(len(arr[j])) * current_scaling_v_to_uA] for j in range(len(arr))] for arr in arrs])
        #SC_illum_stats, PD_illum_stats,  = ([[[np.mean(arr[j]) * current_scaling_v_to_uA, np.std(arr[j]) / np.sqrt(len(arr[j])) * current_scaling_v_to_uA] for j in range(len(arr))] for arr in arrs])
        PD_x_vals = [PD_x_vals[j] for j in range(len(PD_x_vals)) if not (j in bad_data_indeces[i]) ]
        SC_x_vals = [SC_x_vals[j] for j in range(len(SC_x_vals)) if not (j in bad_data_indeces[i]) ]
        true_wavelengths = [true_wavelengths[j] for j in range(len(true_wavelengths)) if not (j in bad_data_indeces[i]) ]
        SC_illum_mean = [SC_illum_stats[j][0] for j in range(len(SC_illum_stats)) if not (j in bad_data_indeces[i]) ]
        SC_illum_std = [SC_illum_stats[j][1] for j in range(len(SC_illum_stats)) if not (j in bad_data_indeces[i]) ]
        SC_dark_mean = [SC_dark_stats[j][0] for j in range(len(SC_dark_stats)) if not (j in bad_data_indeces[i]) ]
        SC_dark_std = [SC_dark_stats[j][1] for j in range(len(SC_dark_stats)) if not (j in bad_data_indeces[i]) ]
        PD_illum_mean = [PD_illum_stats[j][0] for j in range(len(PD_illum_stats)) if not (j in bad_data_indeces[i]) ]
        PD_illum_std = [PD_illum_stats[j][1] for j in range(len(PD_illum_stats)) if not (j in bad_data_indeces[i]) ]
        PD_dark_mean = [PD_dark_stats[j][0] for j in range(len(PD_dark_stats)) if not (j in bad_data_indeces[i]) ]
        PD_dark_std = [PD_dark_stats[j][1] for j in range(len(PD_dark_stats)) if not (j in bad_data_indeces[i]) ]
        #SC_diff = (np.array(SC_illum_mean) - 0.0)
        SC_diff = (np.array(SC_illum_mean) - np.array(SC_dark_mean))
        SC_diff_std = np.sqrt(np.array(SC_illum_std) ** 2.0 + np.array(SC_dark_std) ** 2.0)
        #SC_diff_std = np.array(SC_illum_std)
        PD_diff = (np.array(PD_illum_mean) - np.array(PD_dark_mean)) 
        PD_diff_std = np.sqrt(np.array(PD_illum_std) ** 2.0 + np.array(PD_dark_std) ** 2.0)
        #PD_diff_std = np.array(PD_illum_std)
        ratios = SC_diff / PD_diff
        ratio_errs = np.sqrt((SC_diff_std ** 2.0 * PD_diff ** -2.0) + (SC_diff ** 2.0 * PD_diff_std ** 2.0 ) * PD_diff ** -4.0)


        #print ('[qe_spline_wavelengths, qe_spline_ratios, PD_QE_interp(qe_spline_wavelengths) * qe_spline_ratios] = ' + str([qe_spline_wavelengths, qe_spline_ratios, PD_QE_interp(qe_spline_wavelengths) * qe_spline_ratios] ))

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
        lp_in_first = axarr_QE[0,0].axvline(lp_filter_wavelength, color = indicator_color, alpha = 0.5, linestyle = 'dashed')
        ratio_plots = ratio_plots + [axarr_QE[0,1].scatter(SC_x_vals, ratios, c = combined_color, marker = '.')]
        axarr_QE[0,1].errorbar(SC_x_vals, ratios, yerr = ratio_errs, fmt = 'none', color = combined_color)
        lp_in_second = axarr_QE[0,1].axvline(lp_filter_wavelength, color = indicator_color, alpha = 0.5, linestyle = 'dashed')
        #PD_QE_plots = PD_QE_plots + [axarr_QE[1,0].plot(SC_x_vals, PD_QE_interp(true_wavelengths), c = PD_color, marker = '.')[0]]
        PD_QE_plots = PD_QE_plots + [axarr_QE[1,0].scatter(SC_x_vals, PD_QE_interp(true_wavelengths), color = PD_color, marker = '.')]
        full_SC_QE_plots = full_SC_QE_plots + [axarr_QE[1,0].scatter(SC_x_vals, PD_QE_interp(true_wavelengths) * np.array(ratios), c = SC_color, marker = '.')]
        #axarr_QE[1,0].scatter(qe_spline_wavelengths, PD_QE_interp(qe_spline_wavelengths) * qe_spline_ratios, marker = 'x', c ='k')
        axarr_QE[1,1].axhline(1.0, 0, 1, color = 'k', linestyle = 'dashed')
        zoomed_SC_QE_plots = zoomed_SC_QE_plots + [axarr_QE[1,1].scatter(SC_x_vals, PD_QE_interp(true_wavelengths) * np.array(ratios), c = SC_color, marker = '.')]
        axarr_QE[1,1].errorbar(SC_x_vals, PD_QE_interp(true_wavelengths) * np.array(ratios), yerr = PD_QE_interp(true_wavelengths) * np.array(ratio_errs), fmt = 'none', color = SC_color)
        #axarr_QE[1,1].scatter(qe_spline_wavelengths, PD_QE_interp(qe_spline_wavelengths) * qe_spline_ratios, marker = 'x', c ='k')

        axarr_QE[0,0].set_xlim(x_lims)
        axarr_QE[0,1].set_xlim(x_lims)
        axarr_QE[1,0].set_xlim(x_lims)
        axarr_QE[1,1].set_xlim(x_lims)

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
        stacked_wavelengths = stacked_wavelengths + wavelengths
        stacked_unique_SC_x_vals = np.unique(stacked_unique_SC_x_vals + unique_SC_x_vals.tolist() ).tolist()
        stacked_ratios = stacked_ratios + ratios.tolist()

    """
    qe_spline_indeces = [np.argmin(abs(np.array(stacked_wavelengths) - spline_wave)) for spline_wave in np.linspace(min(stacked_wavelengths), max(stacked_wavelengths), n_qe_spline_points) ]
    qe_spline_wavelengths = [stacked_wavelengths[index] for index in qe_spline_indeces]
    qe_spline_ratios = [np.mean([stacked_ratios[i] for i in range(len(stacked_ratios)) if stacked_wavelengths[i] == wave]) for wave in qe_spline_wavelengths]
    print ('qe_spline_wavelengths = ' + str(qe_spline_wavelengths))
    print ('qe_spline_ratios = ' + str(qe_spline_ratios))
    ratio_cubic_interp = scipy.interpolate.interp1d(qe_spline_wavelengths, qe_spline_ratios, kind = 'cubic')
    """

    #axarr_QE[1,0].fill_between(stacked_unique_SC_x_vals, PD_QE_interp(stacked_unique_SC_x_vals) * ratio_cubic_interp(stacked_unique_SC_x_vals) * 0.995, PD_QE_interp(stacked_unique_SC_x_vals) * ratio_cubic_interp(stacked_unique_SC_x_vals) * 1.005, color ='k', alpha = 0.2)
    #axarr_QE[1,1].fill_between(stacked_unique_SC_x_vals, PD_QE_interp(stacked_unique_SC_x_vals) * ratio_cubic_interp(stacked_unique_SC_x_vals) * 0.995, PD_QE_interp(stacked_unique_SC_x_vals) * ratio_cubic_interp(stacked_unique_SC_x_vals) * 1.005, color ='k', alpha = 0.2)

    if insert_longpass:
        axarr_QE[0,0].legend(PD_plots + SC_plots + [lp_in_first], ['PD photocurrent: ' + cell_id for cell_id in cell_ids] + ['SC photocurrent: ' + cell_id for cell_id in cell_ids] + ['>600nm Longpass Inserted'] , fontsize = labelsize)
        axarr_QE[0,1].legend(ratio_plots + [lp_in_second], ['SC/PD photocurrent ratio: ' + cell_id for cell_id in cell_ids] + ['>600nm Longpass Inserted'], fontsize = labelsize)
#axarr_QE[2].legend([SC_QE_plots[0], PD_QE_plots[0], lp_in_third], ['Inferred Solar Cell QE', 'NIST Photodiode QE', 'Longpass Inserted'])
        axarr_QE[1,0].legend([PD_QE_plots[-1]] + full_SC_QE_plots, ['NIST Photodiode QE (NOT our data)'] + ['Inferred SC QE: ' + cell_id for cell_id in cell_ids], fontsize = labelsize)
        axarr_QE[1,1].legend(zoomed_SC_QE_plots + [lp_in_forth], ['Inferred SC QE: ' + cell_id for cell_id in cell_ids] + ['>600nm Longpass Inserted'], fontsize = labelsize)
    else:
        axarr_QE[0,0].legend(PD_plots + SC_plots,  ['PD photocurrent: ' + cell_id for cell_id in cell_ids] + ['SC photocurrent: ' + cell_id for cell_id in cell_ids]  , fontsize = labelsize)
        axarr_QE[0,1].legend(ratio_plots  , ['SC/PD photocurrent ratio: ' + cell_id for cell_id in cell_ids], fontsize = labelsize )
#axarr_QE[2].legend([SC_QE_plots[0], PD_QE_plots[0], lp_in_third], ['Inferred Solar Cell QE', 'NIST Photodiode QE', 'Longpass Inserted'])
        axarr_QE[1,0].legend([PD_QE_plots[-1]] + full_SC_QE_plots, ['NIST PD QE (NOT our data)'] + ['Inferred SC QE: ' + cell_id for cell_id in cell_ids], fontsize = labelsize)
        axarr_QE[1,1].legend(zoomed_SC_QE_plots , ['Inferred SC QE: ' + cell_id for cell_id in cell_ids], fontsize = labelsize)
    f_QE.suptitle('SC QE measurements: ' + extra_title_str , fontsize = title_font_size)
    f_QE.subplots_adjust(wspace=0.4)
    plt.tight_layout()
    print ('Saving QE plots...')
    f_QE.savefig(data_dir + save_file_names[-1])
    can.saveListsToColumns([SC_x_vals, PD_QE_interp(true_wavelengths) * np.array(ratios), PD_QE_interp(true_wavelengths) * np.array(ratio_errs)],  save_data_file_names[-1], data_dir, header = ['Wavelength (nm), QE (unitless), QE Err (unitless)'], sep = ', ')
    print ('Done.')
