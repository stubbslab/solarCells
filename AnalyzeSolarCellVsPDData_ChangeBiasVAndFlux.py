import numpy as np
import matplotlib.pyplot as plt
import cantrips as can
import sys

def readInSCandPDDataFile(data_file, data_dir):
    with open (data_dir + data_file, 'r') as open_file:
        data = open_file.readlines()
    data = [row.split('|') for row in data]
    zaber_positions = [float(row[0])  for row in data]
    print ('zaber_positions = ' + str(zaber_positions))
    bias_voltages = [float(row[1].split('uA/V')[0])  for row in data]
    print ('bias_voltages = ' + str(bias_voltages))
    data = [row[2:] for row in data]
    print('can.recursiveStrToListOfLists(data[0][0]) = ' + str(can.recursiveStrToListOfLists(data[0][0])))
    print ('[[ (elem[1:-1]) for elem in can.recursiveStrToListOfLists(row[0])] for row in data] = ' + str([[ (elem[1:-1]) for elem in can.recursiveStrToListOfLists(row[0])] for row in data]))
    SC_illum = [[ (float(elem[1:-1])) for elem in can.recursiveStrToListOfLists(row[0])] for row in data]
    SC_dark = [[(float(elem[1:-1])) for elem in can.recursiveStrToListOfLists(row[1])] for row in data]
    PD_illum = [[(float(elem[1:-1])) for elem in can.recursiveStrToListOfLists(row[2])] for row in data]
    PD_dark = [[(float(elem[1:-1])) for elem in can.recursiveStrToListOfLists(row[3])] for row in data]
    #SC_illum_bv = [[ (float(elem[1:-1])) for elem in can.recursiveStrToListOfLists(row[4])] for row in data]
    #SC_dark_bv = [[(float(elem[1:-1])) for elem in can.recursiveStrToListOfLists(row[5])] for row in data]
    #PD_illum_bv = [[(float(elem[1:-1])) for elem in can.recursiveStrToListOfLists(row[6])] for row in data]
    #PD_dark_bv = [[(float(elem[1:-1])) for elem in can.recursiveStrToListOfLists(row[7])] for row in data]

    #dial_positions, bias_voltages, SC_illum, SC_dark, PD_illum_bv, PD_dark_bv, SC_illum_bv, SC_dark_bv, PD_illum_bv, PD_dark_bv = can.safeSortOneListByAnother(bias_voltages, [zaber_positions, bias_voltages, SC_illum, SC_dark, PD_illum_bv, PD_dark_bv, SC_illum_bv, SC_dark_bv, PD_illum_bv, PD_dark_bv])
    #return [dial_positions, bias_voltages, SC_illum, SC_dark, PD_illum_bv, PD_dark_bv, SC_illum_bv, SC_dark_bv, PD_illum_bv, PD_dark_bv]
    dial_positions, bias_voltages, SC_illum, SC_dark, PD_illum, PD_dark = can.safeSortOneListByAnother(bias_voltages, [zaber_positions, bias_voltages, SC_illum, SC_dark, PD_illum, PD_dark])
    return [dial_positions, bias_voltages, SC_illum, SC_dark, PD_illum, PD_dark]

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
    data_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/solarCell/linearity/20210408/'
    data_files = ['SCLinearity_FluxAndBiasChanges_20210408.txt']
    bad_data_indeces = [[ ]]
    index_break = []
    file_legend_ids = ['']
    plot_colors = ['r','b']
    plot_styles = ['-','--']
    SC_plots = []
    PD_plots = []
    #data_files = ['Fiber_Config_1Cal_2Uncal_Hammamatsu2.csv' ]
    #file_legend_ids = ['CalxTip1 / UncalxTip2']
    #plot_colors = ['r']

    f, axarr = plt.subplots(1)
    for i in range(len(data_files)):
        plot_style = plot_styles[i]
        data_file = data_files[i]

        #z_positions, bias_voltages, SC_illum, SC_dark, PD_illum_bv, PD_dark_bv, SC_illum_bv, SC_dark_bv, PD_illum_bv, PD_dark_bv = readInSCandPDDataFile(data_file, data_dir)
        z_positions, bias_voltages, SC_illum, SC_dark, PD_illum, PD_dark = readInSCandPDDataFile(data_file, data_dir)
        n_bias_voltages = len([volt for volt in bias_voltages if volt == bias_voltages[0]])
        n_z_pos = len(bias_voltages) // n_bias_voltages

        print ('[SC_illum, SC_dark, PD_illum,PD_dark] = ' + str([SC_illum, SC_dark, PD_illum ,PD_dark]))
        arrs = [SC_illum, SC_dark, PD_illum, PD_dark]
        SC_illum_stats, SC_dark_stats, PD_illum_stats, PD_dark_stats = ([[[np.mean(arr[j]), np.std(arr[j])] for j in range(len(arr))] for arr in arrs])
        SC_illum_mean = [SC_illum_stats[j][0] for j in range(len(SC_illum_stats)) if not (j in bad_data_indeces[i]) ]
        SC_dark_mean = [SC_dark_stats[j][0] for j in range(len(SC_dark_stats)) if not (j in bad_data_indeces[i]) ]
        PD_illum_mean = [PD_illum_stats[j][0] for j in range(len(PD_illum_stats)) if not (j in bad_data_indeces[i]) ]
        PD_dark_mean = [PD_dark_stats[j][0] for j in range(len(PD_dark_stats)) if not (j in bad_data_indeces[i]) ]
        print ('SC_illum_stats = ' + str(SC_illum_stats))
        ratios = (np.array(SC_illum_mean) - np.array(SC_dark_mean)) / (np.array(PD_illum_mean) - np.array(PD_dark_mean))

        z_pos_mesh = np.reshape(z_positions, (n_z_pos, n_bias_voltages))
        print ('z_pos_mesh = ' + str(z_pos_mesh))
        bv_mesh = np.reshape(bias_voltages, (n_z_pos, n_bias_voltages))
        print ('bv_mesh  = ' + str(bv_mesh ))
        ratios_mesh = np.reshape(ratios, (n_z_pos, n_bias_voltages))
        print ('ratios_mesh = ' + str(ratios_mesh))

        print ('[len(SC_illum_mean), len(PD_illum_mean), len(PD_dark_mean)] = ' + str([len(SC_illum_mean), len(PD_illum_mean), len(PD_dark_mean)]))
        #low_x = (np.array(PD_illum_mean) - np.array(PD_dark_mean))[0:4].tolist() + (np.array(PD_illum_mean) - np.array(PD_dark_mean))[5:7].tolist()
        #low_y = SC_illum_mean[0:4] + SC_illum_mean[5:7]
        #high_x = [float((np.array(PD_illum_mean) - np.array(PD_dark_mean)) [4])] + (np.array(PD_illum_mean) - np.array(PD_dark_mean))[7:].tolist()
        #high_y = [SC_illum_mean[4]] + SC_illum_mean[7:]
        #print ('[len(low_x), len(low_y), len(high_x), len(high_y)] = ' + str([len(low_x), len(low_y), len(high_x), len(high_y)] ))
        #low_fit = np.polyfit(low_x, low_y, 1)
        #high_fit = np.polyfit(high_x, high_y, 1)
        SC_plots = SC_plots + [axarr.contour(z_pos_mesh, bv_mesh, ratios_mesh)]
        #axarr.scatter(high_x, high_y, linestyle = plot_style, color = 'k' )
        #x_lims = axarr.get_xlim()
        #low_fit_plot = axarr.plot(x_lims, np.poly1d(low_fit)(x_lims), c = 'r')[0]
        #high_fit_plot = axarr.plot(x_lims, np.poly1d(high_fit)(x_lims), c = 'k')[0]
        #ratios = (np.array(SC_illum_mean) - np.array(SC_dark_mean)) / (np.array(PD_illum_mean) - np.array(PD_dark_mean))
        #axarr.scatter((np.array(PD_illum_mean) - np.array(PD_dark_mean)), ratios)
        #plt.errorbar(wavelengths, SC_diffs, yerr = SC_diff_stds, fmt = 'none', c = plot_colors[0])
        #PD_plots = PD_plots + [plt.plot(wavelengths, PD_diffs, c = plot_colors[1], linestyle = plot_style)[0]]
        #plt.errorbar(wavelengths, PD_diffs, yerr = PD_diff_stds, fmt = 'none', c = plot_colors[0])

    #plt.legend(SC_plots + PD_plots, file_legend_ids)
    #plt.axvline(602, c = 'grey', alpha = 0.5)
    #plt.text(602, 0.91, '500nm longpass \n inserted', rotation = 90)
    #axarr.set_xlabel(r'PD photo current (total - dark) ($\mu A$)', fontsize = 12)
    #axarr.set_ylabel(r'SC total current ($\mu A$)', fontsize = 12)
    #axarr.legend([low_fit_plot, high_fit_plot],
    #                 [r'$R_D / (R_D + R_{ser}) \times QE_{SC}/QE_{PD} =$'+ str(can.round_to_n(low_fit[0], 3)) + '\n' + r'$(V_{ext} - V_D)/(R_{D} + R_{ser}) =$' + str(can.round_to_n(-low_fit[1], 3))  + r'$\mu A$', r'$QE_{SC}/QE_{PD}=$' + str(can.round_to_n(high_fit[0], 3)) +  '\n' + r'$V_{ext}/R_{sh}=$' + str(can.round_to_n(-high_fit[1], 3)) + r'$\mu A$' ])
    #axarr[1].set_xlabel(r'PD photo current (total - dark) ($\mu A$)', fontsize = 12)
    #axarr[1].set_ylabel(r'SC total current ($\mu A$)', fontsize = 12)
    plt.show()
