import cantrips as can
import numpy as np
import matplotlib.pyplot as plt

def plotSCDarkCurrents(dark_current_data_file, save_file_name,
                           data_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/solarCell/darkCurrent/20210430/',
                            delimiter = ',', n_header_lines = 1, n_scs_per_row = 40, xlabel = r'$\Delta t$ (s)', ylabel = r'photocurrent $(\mu A)$', title = 'Dark Current',
                            fig_height = 10, fig_width = 10, plot_color = 'k', time_sep = 0.1 + 100 / 60):
    dark_data = can.readInColumnsToList(dark_current_data_file, file_dir = data_dir, n_ignore = n_header_lines, n_ignore_end = 0, delimiter = delimiter)
    print ('dark_data = ' + str(dark_data))
    dark_currents = [10.0 ** 6.0 * float(dark_c) for dark_c in dark_data[0]]
    print ('dark_currents = ' + str(dark_currents))
    time_seps = [i * time_sep for i in range(len(dark_currents))]
    print ('np.max(np.abs(dark_currents)) = ' + str(np.max(np.abs(dark_currents)) ))
    print ('[time_seps, dark_currents] = ' + str([time_seps, dark_currents]))
    plt.scatter(time_seps, dark_currents, c = plot_color )
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.savefig(data_dir + save_file_name)

    return 1

def plotSCShuntResistances(save_file_name,
                           data_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/solarCell/MeasuredElectricalProperties/',
                           rshunt_data_file = 'ShuntResistances.csv', delimiter = ',', n_header_lines = 1, n_scs_per_row = 40, cut_off_R = 300,
                           fig_width = 15, fig_height_per_row = 3, batch_colors = {1:'green', 2:'blue', 3:'magenta'}, labelsize = 14,
                           yticks = [0, 1, 2, 3]):
    rshunt_data = can.readInColumnsToList(rshunt_data_file, file_dir = data_dir, n_ignore = n_header_lines, n_ignore_end = 0, delimiter = delimiter)
    #print ('[[elem[0], len(elem)] for elem in rshunt_data] = ' + str([[elem[0], len(elem)] for elem in rshunt_data]))
    cell_ids = rshunt_data[0]
    rshunts = [float(elem) for elem in rshunt_data[1]]
    max_shunt = max(rshunts)
    cshunts = [float(elem) for elem in rshunt_data[2]]
    good_data = [int(elem) for elem in rshunt_data[3]]
    batch_ids = [int(elem) for elem in rshunt_data[4]]
    unique_batches = np.unique(batch_ids)
    n_batches = len(unique_batches)
    n_cells = len(cell_ids)
    n_rows = int(np.ceil(n_cells / n_scs_per_row))
    f, axarr = plt.subplots(n_rows, 1, figsize = (fig_width, (fig_height_per_row * n_rows)) )
    scat_id_dict = {}
    for i in range(n_rows):
        to_plot_indeces = [ i * n_scs_per_row, min((i + 1) * n_scs_per_row, len(cell_ids)) ]
        ids_to_plot = cell_ids[to_plot_indeces[0]:to_plot_indeces[1]] + [''  for i in range(n_scs_per_row - (to_plot_indeces[1] - to_plot_indeces[0])) if i > 0]
        rshunts_to_plot = rshunts[to_plot_indeces[0]:to_plot_indeces[1]] +  [-100 for i in range(n_scs_per_row - (to_plot_indeces[1] - to_plot_indeces[0])) if i > 0]
        good_data_to_plot = good_data[to_plot_indeces[0]:to_plot_indeces[1]] +  [0 for i in range(n_scs_per_row - (to_plot_indeces[1] - to_plot_indeces[0])) if i > 0]
        batch_ids_to_plot = batch_ids[to_plot_indeces[0]:to_plot_indeces[1]] +  [0 for i in range(n_scs_per_row - (to_plot_indeces[1] - to_plot_indeces[0])) if i > 0]
        for batch in unique_batches:
            scat_id_dict[(batch, 0)] = axarr[i].scatter([j for j in range(len(rshunts_to_plot)) if good_data_to_plot[j] == 0 and batch_ids_to_plot[j] == batch], [rshunts_to_plot[j] / 1000 for j in range(len(rshunts_to_plot)) if good_data_to_plot[j] == 0 and batch_ids_to_plot[j] == batch], marker = 'o', c = batch_colors[batch])
            scat_id_dict[(batch, 1)] = axarr[i].scatter([j for j in range(len(rshunts_to_plot)) if good_data_to_plot[j] == 1 and batch_ids_to_plot[j] == batch], [rshunts_to_plot[j] / 1000 for j in range(len(rshunts_to_plot)) if good_data_to_plot[j] == 1 and batch_ids_to_plot[j] == batch], marker = 'x', c = batch_colors[batch])
            [axarr[i].axvline(j , linestyle = '--', c = 'gray',alpha = 0.25) for j in range(len(rshunts_to_plot)) ]
            axarr[i].axhline(cut_off_R / 1000, linestyle = ':', c = 'red',alpha = 0.5)
        axarr[i].set_xticks(range(len(rshunts_to_plot)))
        axarr[i].set_xticklabels(ids_to_plot, fontsize = labelsize)
        axarr[i].set_yticks(yticks)
        axarr[i].set_yticklabels(yticks, fontsize = labelsize)
        axarr[i].set_ylabel(r'SC Shunt Resistance (k$\Omega$)', fontsize = labelsize)
        #axarr[i].set_yticklabels(ids_to_plot + ['-' + str(i) for i in range(n_scs_per_col - len(ids_to_plot)) if i > 0] )
        axarr[i].set_ylim(0.0, max_shunt * 1.1 / 1000.0 )

    legend_types = list(scat_id_dict.keys())
    plt.legend([scat_id_dict[key] for key in legend_types ], ['Batch ' + str(key[0]) + ' - ' + ('BAD' if key[1] else 'GOOD') for key in legend_types ], fontsize = labelsize)
    plt.tight_layout()

    plt.savefig(data_dir + save_file_name)

    return 1


if __name__=="__main__":
    save_file_name = 'C60ShuntResistances.pdf'
    plotSCShuntResistances(save_file_name)
