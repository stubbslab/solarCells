import numpy as np 
import cantrips as c 
import matplotlib.pyplot as plt
import matplotlib  
import scipy.interpolate as interpolate 
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def plotFullLinearitySet(LEDIDs, target_dirs, LI_sens_strs_set, ammeter_A_to_Vs_set, fiber_strs_set, 
                         save_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/solarCell/linearity/',
                         wave_info_str = '_13HzVpp1200mVOffset600mVSqr', filt_info_str = 'NDFiltnone', 
                         file_funct = lambda LED_name, fiber_str, LI_sens_str, ammeter_str, wave_stuff_str, filt_stuff_str: 'linearity_SCandPD_byHand_LEDid' + LED_name + wave_stuff_str + '_fiber' + fiber_str + '_LISens' + LI_sens_str + 'mV_ammeterGain' + ammeter_str + '_' + filt_stuff_str + '.txt',
                         fig_size_units = [4.0, 2.0] , LED_point_colors = ['b', 'r','g','orange','m'], LED_err_colors =  ['b', 'r','g','orange','m'], 
                         markers = ['v', '^'], stacked = 0,
                         pad = 5.0, w_pad = 5.0, h_pad = 5.0, ignore_line_char = None, 
                         save_plot_name = 'solarCellLinearities.png', show = 0, save = 1, fiber_save_strs_set = 'default' ):

    max_nFibers = max([len(fiber_strs) for fiber_strs in fiber_strs_set]) 
    n_LEDs = len(LEDIDs) 
    f, axarr = plt.subplots(n_LEDs, 1, figsize = [max_nFibers * fig_size_units[0], n_LEDs * fig_size_units[1]], squeeze = False, sharex = True)
    plt.tight_layout(pad = pad, w_pad = w_pad, h_pad = h_pad)
    f.subplots_adjust(hspace=0.05)
    normalized_ylims = [0.97, 1.03]
    fill_lims = [0.99, 1.01]
    x_lims = [np.inf, -np.inf]
    solar_cell_area = (125 ** 2.0 - 15 ** 2.0) * 10.0 ** -6.0 #solar cell area in m^2 
    hard_x_lims = [0.0002 * 1000.0, 25.0 * 1000.0] # in nanoamps
    hard_x_lims = [hard_x_lims[0] / (solar_cell_area * 100.0 ** 2.0), hard_x_lims[1] / (solar_cell_area * 100.0 ** 2.0)] #in microamps / square centimeter
    hard_xticks = c.flattenListOfLists([[i * 10.0 ** power for i in range(10)] for power in range(-4, 2)])
    hard_y_lims = [0.91, 1.09]
    hard_yticks = [0.92, 0.96, 1.0, 1.04, 1.08]
    if fiber_save_strs_set in ['default','DEFAULT','Default']:
        fiber_save_strs_set = fiber_strs_set[:] 
    for led_num in range(len(LEDIDs)):
        LEDid = LEDIDs[led_num ]
        color = LED_point_colors[led_num]
        err_color = LED_err_colors[led_num]
        fiber_strs = fiber_save_strs_set[led_num]
        nFibers = len(fiber_strs) 
        ammeter_A_to_Vs = ammeter_A_to_Vs_set[led_num ]
        ammeter_strs = [str(A_to_V) + 'uAperV' for A_to_V in ammeter_A_to_Vs] 
        LI_sens_strs =  LI_sens_strs_set[led_num ]
        LI_sensis = [int(sens) if len(sens.split('_')) == 1 else int(sens.split('_')[0]) / int(sens.split('_')[1]) for sens in LI_sens_strs]
        target_dir =  target_dirs[led_num ]
        target_files = [file_funct(LEDid, fiber_strs[i], LI_sens_strs[i], ammeter_strs[i], wave_info_str, filt_info_str) for i in range(len(fiber_strs))]
        byHandFileObjects = [open(target_dir + target_file, 'r') for target_file in target_files]
        
        app_voltages_dict = {}
        SC_voltages_dict = {}
        PD_voltages_dict = {} 
        SC_voltages_med_dict = {} 
        PD_voltages_med_dict = {} 
        SC_voltages_ste_dict = {} 
        PD_voltages_ste_dict = {} 
        ratios_dict = {} 
        ratio_stds_dict = {}  
        SC_currents_med_dict = {} 
        SC_currents_ste_dict = {} 

        app_voltages_trunc_dict = {}
        SC_voltages_st_dict = {}
        PD_voltages_st_dict = {} 
        SC_voltages_med_st_dict = {} 
        PD_voltages_med_st_dict = {} 
        SC_voltages_ste_st_dict = {} 
        PD_voltages_ste_st_dict = {} 
        ratios_st_dict = {} 
        ratio_stds_st_dict = {}  
        SC_currents_med_st_dict = {} 
        SC_currents_ste_st_dict = {}

        for i in range(len(fiber_strs)): 
            fiber_str = fiber_strs[i] 
            print ('Reading in and managing data for fiber ' + str(fiber_str) + '...') 
            ammeter_A_to_V = ammeter_A_to_Vs[i] 
            byHandFileObject = byHandFileObjects[i] 
            LI_sens = LI_sensis[i]
            LI_sens_str = LI_sens_strs[i] 
            lines = byHandFileObject.readlines()
            if not(ignore_line_char is None):
                lines = [line for line in lines if not(line[0:len(ignore_line_char)] == ignore_line_char)]
            parsed_lines = [line[:-1].split(' ') for line in lines]
            app_voltages = [float(line[0]) for line in parsed_lines[1:]]
            SC_voltages = [[float(elem) for elem in line[2][1:-1].split(',')] for line in parsed_lines[1:]]
            PD_voltages = [[float(elem) for elem in line[1][1:-1].split(',')] for line in parsed_lines[1:]]
            truncated_app_voltages = []
            stacked_SC_voltages = [] 
            stacked_PD_voltages = [] 
            for j in range(len(app_voltages)): 
                app_voltage = app_voltages[j]
                if not(app_voltage in truncated_app_voltages): 
                    truncated_app_voltages = truncated_app_voltages + [app_voltage] 
                    stacked_SC_voltages = stacked_SC_voltages + [SC_voltages[j]]
                    stacked_PD_voltages = stacked_PD_voltages + [PD_voltages[j]] 
                else: 
                    app_volt_index = truncated_app_voltages.index(app_voltage) 
                    stacked_SC_voltages[app_volt_index] = stacked_SC_voltages[app_volt_index] + SC_voltages[j] 
                    stacked_PD_voltages[app_volt_index] = stacked_PD_voltages[app_volt_index] + PD_voltages[j] 
            SC_voltages_med = [np.median(single_meas) for single_meas in SC_voltages]
            PD_voltages_med = [np.median(single_meas) for single_meas in PD_voltages]
            SC_voltages_std = [np.std(single_meas) / np.sqrt(len(single_meas)) for single_meas in SC_voltages]
            PD_voltages_std = [np.std(single_meas) / np.sqrt(len(single_meas)) for single_meas in PD_voltages]
            stacked_SC_voltages_med = [np.median(single_meas) for single_meas in stacked_SC_voltages]
            stacked_PD_voltages_med = [np.median(single_meas) for single_meas in stacked_PD_voltages]
            stacked_SC_voltages_std = [np.std(single_meas) / np.sqrt(len(single_meas)) for single_meas in stacked_SC_voltages]
            stacked_PD_voltages_std = [np.std(single_meas) / np.sqrt(len(single_meas)) for single_meas in stacked_PD_voltages]
            ratios = np.array(SC_voltages_med) / np.array(PD_voltages_med) 
            stacked_ratios = np.array(stacked_SC_voltages_med) / np.array(stacked_PD_voltages_med) 
            ratio_stds = ((np.array(SC_voltages_std) / np.array(PD_voltages_med) ) ** 2.0 + ((np.array(PD_voltages_std) * np.array(SC_voltages_med) ) / (np.array(PD_voltages_med) ** 2.0 )) ** 2.0 ) ** 0.5
            stacked_ratio_stds = ((np.array(stacked_SC_voltages_std) / np.array(stacked_PD_voltages_med) ) ** 2.0 + ((np.array(stacked_PD_voltages_std) * np.array(stacked_SC_voltages_med) ) / (np.array(stacked_PD_voltages_med) ** 2.0 )) ** 2.0 ) ** 0.5 
            V_to_A = 10.0 ** (-6.0) * ammeter_A_to_V * (LI_sens * 10.0 ** -3.0 / 10)
            print ('[ammeter_A_to_V, LI_sens, V_to_A] = ' + str([ammeter_A_to_V, LI_sens, V_to_A])) 
            SC_currents_med = np.array(SC_voltages_med) * V_to_A
            SC_currents_std = np.array(SC_voltages_std) * V_to_A  
            stacked_SC_currents_med = np.array(stacked_SC_voltages_med) * V_to_A
            stacked_SC_currents_std = np.array(stacked_SC_voltages_std) * V_to_A  
            app_voltages_dict[LI_sens_str] = app_voltages[:] 
            SC_voltages_dict[LI_sens_str] = SC_voltages[:]
            PD_voltages_dict[LI_sens_str] = PD_voltages[:]
            SC_voltages_med_dict[LI_sens_str] = SC_voltages_med[:] 
            PD_voltages_med_dict[LI_sens_str] = PD_voltages_med[:] 
            SC_voltages_ste_dict[LI_sens_str] = SC_voltages_std[:] 
            PD_voltages_ste_dict[LI_sens_str] = PD_voltages_std[:]
            ratios_dict[LI_sens_str] = ratios[:]
            ratio_stds_dict[LI_sens_str] = ratio_stds[:]   
            SC_currents_med_dict[LI_sens_str] = SC_currents_med[:] 
            SC_currents_ste_dict[LI_sens_str] = SC_currents_std[:]
            app_voltages_trunc_dict[LI_sens_str] = truncated_app_voltages[:]
            SC_voltages_st_dict[LI_sens_str] = stacked_SC_voltages[:]
            PD_voltages_st_dict[LI_sens_str] = stacked_PD_voltages[:]
            SC_voltages_med_st_dict[LI_sens_str] = stacked_SC_voltages_med[:] 
            PD_voltages_med_st_dict[LI_sens_str] = stacked_PD_voltages_med[:] 
            SC_voltages_ste_st_dict[LI_sens_str] = stacked_SC_voltages_std[:] 
            PD_voltages_ste_st_dict[LI_sens_str] = stacked_PD_voltages_std[:]
            ratios_st_dict[LI_sens_str] = stacked_ratios[:]
            ratio_stds_st_dict[LI_sens_str] = stacked_ratio_stds[:]   
            SC_currents_med_st_dict[LI_sens_str] = stacked_SC_currents_med[:] 
            SC_currents_ste_st_dict[LI_sens_str] = stacked_SC_currents_std[:]
        
        scats = [ 0 for fiber_str in fiber_strs ]
        print ('np.shape(axarr) = ' + str(np.shape(axarr)))
        weights = np.array(c.flattenListOfLists([ratio_stds_dict[LI_sens_str] for LI_sens_str in LI_sens_strs])) ** (-2.0)
        print('weights = ' + str(weights)) 
        all_ratios = np.array(c.flattenListOfLists([ratios_dict[LI_sens_str] for LI_sens_str in LI_sens_strs])) 
        print ('[np.shape(weights[0]), np.shape(weights[1]), np.shape(all_ratios[0]), np.shape(all_ratios[1])] = ' + str([np.shape(weights[0]), np.shape(weights[1]), np.shape(all_ratios[0]), np.shape(all_ratios[1])])) 
        normalization = np.sum(all_ratios * weights) / np.sum(weights)
        print ('normalization = ' + str(normalization))
        ax1 = axarr[led_num][0]
        #ax2 = ax1.twinx()
        ylims = normalized_ylims[:]
        for i in range(len(LI_sens_strs)): 
            fiber_str = fiber_strs[i]
            LI_sens_str = LI_sens_strs[i] 
            print ('Plotting data for Lock-in sensitivity ' + str(LI_sens_str) + '...')
            marker = markers[i] 
            app_voltages = app_voltages_dict[LI_sens_str] [:-1]
            if stacked: 
                app_voltages = app_voltages_dict[LI_sens_str] [:-1]
                ratios = ratios_st_dict[LI_sens_str]  [:-1]
                ratio_stds = ratio_stds_st_dict[LI_sens_str]  [:-1]
                SC_currents_std = SC_currents_ste_st_dict[LI_sens_str]  [:-1]
                SC_currents_med = SC_currents_med_st_dict[LI_sens_str]  [:-1]
            else: 
                app_voltages = app_voltages_trunc_dict[LI_sens_str] [:-1] 
                ratios = ratios_dict[LI_sens_str]  [:-1]
                ratio_stds = ratio_stds_dict[LI_sens_str]  [:-1]
                SC_currents_std = SC_currents_ste_dict[LI_sens_str]  [:-1]
                SC_currents_med = SC_currents_med_dict[LI_sens_str]  [:-1]
            SC_current_densities = SC_currents_med / (solar_cell_area * 100.0 ** 2.0)
            SC_current_density_stds = SC_currents_std / (solar_cell_area * 100.0 ** 2.0)
            scats[i] = ax1.scatter([current * 10.0 ** 9.0 for current in SC_current_densities], np.array(ratios) / normalization, c = color, marker = marker)
            ylims = [min(ylims[0], min((np.array(ratios) - np.array(ratio_stds))/ normalization) * 0.95), max(ylims[1], max((np.array(ratios) + np.array(ratio_stds))/ normalization) * 1.05)]
            ax1.errorbar([current * 10.0 ** 9.0 for current in SC_current_densities], np.array(ratios) / normalization, xerr = [current * 10.0 ** 9.0 for current in SC_current_density_stds],  yerr = np.array(ratio_stds) / normalization, fmt = 'none', color = err_color) 
            #upOneSigRatios = np.array(ratios) + np.array(ratio_stds)
            #downOneSigRatios = np.array(ratios) - np.array(ratio_stds)
            #delta_ratios = (max(upOneSigRatios ) - min(downOneSigRatios) )  
            #y_lims = [min(downOneSigRatios) - delta_ratios * 0.05, max(upOneSigRatios) + delta_ratios * 0.05]

        ylims = hard_y_lims
        ax1.set_ylim(ylims)
        ax1.set_yticks(hard_yticks)
        #ax1.set_xlabel('Vpp of 10Hz TTL Driving LED (mV)')
        new_x_lims = ax1.get_xlim()
        x_lims = [min(new_x_lims[0], x_lims[0]), max(new_x_lims[1], x_lims[1])]
        x_lims = [new_x_lims[0] * 0.9, new_x_lims[1] * 1.1] 
        print ('x_lims = ' + str(x_lims)) 
        #ax2.set_xlim(x_lims) 
        ax1.set_xscale('log')
        #ax1.set_xticks(hard_xticks) 
        #ax2.set_xscale('log') 
        ax1.set_ylabel('(SC LI) / (PD LI) \n (normalized)', fontsize = 12.0) 
        #ax1.set_xlabel(r'SC current while LED on ($\mu A$)')
        #ax1.set_title_label('right') 
        #ax1.set_title('Linearity: LED ' + LEDid)
        #ax2.set_ylabel('Linearity: LED ' + LEDid)
        #'LED ' + LEDid,
        print ('LEDid = ' + str(LEDid))
        print ('x_lims = ' + str(hard_x_lims)) 
        print ('[10.0 ** ((np.log10(hard_x_lims[1]) + np.log10(hard_x_lims[0])) / 2.0), (ylims[1] + ylims[0]) / 2.0] = ' + str([10.0 ** ((np.log10(hard_x_lims[1]) + np.log10(hard_x_lims[0])) / 2.0), (ylims[1] + ylims[0]) / 2.0])) 
        ax1.text(10.0 ** ((np.log10(hard_x_lims[1]) + np.log10(hard_x_lims[0])) / 2.0), (ylims[1] + ylims[0]) / 2.0 + (ylims[1] - ylims[0]) * 0.3, 'LED ' + LEDid, ) 
        ax1.legend(scats, ['LI Sens ' + LI_sens_str + 'mV' for LI_sens_str in LI_sens_strs]) 
        ax1.fill_between(hard_x_lims, [fill_lims[0], fill_lims[0]], [fill_lims[1], fill_lims[1]], color = 'grey', alpha = 0.5)
        ax1.set_xlim(hard_x_lims)

    axarr[-1][0].set_xlabel(r'SC current density (nA$/$cm$^2$), measured by lock-in amplifier', fontsize = 12.0)
    #plt.tight_layout() 
    if save:
        plt.savefig(save_dir + save_plot_name)
    if show:
        plt.show() 
