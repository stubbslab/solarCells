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
    wavelengths = [(row[0])  for row in data]
    stage_pos = [(row[1])  for row in data]
    PD_ill_times = [[ ( float(elem)) for elem in can.recursiveStrToListOfLists(row[2]) if elem != ignore_str] for row in data]
    PD_illum = [[ ( float(elem)) for elem in can.recursiveStrToListOfLists(row[3]) if elem != ignore_str] for row in data]
    PD_illum_bv = [[(float(elem)) for elem in can.recursiveStrToListOfLists(row[4]) if elem != ignore_str] for row in data]
    PD_dk_times = [[ ( float(elem)) for elem in can.recursiveStrToListOfLists(row[5]) if elem != ignore_str] for row in data]
    PD_dark = [[(float(elem)) for elem in can.recursiveStrToListOfLists(row[6]) if elem != ignore_str] for row in data]
    PD_dark_bv = [[(float(elem)) for elem in can.recursiveStrToListOfLists(row[7]) if elem != ignore_str] for row in data]
    SC_ill_times = [[ (float(elem)) for elem in can.recursiveStrToListOfLists(row[8]) if elem != ignore_str] for row in data]
    SC_illum = [[(float(elem)) for elem in can.recursiveStrToListOfLists(row[9]) if elem != ignore_str] for row in data]
    SC_illum_bv = [[(float(elem)) for elem in can.recursiveStrToListOfLists(row[10]) if elem != ignore_str] for row in data]
    SC_dk_times = [[ (float(elem)) for elem in can.recursiveStrToListOfLists(row[11]) if elem != ignore_str] for row in data]
    SC_dark = [[(float(elem)) for elem in can.recursiveStrToListOfLists(row[12]) if elem != ignore_str] for row in data]
    SC_dark_bv = [[(float(elem)) for elem in can.recursiveStrToListOfLists(row[13]) if elem != ignore_str] for row in data]
    #SC_dark_bv = [[(float(elem)) for elem in can.recursiveStrToListOfLists(row[7]) if elem != ignore_str] for row in data]
    SC_temps = [float(row[14]) for row in data]
    #SC_illum_bv = [[ (float(elem[1:-1])) for elem in can.recursiveStrToListOfLists(row[4])] for row in data]
    #SC_dark_bv = [[(float(elem[1:-1])) for elem in can.recursiveStrToListOfLists(row[5])] for row in data]
    #PD_illum_bv = [[(float(elem[1:-1])) for elem in can.recursiveStrToListOfLists(row[6])] for row in data]
    #PD_dark_bv = [[(float(elem[1:-1])) for elem in can.recursiveStrToListOfLists(row[7])] for row in data]

    #dial_positions, bias_voltages, SC_illum, SC_dark, PD_illum_bv, PD_dark_bv, SC_illum_bv, SC_dark_bv, PD_illum_bv, PD_dark_bv = can.safeSortOneListByAnother(bias_voltages, [zaber_positions, bias_voltages, SC_illum, SC_dark, PD_illum_bv, PD_dark_bv, SC_illum_bv, SC_dark_bv, PD_illum_bv, PD_dark_bv])
    #return [dial_positions, bias_voltages, SC_illum, SC_dark, PD_illum_bv, PD_dark_bv, SC_illum_bv, SC_dark_bv, PD_illum_bv, PD_dark_bv]
    #wavelengths, stage_pos, SC_ill_times, SC_illum, SC_illum_bv, SC_dk_times, SC_dark, SC_dark_bv, PD_ill_times, PD_illum, PD_illum_bv, PD_dk_times, PD_dark, PD_dark_bv, SC_temps = can.safeSortOneListByAnother(wavelengths, [wavelengths, stage_pos, SC_ill_times, SC_illum, SC_illum_bv, SC_dk_times, SC_dark, SC_dark_bv, PD_ill_times, PD_illum, PD_illum_bv, PD_dk_times, PD_dark, PD_dark_bv, SC_temps ])
    return [wavelengths, stage_pos, SC_ill_times, SC_illum, SC_illum_bv, SC_dk_times, SC_dark, SC_dark_bv, PD_ill_times, PD_illum, PD_illum_bv, PD_dk_times, PD_dark, PD_dark_bv, SC_temps]


def RC_funct(ts, freq, time_const, floor, asym_vpp, phi):
    #delta_t = np.mean([(ts[i + 1] + ts[i]) / 2.0 for i in range(0, len(ts) - 1)])
    #print ('[freq, time_const, floor, asym_vpp, phi] = ' + str([freq, time_const, floor, asym_vpp, phi]))
    ceiling = floor + asym_vpp
    charging_time = 1 / freq / 2
    average = (floor + ceiling) / 2.0
    lower_stop = (np.exp(charging_time / time_const) - np.exp(-charging_time / time_const)) ** -1.0 * (ceiling * (1 - np.exp(-charging_time / time_const)) + floor * (np.exp(charging_time / time_const) - 1) )
    upper_stop = lower_stop * np.exp(charging_time / time_const) - floor * (np.exp(charging_time / time_const) - 1)
    t0 = -phi / (freq * 2.0 * np.pi)
    delta_ts = ts - t0
    waveform_args = delta_ts * freq * 2.0 * np.pi
    high_or_low = np.where( (delta_ts - t0 ) // charging_time  % 2 == 0, 0, 1)
    local_delta_ts = (delta_ts - t0 ) % (1 / freq / 2)
    RC_wavefunct = np.where(high_or_low, upper_stop - (upper_stop - floor) * (1 - np.e ** (-local_delta_ts / time_const) ), (ceiling - lower_stop) * (1 - np.e ** (-local_delta_ts / time_const)) + lower_stop)
    #freq_shifted_ts = ( ts - 1/freq * phi / (2.0 * np.pi) )/ freq
    return RC_wavefunct

if __name__ == "__main__":
    dir_root = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/solarCell/timeResponse/'
    #data_files = [ stem + date_str + '.txt' for stem in ['SC_vs_PD_QE_with_B2987A_cellPosA_inten2_10_', 'SC_vs_PD_QE_with_B2987A_cellPosA_inten4_10_', 'SC_vs_PD_QE_with_B2987A_cellPosB_inten2_10_', 'SC_vs_PD_QE_with_B2987A_cellPosB_inten4_10_'] ]

    #freqs = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500]
    freqs = [200, 250, 300, 350]
    stems = ['PD_M625F2LED_ChopFreq' + str(freq) + 'Hz_' for freq in freqs]
    #stems = ['SC_QE_from_mono_SC_ED_']
    suffixes = ['' for freq in freqs]
    #suffixes = ['_singleStream1_straightSource']
    date_strs = ['20210615' for freq in freqs]
    #date_strs = ['20210607']
    extra_title_str = 'Cell ID ED, Steady Light'
    cell_ids = ['Cell ED - ' + date_strs[0] + ' ' + str(freq) + 'Hz' for freq in freqs]
    #cell_ids = ['Cell ED - ' + date_strs[0] + ' AM']
    shunt_resistances = [610 for freq in freqs]
    cap_guess = 10.0 ** -6.0
    #shunt_resistances = [610]
    plot_wavelengths = 0
    xlabel = r'$\Delta t$ (ms)'
    full_ylabel = r'$V_{out}$ (mV)'
    data_ylabel = r'$V_{out}$ (mV)'
    model_ylabel = r'$V_{Model}$ (mV)'
    resid_ylabel = r'$\Delta V_{out}$ (mV)'

    data_color, data_style = ['b', '--']
    model_color, model_style = ['r', ':']
    resid_color, resid_style = ['purple', '-']
    #if we don't plot wavelengths, we plot in times
    x_lims = [0.0, 0.055]

    data_files = [ stems[i] + date_strs[i] + suffixes[i] + '.txt' for i in range(len(freqs))]
    file_roots = [f[0:-len('.txt')] for f in data_files]
    save_file_name = 'PD_M625F2LED_ChopFreqs_' + '_'.join([str(freq) for freq in freqs]) +'.pdf'
    title_font_size = 24
    labelsize = 18
    ticksize = 16
    plt.rc('xtick',labelsize=ticksize)
    plt.rc('ytick',labelsize=ticksize)
    time_sep = 0.1

    data_dirs = [dir_root + date_str + '/' for date_str in date_strs]

    f, full_axarr = plt.subplots(len(data_files) * 3, 1, sharex = True, figsize = (15, 3 * len(data_files) * 2))

    for i in range(len(data_files)):
        data_ax = full_axarr[i * 3]
        model_ax = full_axarr[i * 3 + 1]
        resid_ax = full_axarr[i * 3 + 2]
        print ('Working on ' + str(i+1) + 'th data file of ' + str(len(data_files)))
        known_freq = freqs[i]
        Rsh = shunt_resistances[i]

        data_file = data_files[i]
        file_root = file_roots[i]
        data_dir = data_dirs[i]
        date_str = date_strs[i]

        #z_positions, bias_voltages, SC_illum, SC_dark, PD_illum_bv, PD_dark_bv, SC_illum_bv, SC_dark_bv, PD_illum_bv, PD_dark_bv = readInSCandPDDataFile(data_file, data_dir)

        time_strs, volt_strs = can.readInColumnsToList(data_file, data_dir, delimiter = ',', n_ignore = 1)
        delta_ts = [float(time) - float(time_strs[0]) for time in time_strs]
        volts = [float(volt) * 1000.0 for volt in volt_strs]
        N_points = len(volts)
        time_sep = np.mean([delta_ts[i+1] - delta_ts[i] for i in range(len(delta_ts) - 1)])
        volt_min, volt_max, volt_median, volt_std = [np.min(volts), np.max(volts), np.median(volts), np.std(volts)]
        ylims = [volt_min - (volt_max - volt_min) * 0.1, volt_max + (volt_max - volt_min) * 0.5]
        if x_lims == None:
            xticks = [can.round_to_n(tick, 2) for tick in np.linspace(delta_ts[0], delta_ts[-1], 11) ][1:-1]
        else:
            xticks = [can.round_to_n(tick, 2) for tick in np.linspace(*x_lims, 6) ][1:-1]
        xticklabels = [str(float(can.round_to_n(tick * 1000,3))) for tick in xticks ]
        print ('xticks = ' + str(xticks))
        print ('xticklabels = ' + str(xticklabels))
        xlims = [delta_ts[0] - time_sep, delta_ts[-1] + time_sep]
        data_fft = [scipy.fft.fftfreq(N_points, time_sep)[1:N_points//2], np.abs(scipy.fft.fft(volts)[0:N_points//2])[1:]]
        peak_freq = data_fft[0][np.argmax(data_fft[1])]
        init_freq, init_time_const, init_floor, init_vpp, init_phi = [peak_freq, Rsh * cap_guess, volt_min, volt_max - volt_min, np.pi * 1]
        init_guess = [init_freq, init_time_const, init_floor, init_vpp, init_phi]
        updated_guess = init_guess[:]
        #Do fit in two parts - first phase and frequency match, then fit amplitude and RC time constant
        amp_locked_RC_funct = lambda ts, f_false, phi: RC_funct(ts, f_false, init_time_const, init_floor, init_vpp, phi)
        freq_phase_matched_RC_curve = scipy.optimize.curve_fit(amp_locked_RC_funct, delta_ts, volts, bounds = [(0.0, -2.0 * np.pi * 0.1), (np.inf, 2.0 * np.pi * 1.1)], p0 = [init_guess[0], init_guess[4]])
        updated_guess[4] = freq_phase_matched_RC_curve [0][1]
        #updated_guess[4] = freq_phase_matched_RC_curve [0][1]
        false_fit_freq = freq_phase_matched_RC_curve [0][0]
        #We need to correct the data sampling based on our known frequency, relative to the innacurate frequency we measure.  That means rescaling the time.
        time_scaling = false_fit_freq / known_freq
        updated_guess[0] = known_freq
        print ('time_scaling = ' + str(time_scaling))
        delta_ts = np.array(delta_ts) * time_scaling
        scaled_xlims = [xlims[0] * time_scaling, xlims[1] * time_scaling]
        print ('freq_phase_matched_RC_curve[0] = ' + str(freq_phase_matched_RC_curve[0] ))
        print ('updated_guess = ' + str(updated_guess))
        fitted_RC_curve = scipy.optimize.curve_fit(RC_funct, delta_ts, volts, bounds = [(0.0, 0.0, -np.inf, 0.0, -2.0 * np.pi * 0.1), (np.inf, np.inf, np.inf, np.inf, 2.0 * np.pi * 1.1)], p0 = updated_guess)
        print ('fitted_RC_curve = ' + str(fitted_RC_curve ))
        fit_freq, fit_time_const, fit_floor, fit_vpp, fit_phi = fitted_RC_curve[0]
        fit_freq_err, fit_time_const_err, fit_floor_err, fit_vpp_err, fit_phi_err = [fitted_RC_curve[1][i,i] for i in range(5)]
        data_plot = data_ax.plot(delta_ts, volts, c = data_color, linestyle = data_style)[0]
        #axarr[i].plot(delta_ts, RC_funct(delta_ts, *init_guess), c = 'r')
        #axarr[i].plot(delta_ts, amp_locked_RC_funct(delta_ts, *(freq_phase_matched_RC_curve[0])), c = 'orange')
        #guess_plot = full_ax.plot(delta_ts, RC_funct(delta_ts, *(updated_guess)), c = 'g', linestyle = '--', alpha = 0.5)[0]
        fit_plot = model_ax.plot(delta_ts, RC_funct(delta_ts, *(fitted_RC_curve[0])), c = model_color, linestyle = model_style, alpha = 1.0)[0]
        resid_plot = resid_ax.plot(delta_ts, np.array(volts) - RC_funct(delta_ts, *(fitted_RC_curve[0])), c = resid_color, linestyle = resid_style, )[0]
        data_ax.legend([data_plot], [r'SC Chopped Data (expected $f=$' + str(known_freq) +'Hz)'], ncol = 2)
        model_ax.legend([fit_plot], [r'Best-fit RC model ($f=$' + str(can.round_to_n(fit_freq, 4)) + 'Hz , $RC=$' + str(can.round_to_n(fit_time_const * 10.0 ** 6, 2)) + r'$\mu$s'], ncol = 2)
        resid_ax.legend([resid_plot], [r'Residual (Data - fit) for expected $f=$' + str(known_freq) +'Hz)'])
        if i == len(data_files) - 1:
            print ('i = ' + str(i))
            print ('xticks = ' + str(xticks))
            print ('delta_ts = ' + str(delta_ts))
            resid_ax.set_xlabel(xlabel, fontsize = labelsize)
            resid_ax.set_xticks(xticks)
            resid_ax.set_xticklabels(xticklabels)
        else:
            resid_ax.set_xticks([])
        #full_ax.set_xticks([])
        data_ax.set_ylabel(data_ylabel, fontsize = labelsize)
        model_ax.set_ylabel(model_ylabel, fontsize = labelsize)
        resid_ax.set_ylabel(resid_ylabel, fontsize = labelsize)
        if not (x_lims == None):
            data_ax.set_xlim(x_lims)
            model_ax.set_xlim(x_lims)
            resid_ax.set_xlim(x_lims)
        else:
            data_ax.set_xlim(scaled_xlims)
            model_ax.set_xlim(scaled_xlims)
            resid_ax.set_xlim(scaled_xlims)
        data_ax.set_ylim(ylims)
        model_ax.set_ylim(ylims)
        #resid_ax.set_ylim(ylims)

    plt.savefig(data_dirs[-1] + save_file_name)
