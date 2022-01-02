import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize
import cantrips as can
import cmath

def readInCurrentAndVoltages(data_file_name, data_dir = '', n_ignore = 1):

    all_data = can.readInColumnsToList(data_file_name, file_dir = data_dir, n_ignore = n_ignore, convert_to_float = 1, delimiter = ',')

    return all_data


def plotCurrentAndVoltages(data_file_name, save_file_name, save = 1, show = 0, data_dir = '', time_sep = 1, figsize = [10, 9], colors = ['r', 'b', 'k'], n_burn_in = 0):
    plt.rc('font', family='serif')
    plt.rc('text', usetex=True)
    dark_currents, dark_voltages = readInCurrentAndVoltages(data_file_name, data_dir = data_dir )
    delta_times = [i * time_sep for i in range(len(dark_currents[n_burn_in:]))]
    dark_currents = [current * 10.0 ** 6.0 for current in dark_currents[n_burn_in:]]
    dark_voltages = [voltage * 10.0 ** 3.0 for voltage in dark_voltages[n_burn_in:]]
    resistances = np.abs(np.array(dark_voltages) / np.array(dark_currents)) 
    f, axarr = plt.subplots(3,1, figsize = figsize, sharex = True)
    xlabel = r'$\Delta t$ (s)'
    ylabels = [r'Dark current ($\mu$A)', r'Dark voltage (mV)', r'Inferred resistance (k$\Omega$)']
    A_scat, v_scat, R_scat = [axarr[0].plot(delta_times, dark_currents, c = colors[0])[0], axarr[1].plot(delta_times, dark_voltages, c = colors[1])[0], axarr[2].plot(delta_times, resistances, c = colors[2])[0]]
    axarr[1].set_xlabel(xlabel)
    [axarr[i].set_ylabel(ylabels[i]) for i in range(len(ylabels))]
    if save:
        plt.savefig(data_dir + save_file_name )
    if show:
        plt.show()
    plt.close('all')
    return 1

if __name__ == "__main__":
    data_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/solarCell/darkCurrent/20210503/'
    data_file = 'SC_dark_current_on_B2987A_20210503_N70000.txt'
    save_img_file = 'SC_dark_current_on_B2987A_20210503_N70000.pdf'
    time_sep = 0.75
    n_burn_in = 1000
    plotCurrentAndVoltages(data_file, save_img_file, save = 1, data_dir = data_dir, time_sep = time_sep, n_burn_in = n_burn_in )
