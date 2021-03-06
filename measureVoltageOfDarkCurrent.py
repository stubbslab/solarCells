import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize
import cantrips as c
import cmath 

def readInVoltages(file_names = 'default', data_dir = '', freqs = [] ):
    if file_names is 'default':
        Vin_files = ['Vin_freq_' + str(freq) + 'Hz.txt'  for freq in freqs]
        Vout_files = ['Vout_freq_' + str(freq) + 'Hz.txt'  for freq in freqs]

    Vin_data = [c.readInColumnsToList(file_name, file_dir = data_dir, n_ignore = 1, convert_to_float = 1) for file_name in Vin_files]
    Vout_data = [c.readInColumnsToList(file_name, file_dir = data_dir, n_ignore = 1, convert_to_float = 1) for file_name in Vout_files]

    return [[data[1] for data in Vin_data], [data[1] for data in Vout_data]] 
        

def plotVoltages(save = 1, show = 0, target_dir = '', save_name = 'Vout_over_Vin_for_cap.png', data_dir = '', names = 'default'):
    if names is 'default':
         freqs = ['AA','AB','AC','AD','AE','AF','AG','AI','AJ','AK','AL']
    plt.rc('font', family='serif')
    plt.rc('text', usetex=True)
    dark_currents = readInVoltages(data_dir = data_dir, freqs = freqs )
    Vin_means = [np.mean(volts) for volts in Vin_volts]
    Vin_std = [np.std(volts) for volts in Vin_volts]
    Vout_means = [np.mean(volts) for volts in Vout_volts]
    Vout_std = [np.std(volts) for volts in Vout_volts]
    xlabel = 'Frequency of applied voltage (Hz)'
    ylabel = '$V_{\mathrm{out}}/V_{\mathrm{in}}$'
    scat = plt.scatter(freqs, np.array(Vout_means) / np.array(Vin_means))
    ZC = lambda f, C : complex(0.0, -(2.0 * cmath.pi * f * C) ** (-1.0) )
    ZR = lambda f, R: complex (R, 0.0)
    Z_R_par_C = lambda f, C, R: (1.0 / ZR (f, R) + 1.0 / ZC(f,C) ) ** -1.0
    C_ref = 10.0 * 10.0 ** -6.0 
    model_Vin_over_Vout = lambda f, R, C: abs(ZC(f, C_ref) / (ZC(f, C_ref) + Z_R_par_C(f, R, C)) ) 
    RC_fit = optimize.curve_fit(lambda freqs_in, R, C: [model_Vin_over_Vout(f, R, C) for f in freqs_in], freqs, np.array(Vout_means) / np.array(Vin_means), p0 = [100.0, 10.0 ** -6.0])
    print ('RC_fit = ' + str(RC_fit)) 
    line = plt.plot(freqs, [model_Vin_over_Vout(f, *[c.round_to_n(fit_param, 3) for fit_param in RC_fit[0]]) for f in freqs] )
    plt.legend(line, [r'$R_{\mathrm{sh}} = ' + str(c.round_to_n(RC_fit[0][1],3)) +  '\mathrm{\Omega}$' + '\n' + r'$C_{\mathrm{sh}} = ' + str(c.round_to_n(RC_fit[0][0] * 10.0 ** 6.0,3)) + ' \mathrm{\mu F}$']) 
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xscale('log')
    if save:
        plt.savefig(target_dir + save_name )
    if show:
        plt.show()
    plt.close('all')
    return 1 

if __name__ == "__main__":
    data_dir = '/Users/sasha/Documents/Harvard/physics/stubbs/solarCell/darkCurrentMeasurements/2019_06_12/' 
    target_dir = '/Users/sasha/Documents/Harvard/physics/stubbs/solarCell/darkCurrentMeasurements/2019_06_12/'
    file_name = 'unsoldered_dark_current_hists.png'
    plotVoltages(save = 1, target_dir = target_dir, save_name = file_name, data_dir = data_dir ) 
