import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize
import cantrips as c
import cmath
from scipy.stats.distributions import  t
import scipy.special

def readInVoltages(file_name = 'default', data_dir = '', freqs = [] ):
    if file_name is 'default':
        file_name = 'capVoltages_SCandPD_byHand__5_100VppSine_LISCSens50mV_LIFullSens50mV_seriesCap1uF.txt'
    cap_data = c.readInColumnsToList(file_name, file_dir = data_dir, n_ignore = 1, convert_to_float = 0, delimiter = ' ', ignore_line_char = '#')
    frequencies = [float(freq) for freq in cap_data[0]]
    Vin_data = [[float(volt) for volt in Vin_str[1:-1].split(',')] for Vin_str in cap_data[1]]
    Vout_data = [[float(volt) for volt in Vin_str[1:-1].split(',')] for Vin_str in cap_data[2]]
    return [frequencies[0:-1], Vin_data[0:-1], Vout_data[0:-1]]


def plotVoltageRatios(save = 1, show = 0, target_dir = '', save_name = 'Vout_over_Vin_for_cap.png', data_dir = '', freqs = 'default', file_name = 'default', max_inclusion = -1 ):
    #if freqs is 'default':
    #    freqs1 = [7 + freq for freq in list(range(0, 100, 10))]
    #    freqs2 = [7 + freq for freq in list(range(100, 200, 50))]
    #    freqs3 = [7 + freq for freq in list(range(200, 1000, 100))]
    #    freqs4 = [7 + freq for freq in list(range(1000, 2000, 500))]
    #    freqs5 = [7 + freq for freq in list(range(2000, 10000, 1000))]
    #    freqs6 = [7 + freq for freq in list(range(10000, 30000, 5000) + [30000])]
    #    #freqs8 = []
    #    #freqs9 = []
    #    #freqs10 = []
    #    freqs = freqs1 + freqs2 + freqs3 + freqs4 + freqs5 + freqs6
    #    print ('sorted(freqs) = ' + str(sorted(freqs)))
    #    print ('len(freqs) = ' + str(len(freqs)) )
    plt.rc('font', family='serif')
    plt.rc('text', usetex=True)
    plt.rc('xtick',labelsize=16)
    plt.rc('ytick',labelsize=16)
    freqs, Vin_volts, Vout_volts = readInVoltages(data_dir = data_dir, freqs = freqs, file_name = file_name )
    print ('[freqs, Vin_volts, Vout_volts] = ' + str([freqs, Vin_volts, Vout_volts]))
    print ('freqs = ' + str(freqs))
    #freqs = freqs[0:max_inclusion]
    print ('freqs = ' + str(freqs))
    #Vin_volts = Vin_volts[0:max_inclusion]
    #Vout_volts = Vout_volts[0:max_inclusion]
    Vin_means = [np.mean(volts) for volts in Vin_volts]
    Vin_ste = [np.std(volts) / np.sqrt(len(volts)) for volts in Vin_volts]
    Vout_means = [np.mean(volts) for volts in Vout_volts]
    Vout_ste = [np.std(volts) / np.sqrt(len(volts))  for volts in Vout_volts]
    xlabel = 'Frequency of applied voltage (Hz)'
    #xlabel = r'$V_{pp}$ (V)'
    ylabel = '$V_{\mathrm{out}}/V_{\mathrm{in}}$'
    title = 'Voltage ratio from cap-cap divider with lock-ins in lock-in mode'
    max_ratios = np.max(np.array(Vout_means) / np.array(Vin_means) )
    ratios = np.array(Vout_means) / np.array(Vin_means)
    ratio_sigs = np.sqrt((np.array(Vout_ste) / np.array(Vin_means)) ** 2.0 + (np.array(Vin_ste) * np.array(Vout_means)  /( np.array(Vin_means) ** 2.0)) ** 2.0)
    scat = plt.scatter(freqs, np.array(Vout_means) / np.array(Vin_means)) #  / max_ratios)
    ZC = lambda f, C : complex(0.0, -(2.0 * cmath.pi * f * C) ** (-1.0) )
    ZR = lambda f, R: complex (R, 0.0)
    Z_R_par_C = lambda f, R, C: (1.0 / ZR (f, R) + 1.0 / ZC(f,C) ) ** -1.0
    C_ref = 1.0 * 10.0 ** -6.0
    R_ser = 0.0
    print ('Here')
    model_Vin_over_Vout = lambda f, R_sh, C_sh:        abs(ZC(f, C_ref) / (ZC(f, C_ref) + Z_R_par_C(f, R_sh, C_sh) + ZR(f, R_ser)) )
    #model_Vin_over_Vout = lambda f, R_sh, C_sh, R_ser: abs(ZC(f, C_ref) / (ZC(f, C_ref) + Z_R_par_C(f, R_sh, C_sh) + ZR(f, R_ser)) )
    #RC_fit = optimize.curve_fit(lambda freqs_in, R_sh, C_sh, R_ser: [model_Vin_over_Vout(f, R_sh, C_sh, R_ser) for f in freqs_in], freqs, np.array(Vout_means) / np.array(Vin_means), p0 = [100.0, 10.0 ** -6.0, 1.0])
    RC_fit = optimize.curve_fit(lambda freqs_in, R_sh, C_sh: [model_Vin_over_Vout(f, R_sh, C_sh) for f in freqs_in], freqs, ratios, p0 = [100.0, 10.0 ** -6.0]) #, sigma = ratio_sigs)
    print ('RC_fit = ' + str(RC_fit))
    cov = RC_fit[1]
    uncertainty = np.sqrt(np.diag(cov))
    n_vars = 2    # number of fit parameters
    n_points = len(freqs) # number of data points
    dof = max(0, n_points - n_vars)
    alpha  = 1-scipy.special.erf(1) #1 sigma confidence value
    t_val = t.ppf(1.0-alpha/2., dof)
    print ('uncertainty = ' + str(uncertainty))
    scaled_uncertainty = [elem * t_val for elem in uncertainty]
    print ('scaled_uncertainty = ' + str(scaled_uncertainty))
    #print ('[c.round_to_n(fit_param, 3) for fit_param in RC_fit[0]] = ' + str([c.round_to_n(fit_param, 3) for fit_param in RC_fit[0]]))
    #print ('freqs[-1] = ' + str(freqs[-1]))
    #print ('model_Vin_over_Vout(freqs[-1], *[c.round_to_n(fit_param, 3) for fit_param in RC_fit[0]]) = '  + str(model_Vin_over_Vout(freqs[-1], *[c.round_to_n(fit_param, 3) for fit_param in RC_fit[0]])))
    #print ('[Z_R_par_C(freqs[-1], *RC_fit[0]), ZC(freqs[-1], C_ref)] = ' + str([Z_R_par_C(freqs[-1], *RC_fit[0]), ZC(freqs[-1], C_ref)]))

    #line = plt.plot(freqs, [model_Vin_over_Vout(f, *[c.round_to_n(fit_param, 3) for fit_param in RC_fit[0]], R_ser) for f in freqs] )

    #print ('[RC_fit[0], R_ser] = ' + str([RC_fit[0], R_ser]))
    line = plt.plot(freqs, [model_Vin_over_Vout(f, *[c.round_to_n(fit_param, 3) for fit_param in RC_fit[0]]) for f in freqs] )
    #plt.legend(line, [r'$R_{\mathrm{sh}} = ' + str(c.round_to_n(RC_fit[0][0],3)) +  '\mathrm{\Omega}$' + '\n' + r'$C_{\mathrm{sh}} = ' + str(c.round_to_n(RC_fit[0][1] * 10.0 ** 6.0,3)) + ' \mathrm{\mu F}$' + '\n' + r'$R_{\mathrm{ser}} = ' + str(c.round_to_n(RC_fit[0][2],3)) + ' \mathrm{\Omega}$'])
    #plt.legend(line, [r'$R_{\mathrm{sh}}$ = ' + str(c.round_to_n(RC_fit[0][0],4)) + r'$\pm$' + str(c.round_to_n(scaled_uncertainty[0],2)) +  r'$\mathrm{\Omega}$' + '\n' + r'$C_{\mathrm{sh}}$ = ' + str(c.round_to_n(RC_fit[0][1] * 10.0 ** 9.0,4))  + r'$\pm$' + str(c.round_to_n(scaled_uncertainty[1] * 10.0 ** 9.0,2)) + r'$\mathrm{nF}$'])
    plt.legend(line, [r'$R_{\mathrm{sh}}$ = ' + str(c.round_to_n(RC_fit[0][0],3)) + r'$\mathrm{\Omega}$' + '\n' + r'$C_{\mathrm{sh}}$ = ' + str(c.round_to_n(RC_fit[0][1] * 10.0 ** 6.0,3))  + r'$\mathrm{\mu F}$'], fontsize = 16)
    plt.xlabel(xlabel, fontsize = 18.0)
    plt.ylabel(ylabel, fontsize = 18.0)
    #plt.title(title)
    plt.xscale('log')
    #plt.yscale('log')
    if save:
        plt.tight_layout()
        plt.savefig(target_dir + save_name )
    if show:
        plt.tight_layout()
        plt.show()
    plt.close('all')
    return 1

def plotVoltages(save = 1, show = 0, target_dir = '', save_name = 'Vout_over_Vin_for_cap.png', data_dir = '', freqs = 'default', file_name = 'default'):
    plt.rc('font', family='serif')
    plt.rc('text', usetex=True)
    freqs, VPD_volts, VSC_volts = readInVoltages(data_dir = data_dir, freqs = freqs, file_name = file_name )
    amp_A_to_V = 10.0 ** 3.0
    LI_volt_scaling = 0.5 / 10.0
    V_to_A = LI_volt_scaling / amp_A_to_V
    VPD_currents = [[volt * V_to_A for volt in voltages] for voltages in VPD_volts]
    VSC_currents = [[volt * V_to_A for volt in voltages] for voltages in VSC_volts]
    VPD_means = [np.mean(volts) for volts in VPD_currents]
    VPD_ste = [np.std(volts)  / np.sqrt(len(volts)) for volts in VPD_currents]
    VSC_means = [np.mean(volts)  / np.sqrt(len(volts)) for volts in VSC_currents]
    VSC_ste = [np.std(volts) for volts in VSC_currents]
    xlabel = 'Frequency of applied voltage (Hz)'
    ylabel = r'Currents ($\mu$A)'
    title = 'Time response to photodiode and solar cell'
    scatPD = plt.scatter(freqs, np.array(VPD_means) * 10.0 ** 6.0)
    plt.errorbar(freqs, np.array(VPD_means) , yerr = VPD_ste, fmt = 'none')
    scatSC = plt.scatter(freqs, np.array(VSC_means) * 10.0 ** 6.0)
    plt.errorbar(freqs, np.array(VSC_means) , yerr = VSC_ste, fmt = 'none')
    plt.legend([scatPD, scatSC], ['Photodiode', 'Solar cell'])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.xscale('log')
    #plt.yscale('log')
    if save:
        plt.savefig(target_dir + save_name )
    if show:
        plt.show()
    plt.close('all')
    return 1



if __name__ == "__main__":
    target_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/solarCell/capMeasurements/20190829/'
    data_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/solarCell/capMeasurements/20190829/'
    file_name = 'Vdiode_and_VSC_for_cap.png'
    #plotVoltageRatios(save = 0, show = 1, target_dir = target_dir, save_name = file_name, data_dir = data_dir )
    plotVoltageRatios(save = 1, show = 1, target_dir = target_dir, save_name = file_name, data_dir = data_dir, max_inclusion = 24 )
