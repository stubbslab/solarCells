import cantrips as can
import numpy as np
import scipy
import matplotlib.pyplot as plt

if __name__=="__main__":
    HG2_ref_wavelengths = [ 546.074, 578.013, 696.543, 706.722 , 714.704 , 727.294 , 738.398, 750.387 , 763.511, 772.376 , 794.818, 800.616, 811.531, 826.452, 842.465, 852.144, 866.794  ]
    #Identify lines by eye - open the fits file in ds9
    HG2_ref_pixel_guesses = [372,     411,     549,    561,       570,     585,      598,     612,       628,     639,      665,     672,    685,     703,     721,     735,      753, ]
    spec_lims = [450, 700]
    fit_funct = lambda xs, A, mu, sig, shift: A * np.exp(-(np.array(xs) - mu) ** 2.0 / sig ** 2.0) + shift
    fit_pix_width = 5
    mono_down_select_pix_width = 50

    data_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/skySpectrograph/data/ut20210628/'
    save_fig_name = 'MonochromatorSpectralSolution_20210628.pdf'
    HG2_ref_file = 'HG2_f23p4_2021_06_28_73.fits'
    mono_files = ['Mono450nm_f23p7_2021_06_28_78.fits', 'Mono500nm_f23p7_2021_06_28_83.fits', 'Mono550nm_f23p7_2021_06_28_89.fits', 'Mono600nm_f23p7_2021_06_28_94.fits', 'Mono650nm_f23p9_2021_06_28_98.fits',
                  'Mono700nm_f23p6_2021_06_28_104.fits', 'Mono750nm_f23p6_2021_06_28_105.fits', 'Mono800nm_f23p5_2021_06_28_113.fits', 'Mono850nm_f23p4_2021_06_28_119.fits', 'Mono900nm_f23p5_2021_06_28_124.fits',
                  'Mono950nm_f23p5_2021_06_28_129.fits', 'Mono1000nm_f23p4_2021_06_28_135.fits']
    ref_wavelengths = np.arange(450, 1001, 50)
    colors = ['b','r','g','orange','c','purple','firebrick','magenta','darkblue', 'slateblue', 'salmon','goldenrod']
    HG2_array, HG2_header = can.readInDataFromFitsFile(HG2_ref_file, data_dir)

    HG2_oneD_array = np.sum(HG2_array[spec_lims[0]:spec_lims[1], :], axis = 0)
    HG2_floor = np.median(HG2_oneD_array)
    zSub_HG2_oneD_array = HG2_oneD_array - HG2_floor
    pixels = range(len(HG2_oneD_array))
    n_pixels = len(pixels)
    fit_indeces = [ [max(0, ref_pix - fit_pix_width), min(n_pixels, ref_pix + fit_pix_width)] for ref_pix in HG2_ref_pixel_guesses ]
    fit_xs, fit_ys = [ [pixels[fit_index[0]:fit_index[1]] for fit_index in fit_indeces],  [HG2_oneD_array[fit_index[0]:fit_index[1]] for fit_index in fit_indeces] ]
    init_guesses = [ [max(fit_ys[i]) - HG2_floor, HG2_ref_pixel_guesses[i], 3.0, HG2_floor ] for i in range(len(HG2_ref_pixel_guesses)) ]

    best_fits = [ scipy.optimize.curve_fit(fit_funct, fit_xs[i], fit_ys[i], init_guesses[i]) for i in range(len(HG2_ref_pixel_guesses)) ]
    HG2_ref_pixels = [best_fit[0][1] for best_fit in best_fits]
    print ('HG2_ref_pixel_guesses = ' + str(HG2_ref_pixel_guesses))
    print ('HG2_ref_pixels = ' + str(HG2_ref_pixels))

    #plt.plot(pixels, HG2_oneD_array, c = 'k')
    #[plt.plot(fit_xs[i], fit_funct(fit_xs[i], *init_guesses[i]), c = 'r') for i in range(len(HG2_ref_pixel_guesses)) ]
    #[plt.plot(fit_xs[i], fit_funct(fit_xs[i], *best_fits[i][0]), c = 'g') for i in range(len(HG2_ref_pixel_guesses)) ]
    pix_to_wave = np.poly1d(np.polyfit(HG2_ref_pixels, HG2_ref_wavelengths, 2))
    wave_to_pix = np.poly1d(np.polyfit(HG2_ref_wavelengths, HG2_ref_pixels, 2))

    mono_data_reads = [can.readInDataFromFitsFile(mono_file, data_dir) for mono_file in mono_files]
    mono_arrays = [data_read[0] for data_read in mono_data_reads]
    mono_headers = [data_read[1] for data_read in mono_data_reads]
    mono_oneD_arrays = [np.sum(mono_array[spec_lims[0]:spec_lims[1], :], axis = 0) for mono_array in mono_arrays]
    mono_oneD_pixels = [ [int(max(wave_to_pix(ref_wavelengths[i]) - mono_down_select_pix_width, 0) ), int(min(wave_to_pix(ref_wavelengths[i]) + mono_down_select_pix_width, n_pixels))] for i in range(len(mono_oneD_arrays)) ]
    mono_oneD_arrays = [mono_oneD_arrays[i][mono_oneD_pixels[i][0]: mono_oneD_pixels[i][1]] for i in range(len(mono_oneD_arrays)) ]
    mono_floors = [np.median(mono_oneD_array) for mono_oneD_array in mono_oneD_arrays]
    zSub_mono_oneD_arrays = [mono_oneD_arrays[i] - mono_floors[i] for i in range(len(mono_floors))]
    f, axarr = plt.subplots(2,1, figsize = (15,8))

    axarr[0].plot(pix_to_wave(pixels), zSub_HG2_oneD_array  / np.max(zSub_HG2_oneD_array ))
    [axarr[1].plot(pix_to_wave(range(mono_oneD_pixels[i][0], mono_oneD_pixels[i][1])), zSub_mono_oneD_arrays[i] / np.max(zSub_mono_oneD_arrays[i]), c = colors[i]) for i in range(len(zSub_mono_oneD_arrays))]
    axarr[1].set_xlabel('Fitted wavelength (nm)', fontsize = 14)
    axarr[0].set_ylabel('Normalized Counts - HG2', fontsize = 14)
    axarr[1].set_ylabel('Normalized Counts - Mono', fontsize = 14)
    [axarr[0].axvline(wave, c = 'k', linestyle = '--', alpha = 0.5) for wave in HG2_ref_wavelengths]
    [axarr[1].axvline(ref_wavelengths[i], c = colors[i], linestyle = '--', alpha = 0.5) for i in range(len(ref_wavelengths))]
    mono_guesses = [ [np.max(zSub_mono_oneD_arrays[i]), wave_to_pix(ref_wavelengths[i]), 10.0, 0.0] for i in range(len(ref_wavelengths))]
    mono_fits = [scipy.optimize.curve_fit(fit_funct, range(mono_oneD_pixels[i][0], mono_oneD_pixels[i][1]), zSub_mono_oneD_arrays[i], mono_guesses[i]) for i in range(len(ref_wavelengths))]
    [axarr[1].plot(pix_to_wave(range(mono_oneD_pixels[i][0], mono_oneD_pixels[i][1])), fit_funct(range(mono_oneD_pixels[i][0], mono_oneD_pixels[i][1]), *mono_fits[i][0]) / np.max(zSub_mono_oneD_arrays[i]), c = 'k', alpha = 0.5) for i in range(len(zSub_mono_oneD_arrays))]

    mono_wave_centers = [pix_to_wave(fit[0][1] ) for fit in mono_fits]
    mono_wave_FWHM = [(pix_to_wave(fit[0][1] + fit[0][2]) - pix_to_wave(fit[0][1] - fit[0][2])) / 2 * 2.355 for fit in mono_fits]
    mono_delta_waves = np.array(mono_wave_centers) - np.array(ref_wavelengths)
    print ('mono_wave_FWHM = ' + str(mono_wave_FWHM))
    print ('mono_delta_waves = ' + str(mono_delta_waves))
    [axarr[1].text(pix_to_wave(mono_oneD_pixels[i][0]), 1.1, r'$\Delta \lambda$=' + str(can.round_to_n(mono_delta_waves[i], 3)) + 'nm', fontsize = 8, color = colors[i]) for i in range(len(mono_delta_waves))]
    [axarr[1].text(pix_to_wave(mono_oneD_pixels[i][0]), 0.95, r'FWHM $\lambda$=' + '\n' + str(can.round_to_n(mono_wave_FWHM[i], 3)) + 'nm', fontsize = 8, color = colors[i]) for i in range(len(mono_delta_waves))]
    axarr[1].set_ylim(-0.05, 1.2)
    axarr[0].set_ylim(-0.05, 1.2)
    axarr[0].set_xlim(400, 1050)
    axarr[1].set_xlim(400, 1050)
    axarr[0].set_xticks(np.arange(450, 1001, 50))
    axarr[1].set_xticks(np.arange(450, 1001, 50))
    plt.savefig(data_dir + save_fig_name)

    plt.show()
