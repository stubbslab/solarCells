import numpy as np

if __name__=="__main__":
    freqs = [25, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500]
    data_files = ['PD_M625F2LED_ChopFreq' + str(freq) + 'Hz_20210615.txt' for freq in freqs]
    data_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/solarCell/timeResponse/20210615/'

    for data_file in data_files:
        open_file = open(data_dir + data_file, 'r')
        lines = open_file.readlines()
        open_file.close()
        print ('lines = ' + str(lines))
        trimmed_lines = [line for line in lines if not(line in ['', ' \n', '\n'])]
        print ('trimmed_lines = ' + str(trimmed_lines))
        write_file = open(data_dir + data_file, 'w')
        write_file.writelines(trimmed_lines)
        print ('Done with file  ' + str(data_file))
