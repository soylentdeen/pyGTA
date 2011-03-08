import monochromator
import srs510
import serial
import time
import numpy

configuration_text = open('grating_test.config.txt', 'r').read().split('\n')
parameters = dict()
for line in configuration_text:
    if(len(line) > 0):
        if line[0] != '#':
            l = line.split()
            try:
                parameters[l[0].upper()] = float(l[2])
                print 'Added float ', parameters[l[0].upper()], ' as entry for key ', l[0].upper()
            except:
                parameters[l[0].upper()] = l[2]
                print 'Added text variable ', parameters[l[0].upper()], ' as entry for key ', l[0].upper()

print "Grating Testing Program"
print "Cross your fingers!"

sim = False
srs = srs510.SRS510(sim, parameters)

a = srs.measure_const_SNR(10)

'''
mono = monochromator.Monochromator(False, parameters)

print mono.read_wavelength()

mono.goto_wavelength(1800)

mono.set_slit_width(500)

#raw_input()

'''
#mono.move_filter_wheel(4)

#raw_input("Filter wheel moved to position 4")

#mono.move_filter_wheel(6)

#raw_input("Filter wheel moved to position 6")
'''

mono.move_filter_wheel(1)

#raw_input("Filter wheel moved to position 6")
time.sleep(1)

mono.set_slit_width(100)

print mono.read_wavelength()

mono.goto_wavelength(1500)

time.sleep(10)

print mono.read_wavelength()
#'''
print "Ta Da!"

