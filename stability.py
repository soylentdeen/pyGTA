import serial
import time
import numpy
import srs510
import monochromator
import motion_control
from testing_modes import *

class InstrumentSuite(object):
    def __init__(self, sim, parameters):
        self.mono = monochromator.Monochromator(sim, parameters)
        self.srs = srs510.SRS510(sim, parameters)
        self.motor = motion_control.Motion_Controller(parameters)
        #success = self.motor.go_home()
        #if success != 1:
        #    print 'Error!  Motor did not go home!'

    def close(self):
        self.mono.close()
        self.srs.close()
        self.motor.close()


ver = 0.3
simulate = False  # Take live data

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

inst = InstrumentSuite(simulate, parameters)

data_dir = parameters['DATA_DIR']
print 'Data will be stored in : ', data_dir
outfile = raw_input('Enter the output data file :')

data_file = data_dir+outfile
for f in [data_file]:
    with open(f, 'w') as file:
        file.write(time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.localtime()))
        file.write('\n')


files = [data_file]

stability_test = Stability(inst, parameters)
stability_test.setup()
stability_test.run(files)

inst.close()

#print asdf
