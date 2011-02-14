import socket
import sys
import numpy
import time
import scipy.optimize
import scipy.interpolate
import scipy.integrate
import matplotlib.pyplot as pyplot
import GT_util

class BlazeEfficiencyReduce( object ):
    def __init__( self ):
        """ Sets up the blaze efficiency mode"""
        self.wl = []
        self.beta = []


    def get_Background(self, bg_file):
        x = []
        y = []
        with open(bg_file, 'r') as f:
            bg = f.read().split('\n')
            for line in bg:
                if (len(line) > 0):
                    l = line.split(', ')
                    if len(l) == 3:
                        x.append(float(l[0]))
                        y.append(float(l[1]))

        return [[x], [y]]

    def get_PSF(self, psf_file):
        xpts = []
        ypts = []
        dypts = []
        wls = []
        strength = []

        with open(psf_file, 'r') as f:
            psf = f.read().split('\n')
            for line in psf:
                if (len(line) > 0):
                    l = line.split(', ')
                    if len(l) == 3:
                        if float(l[0]) < 100:
                            xpts.append(float(l[0]))
                            ypts.append(float(l[1]))
                        else:
                            wls.append(float(l[0]))
                            strength.append(float(l[1]))
        retval = [[xpts, ypts], wls, strength]
        return retval

    def sweep(self, psf, scale_factors, wls, raw_file):
        x = []
        y = []
        raw_wl = None
        eff = []
        deff = []

        with open(raw_file, 'r') as f:
            raw = f.read().split('\r\n')
            for line in raw:
                if (len(line) > 0):
                    l = line.split(', ')
                    #print l
                    if len(l) == 1:
                        if (raw_wl):
                            observations = scipy.integrate.simps(y, x)
                            bm = scipy.where(numpy.array(wls) == raw_wl)
                            reference = scipy.integrate.simps(psf[1], psf[0])*scale_factors[bm[0]]/scale_factors[0]
                            eff.append( (observations*1.0)/(reference*1.0))
                            deff.append(0.02)
                            print scale_factors[bm[0]]
                            print raw_wl
                            print '    ---------------- '
                            
                        raw_wl = float(l[0])
                        x = []
                        y = []
                    elif len(l) == 3:
                        if ((len(x) > 0) and (x[-1] == float(l[0]))):
                            old = y[-1]
                            new = (old+float(l[1]))/2.0
                            y[-1] = new
                        else:
                            x.append(float(l[0]))
                            y.append(float(l[1]))
                    #raw_input()

        observations = scipy.integrate.simps(y, x)
        bm = scipy.where(numpy.array(wls) == raw_wl)
        reference = scipy.integrate.simps(psf[1], psf[0])*scale_factors[bm[0]]/scale_factors[0]
        eff.append( (observations*1.0)/(reference*1.0))
        deff.append( 0.02)
        print scale_factors[bm[0]]
        print raw_wl
        return [wls, eff, deff]

        '''
        with open(raw_file, 'a') as f:
            f.write(str(wl)+'\n')
            for angle, reading, dreading in zip(x, zip(*y)[0], zip(*y)[1]):
                f.write(str(angle)+', '+str(reading)+', '+str(dreading)+'\n')
                    
        return retval, x, zip(*y)[0], zip(*y)[1]
        '''

    def run(self, files):

        # Sets up the data files
        raw_file = files[0]
        bg_file = files[1]
        psf_file = files[2]
        data_file = files[3]

        # Measures the background
        background = self.get_Background(bg_file)

        psfs = self.get_PSF(psf_file)

        psf = psfs[0]
        wls = psfs[1]
        scale_factors = psfs[2]

        effs = self.sweep(psf, scale_factors, wls, raw_file)
        
        print len(effs[0])
        print effs[0]
        print len(effs[1])
        print effs[1]
        print len(effs[2])
        print effs[2]
        GT_util.save_figure([effs[0]], [effs[1]], [effs[2]], ['54829'], 'Blaze Efficiency', 'Wavelength (nm)', 'Efficiency', data_file+'.png')


test = BlazeEfficiencyReduce()

raw_file = 'kosi_54829_100710_2.out.raw'
bg_file = 'kosi_54829_100710_2.out.bkgnd'
psf_file = 'kosi_54829_100710_2.out.psf'
data_file = 'kosi_54829_100710_2.out'

directory = '/home/deen/Data/Instrumentation/IGRINS/VPH_Testing/Efficiency_Testing/K_Band_VPH/kosi_54829/100710/'

datafiles = [directory+raw_file, directory+bg_file, directory+psf_file, directory+data_file]

test.run(datafiles)
