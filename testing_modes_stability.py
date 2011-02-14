import socket
import sys
import numpy
import time
import scipy.optimize
import scipy.interpolate
import matplotlib.pyplot as pyplot
import IGRINS_util

class Grating( object ):
    def __init__(self, d, n):
        """ Initializes a grating object with grating constant d, and index of refraction n"""
        self.d = d   # grooves/mm
        self.n = n   # index of refraction

    def calc_beta(self, wavelength, alpha, order):
        """ wavelength = microns, alpha = degrees, returns degrees """
        #print wavelength
        #print alpha
        #print order
        #print self.d
        retval = numpy.arcsin(order*wavelength*self.d/1000.0 - numpy.sin(alpha/180.0*3.14159))
        #print retval
        return numpy.degrees(retval)

class Mode( object ):
    """A testing mode (super class)"""
    def run(self):
        """Runs the grating testing mode"""

    def setup(self):
        """Sets up the mode"""

    def reference(self):
        """Perorms a reference measurement"""

class Stability( Mode ):
    def __init__(self, file, inst_suite):
        print "Initialization"
        self.inst_suite = inst_suite
        self.wl = None
        setl.test_time = None

    def setup(self, wl):
        # do set up stuff here
        while (not self.wl):
            try:
                self.wl = float(raw_input('Enter the wavelength to test (in nm): '))
            except:
                print 'Error! Enter a valid wavelength!'
                self.wl = None

        while (not self.test_time):
            try:
                self.test_time = float(raw_input('How long do you want to run the test? (in seconds)'))
                if self.test_time < 0:
                    self.test_time = None
                    print 'Error!  You must enter a positive time!'
            except:
                print 'Error! Enter a valid test duration (in seconds)'
                self.test_time = None
        

    def run(self, files):
        data_file = files[0]
        raw_file = files[1]



class BlazeEfficiency( Mode ):
    def __init__( self, inst_suite, parameters ):
        """ Sets up the blaze efficiency mode"""
        self.inst_suite = inst_suite
        self.wl = []
        self.beta = []
        self.alpha = None
        self.nwl = len(self.wl)
        self.nbeta = len(self.beta)
        self.grating = None
        self.parameters = parameters

    def setup(self):
        """ Gets necessary parameters from the user to start the observation """
        wl_start = -1
        while (wl_start < 0):
            try:
                wl_start = float(raw_input('Enter Wavelength Start (in nm) :'))
            except:
                wl_start = -1.0
        wl_stop = -1.0
        while (wl_stop < wl_start):
            try:
                wl_stop = float(raw_input('Enter Wavelength Stop (in nm) :'))
            except:
                wl_stop = -1.0
        n_wl_pts = -1
        while (n_wl_pts <=0):
            try:
                n_wl_pts = int(raw_input('How many intermediate wavelengths do you wish to measure? :'))
            except:
                n_wl_pts = -1

        # Generates the list of wavelength points to test
        self.wl = numpy.linspace(wl_start, wl_stop, n_wl_pts)


        # For the purposes of this test, the blaze angle is also the incidence angle alpha
        while( not self.alpha ):
            try:
                self.alpha = float(raw_input('Enter Blaze angle (alpha) :'))
            except:
                self.alpha = None

        # Gets the grating parameters and creates the grating object
        a = True
        while a == True:
            try:
                gconst = float(raw_input('Enter the Grating constant (g/mm) : '))
                n = float(raw_input('Enter the index of refraction : '))
                self.grating = Grating(gconst, n)
                a = False
            except:
                print "Error! Please enter a sensible number!"
        # Calculates the beta angles for all the wavelengths to be tested
        self.beta = self.grating.calc_beta(self.wl/1000.0, self.alpha, 1)

        # Gets the user's desired Signal-to-Noise ratio
        a = True
        while a == True:
            try:
                self.SNR = float(raw_input('Enter your desired Signal-to-Noise :'))
                a = False
            except:
                print "Error! Please enter a sensible number!"

        # Sets the self.background variable to an empty list
        self.background = []

    def get_background(self, bg_file):
        x = numpy.linspace(self.parameters['REF_START_BETA'], self.parameters['REF_STOP_BETA'], self.parameters['N_REF_PTS'])
        y = []
        ch = raw_input('Please disconnect light source.  Press Enter to continue :')
        for b in x:
            print b
            self.inst_suite.motor.goto(b)
            time.sleep(self.parameters['SETTLE_TIME'])
            y.append(self.inst_suite.srs.measure_const_SNR(self.SNR))

        self.background = [x, y]

        IGRINS_util.save_figure([x], [zip(*y)[0]], [zip(*y)[1]], ['Background Signal'], 'Background', 'Angle (degrees)', 'Signal (V)', bg_file+'.png')
        with open(bg_file, 'a') as f:
            for angle, reading, dreading in zip(x, zip(*y)[0], zip(*y)[1]):
                f.write(str(angle)+', '+str(reading)+', '+str(dreading)+'\n')
        ch = raw_input('Replace light source.  Press Enter to continue :')
        

    def get_PSF(self, psf_file):
        x = numpy.linspace(self.parameters['REF_START_BETA'], self.parameters['REF_STOP_BETA'], self.parameters['N_REF_PTS'])
        xpts = []
        ypts = []
        dypts = []
        names = []
        retval = []

        self.inst_suite.mono.goto_wavelength(self.wl[0])
        y = []

        print "Starting PSF Scan"
        for b, y_bkgnd in zip(self.background[0], self.background[1]):
            self.inst_suite.motor.goto(b)
            time.sleep(self.parameters['SETTLE_TIME'])
            y.append(self.inst_suite.srs.measure_const_SNR(self.SNR) - y_bkgnd[0])
        retval.append([x, y])

        # Find the maximum value of y, go to that beta angle
        print y
        print self.background[0]
        order = numpy.array(zip(*y)[0]).argsort()
        self.inst_suite.motor.goto(self.background[0][order[-1]])
        bkgnd = self.background[1][order[-1]]
        for wl in self.wl:
            print "wavelength = "+str(wl)+"\n"
            self.inst_suite.mono.goto_wavelength(wl)
            time.sleep(self.parameters['SETTLE_TIME'])
            reading = self.inst_suite.srs.measure_const_SNR(self.SNR)
            retval.append([wl, reading[0]-bkgnd[0], (reading[1]**2.0+bkgnd[1]**2.0)**(0.5)])

        with open(psf_file, 'a') as f:
            for angle, reading in zip(*retval[0]):
                f.write(str(angle)+', '+str(reading[0])+', '+str(reading[1])+'\n')
            for wavelength, reading, dreading in retval[1:]:
                f.write(str(wavelength)+ ', '+str(reading)+', '+str(dreading)+'\n')
        IGRINS_util.save_figure([retval[0][0]], [zip(*retval[0][1])[0]], [zip(*retval[0][1])[1]], [str(self.wl[0])], 'Background Subtracted Reference PSF', 'Angle (degrees)', 'Signal (V)', psf_file+'.png')
        return retval

    def sweep(self, beta, psf, scale_factor, wl, raw_file):
        """ sweeps the detector through a range of angles centered on the predicted angle.
            Fits the observed data to the reference PSF """

        # creates array of beta angles to check
        x = numpy.linspace(beta+self.parameters['SWEEP_START_BETA'], beta+self.parameters['SWEEP_STOP_BETA'], self.parameters['N_SWEEP_PTS'])
        y = []
        for b in x:
            print b
            self.inst_suite.motor.goto(b)
            time.sleep(self.parameters['SETTLE_TIME'])
            y.append(self.inst_suite.srs.measure_const_SNR(self.SNR))

        #print psf
        #print scale_factor
        # Scales the PSF by the scale factor
        y_psf = [s[0] for s in psf[1]]*scale_factor
        y_obs = [s[0] for s in y]
        
        # beam and observations are interpolation objects.
        beam = scipy.interpolate.interpolate.interp1d(psf[0], y_psf)
        observations = scipy.interpolate.interpolate.interp1d(x, y_obs)

        # interpolates the observed data with a sampling of 10x
        xnew = numpy.linspace(beta-self.parameters['SWEEP_START_BETA'], beta+self.parameters['SWEEP_STOP_BETA'], self.parameters['N_SWEEP_PTS']*10.0)
        y_new = observations(xnew)

        # Initial guess for efficiency is just ratio of observaed maximum to PSF maximum
        pguess = [float(max(y_new))/float(max(y_psf)), 0.0]

        # Attempts to scale the observations to the PSF using two varibles
        # p[0] - efficiency (1.0 = 100%, 0.0 = 0%)
        # p[1] - x_slop (offsets the measured points in x)
        try:
            fitfunc = lambda p, beta: p[0]*beam(xnew-beta+p[1])
            errfunc = lambda p, beta, y: fitfunc(p, beta) - y
            p1, success = scipy.optimize.leastsq(errfunc, pguess, args = (beta, y_new))
            #p1 = [0.3, 0.1]
            print "Fitting Coefficients are : ", p1
        except:
            print "ERROR!  Fit did not converge!"
            p1 = [0.0, 0.0]

        with open(raw_file, 'a') as f:
            f.write(str(wl)+'\n')
            for angle, reading, dreading in zip(x, zip(*y)[0], zip(*y)[1]):
                f.write(str(angle)+', '+str(reading)+', '+str(dreading)+'\n')
                    
        return p1[0], x, zip(*y)[0], zip(*y)[1]

    def run(self, files):

        # Sets up the data files
        data_file = files[0]
        raw_file = files[1]
        bg_file = files[2]
        psf_file = files[3]
        xpts = []
        ypts = []
        dypts = []
        names = []
        with open(data_file, 'a') as f:
            f.write('Grating Constant : ' + str(self.grating.d)+'\n')
            f.write('Blaze/Alpha Angle : ' + str(self.alpha)+'\n')

        # Zeros out the readings variable
        self.readings = []

        # Ensure that all necessary variables are defined
        if (len(self.beta) <= 0) or (len(self.wl) <= 0):
            print "Error!  Set up the wavelength and alpha points!"
        else:

            # Measures the background
            self.get_background(bg_file)

            # Measures the reference PSF
            psfs = self.get_PSF(psf_file)

            # Gets ready to take grating data
            ch = raw_input('Place Grating in beam!')
            for  in zip(self.wl, self.beta, psfs[1:]):
                wl = int(pair[0])
                beta = float(pair[1])
                psf_scale_factor = pair[2]
                print "Pair[2] = ", pair[2]
                print "psf_scale_factor = ", psf_scale_factor
                print "psfs[0] = ", psfs[0]
                print "psfs[1][1] = ", psfs[1][1]
                print "Monochromator, please go to ", wl, " nanometers"
                self.inst_suite.mono.goto_wavelength(wl)
                time.sleep(10.0)
                while (self.inst_suite.mono.read_wavelength() != wl):
                    time.sleep(2.0)
                print "Motor Controller, please go to ", beta+self.alpha, " degrees"
                coeff, x, y, dy = self.sweep(beta+self.alpha, psfs[0], (psf_scale_factor[1]/psfs[1][1]), wl, raw_file)
                self.readings.append(coeff)
                xpts.append(x)
                ypts.append(y)
                dypts.append(dy)
                names.append(str(wl))
                with open(data_file, 'a') as f:
                    f.write(str(wl)+', '+str(beta)+', '+str(coeff)+'\n')
                print self.readings[-1]
            IGRINS_util.save_figure(xpts, ypts, dypts, names, 'Raw Data', 'Angle (degrees)', 'Signal (V)', raw_file+'.png')
        print "Done!"

class ScatteredLight( Mode ):
    """ """
    def __init__( self, inst_suite ):
        """Sets up the scatter-scan mode"""
        self.inst_suite = inst_suite
        self.wl = []
        self.beta = []
        self.nwl = len(self.wl)
        self.nbeta = len(self.beta)

    def setup(self):
        wl_start = -1.0
        while (wl_start < 0):
            try:
                wl_start = float(raw_input('Enter Wavelength Start (in nm) :'))
            except:
                wl_start = -1.0
        wl_stop = -1.0
        while (wl_stop < wl_start):
            try:
                wl_stop = float(raw_input('Enter Wavelength Stop (in nm) :'))
            except:
                wl_stop = -1.0
        n_wl_pts = -1
        while (n_wl_pts <=0):
            try:
                n_wl_pts = int(raw_input('Enter the number of datapoints at each angle :'))
            except:
                n_wl_pts = -1
        self.wl = numpy.linspace(wl_start, wl_stop, n_wl_pts)

        beta_start = None
        while (not (beta_start) ):
            try:
                beta_start = float(raw_input('Enter the starting Beta Angle : '))
            except:
                beta_start = None
        beta_stop = None
        while (not (beta_stop) ):
            try:
                beta_stop = float(raw_input('Enter the stopping Beta Angle : '))
            except:
                beta_stop = None
        n_beta_pts = -1
        while (n_beta_pts <= 0):
            try:
                n_beta_pts = int(raw_input('Enter the number of Beta points : '))
            except:
                n_beta_pts = -1
        self.beta = numpy.linspace(beta_start, beta_stop, n_beta_pts)

    def run(self):
        self.readings = []
        if (len(self.beta) <= 0) or (len(self.wl) <= 0):
            print "Error!  Set up the wavelength and beta points!"
        else:
            for x in self.wl:
                x = int(x)
                readings_wl_cut = []
                print "Goto ", x
                self.inst_suite.mono.goto_wavelength(x)
                time.sleep(5.0)
                while (self.inst_suite.mono.read_wavelength() != x):
                    print "waiting..."
                    time.sleep(2.0)
                for b in self.beta:
                    self.inst_suite.motor.goto(b)
                    time.sleep(1.0)
                    readings_wl_cut.append(self.inst_suite.srs.measure_const_SNR(100))
                    print readings_wl_cut[-1]
                self.readings.append(readings_wl_cut)
        return self.wl, self.beta, self.readings
