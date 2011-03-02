import socket
import sys
import numpy
import time
import scipy.optimize
import scipy.interpolate
import scipy.integrate
import matplotlib.pyplot as pyplot
import GT_util

class DispersiveElement( object ):
    def calc_beta(self):
        '''Returns the beta angle for a given measurement'''

class Grating( DispersiveElement ):
    def __init__(self, n, blaze, **kwargs):
        """ Initializes a grating object with grating constant d, and index of refraction n"""
        if ("d" in kwargs):
            self.d = kwargs["d"]   # grooves/mm
            self.sigma = 1000.0/self.d   #microns/groove
        elif ("sigma" in kwargs):
            self.sigma = kwargs["sigma"]
            self.d = 1000.0/self.sigma
        if ("name" in kwargs):
            self.name = kwargs["name"]
        else:
            self.name = 'Grating'
        self.n = n   # index of refraction
        self.blaze = blaze

    def calc_beta(self, wavelength, alpha, order):
        """ wavelength = microns, alpha = degrees, returns degrees """
        retval = numpy.arcsin(order*wavelength*self.d/1000.0 - numpy.sin(alpha/180.0*3.14159))
        return numpy.degrees(retval)

class Grism( DispersiveElement ):
    def __init__(self, n, blaze, delta, **kwargs):
        """ Initializes a grating object with grating constant d, and index of refraction n"""
        #print kwargs
        if ("d" in kwargs):
            self.d = kwargs["d"]   # grooves/mm
            self.sigma = 1000.0/self.d    # microns/groove
        elif ("sigma" in kwargs):
            self.sigma = kwargs["sigma"]   #microns/groove
            self.d = 1000.0/self.sigma
        if ("name" in kwargs):
            self.name = kwargs["name"]
        else:
            self.name = 'Grism'
        self.n = n   # index of refraction
        self.delta = delta   # Wedge angle of prism
        self.blaze = blaze   # angle of groove faces wrt grating surface

    def calc_beta(self, wavelength, alpha, order):
        """ wavelength = microns, alpha = degrees, returns degrees """
        retval = numpy.arcsin(order*wavelength/self.sigma - self.n*numpy.sin(self.delta/180.0*3.14159- numpy.arcsin(numpy.sin(alpha/180.0*3.14159/self.n))))
        return (numpy.degrees(retval)+self.delta)     # -1 b/c the grism mount is reversed


class Mode( object ):
    """A testing mode (super class)"""
    def run(self):
        """Runs the grating testing mode"""

    def setup(self):
        """Sets up the mode"""

class Stability( Mode ):
    def __init__(self, inst_suite, parameters):
        self.inst_suite = inst_suite
        self.wl = None
        self.test_time = None
        self.parameters = parameters
        self.mailer = None
        self.twitterer = GT_util.Twitterer()

    def setup(self):
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
                    print 'Error!  Unless entropy has reversed itself, you must enter a positive time!'
            except:
                print 'Error! Enter a valid test duration (in seconds)'
                self.test_time = None

	    #Sets up the email
        while (not self.mailer):
            try:
                emails = raw_input('Enter the email(s) you wish to notify when user intervention is required (separate with commas) :')
                self.mailer = GT_util.Mailer(self.parameters["MAILSERVER"], self.parameters["SENDER"], emails.split(','))
            except:
                print 'Error! Enter a valid email address!'

        # Gets the user's desired Signal-to-Noise ratio
        a = True
        while a == True:
            try:
                self.SNR = float(raw_input('Enter your desired Signal-to-Noise :'))
                a = False
            except:
                print "Error! Please enter a sensible number!"
                
        self.filter = -1
        while ( (self.filter <= 0) or (self.filter >= 7) ):
            try:
                self.filter = int(raw_input('Which Filter do you want to use?'))
            except:
                self.filter = -1
        self.inst_suite.mono.set_filter(self.filter)

        self.slit_width = -1
        while ( (self.slit_width < 50) or (self.slit_width > 2000) ):
            try:
                self.slit_width = int(raw_input("Enter the monochromator slit width :"))
            except:
                self.slit_width = -1

        self.inst_suite.mono.set_slit_width(self.slit_width)

    def get_Background(self):
        self.background = None

        # Set monochromator to correct wavelength
        self.inst_suite.mono.goto_wavelength(self.wl)

        # Open filter wheel
        self.inst_suite.mono.open_shutter()

        # Find middle of PSF
        print "Finding Middle of beam"
        angles = numpy.linspace(-0.5, 0.5, 11)
        y = []
        for angle in angles:
            txt_out = "Angle = "+str(angle)
            sys.stdout.write(txt_out)
            self.inst_suite.motor.goto(angle)
            time.sleep(self.parameters['SETTLE_TIME'])
            yval = self.inst_suite.srs.measure_const_SNR(10)
            y.append([yval[0], yval[1]])
            for i in numpy.arange(len(txt_out)):
                sys.stdout.write('\b')
            sys.stdout.flush()

        # Find the maximum value of y, go to that beta angle
        order = numpy.array(zip(*y)[0]).argsort()
        self.inst_suite.motor.goto(angles[order[-1]])
        
        # Close filter wheel
        self.inst_suite.mono.close_shutter()

        # Take background
        self.background = self.inst_suite.srs.measure_const_SNR(self.SNR)

    def get_TimeSequence(self, out_file):
        self.time_sequence = []
        self.elapsed_time = []
        
        # Open filter wheel
        self.inst_suite.mono.open_shutter()
        #self.inst_suite.mono.goto_wavelength(self.wl)
        time.sleep(self.parameters['SETTLE_TIME'])

        self.start_time = time.time()

        self.time_sequence.append(self.inst_suite.srs.measure_const_SNR(self.SNR))
        self.elapsed_time.append(time.time()-self.start_time)
        txt_out = ''

        n_milestones = 4

        milestones = numpy.linspace(0, self.test_time, n_milestones+1)
        last_milestone = 0

        wave = self.inst_suite.mono.read_wavelength()
        
        while self.elapsed_time[-1] < self.test_time:
            if (len(milestones) > last_milestone+1):
                if (self.elapsed_time[-1] > milestones[last_milestone+1]):
                    last_milestone += 1
                    fraction = str(((1.0*last_milestone)/(1.0*len(milestones)-1.0)*100.0))
                    mesg = 'Most momentus! I am '+fraction+'% done with the task! Huzzah!'
                    status = self.twitterer.tweet(mesg)
                    
            # Take reference measurement
            self.time_sequence.append(self.inst_suite.srs.measure_const_SNR(self.SNR))
            self.elapsed_time.append(time.time()-self.start_time)
            for i in numpy.arange(len(txt_out)):
                sys.stdout.write('\b')
            txt_out = str(self.test_time - self.elapsed_time[-1])+'  '+str(wave)
            sys.stdout.write(txt_out)

        with open(out_file, 'a') as f:
            f.write('Background : '+str(self.background[0])+' +/- '+str(self.background[1])+'\n')
            for seq, t in zip(self.time_sequence, self.elapsed_time):
                f.write(''+str(t)+' '+str(seq[0])+' '+str(seq[1])+'\n')

        GT_util.save_figure([self.elapsed_time], [zip(*self.time_sequence)[0]], [zip(*self.time_sequence)[1]], [str(self.wl)], 'Temporal Stability', 'Elapsed Time (seconds)', 'Signal (V)', out_file+'.png')


    def run(self, files):
        mesg = 'Bully! I should be delighted to undertake this endeavor for you!  Haunting optics benches gets so dreadfully boring'
        self.mailer.send_message(mesg)
        status = self.twitterer.tweet(mesg)
        print status
        
        data_file = files[0]
        self.get_Background()

        self.get_TimeSequence(data_file)
        mesg = 'Zounds!  I have finished the stability test!'
        self.mailer.send_message(mesg)
        status = self.twitterer.tweet(mesg)
        print status

class TransmissionEfficiency( Mode ):
    ''' Tests the straight-through transmission efficiency of a slab of non-dispersive material. '''
    def __init__(self, inst_suite, parameters):
        self.inst_suite = inst_suite
        self.wl = []
        self.parameters = parameters
        self.mailer = None
        self.twitterer = GT_util.Twitterer()

    def setup(self):
        wl_start = None
        while not wl_start:
            try:
                wl_start = float(raw_input('Enter the wavelength starting point (in nm): '))
            except:
                print 'Error!  Enter a valid wavelength!'

        wl_stop = None
        while not wl_stop:
            try:
                wl_stop = float(raw_input('Enter the wavelength stopping point (in nm): '))
            except:
                print 'Error!  Enter a valid wavelength!'

        nwl = None
        while not nwl:
            try:
                nwl = float(raw_input('Enter the number of intermediate wavelength points you wish to test : '))
                if nwl < 0:
                    nwl = None
                    print 'Error!  You must enter a positive number!'
            except:
                print 'Error!  Please enter a sensible number!'
        
        self.wl = numpy.linspace(wl_start, wl_stop, nwl)

	    #Sets up the email
        while (not self.mailer):
            try:
                emails = raw_input('Enter the email(s) you wish to notify when user intervention is required (separate with commas) :')
                self.mailer = GT_util.Mailer(self.parameters["MAILSERVER"], self.parameters["SENDER"], emails.split(','))
            except:
                print 'Error! Enter a valid email address!'

        # Gets the user's desired Signal-to-Noise ratio
        a = True
        while a == True:
            try:
                self.SNR = float(raw_input('Enter your desired Signal-to-Noise :'))
                a = False
            except:
                print "Error! Please enter a sensible number!"
                
        self.filter = -1
        while ( (self.filter <= 0) or (self.filter >= 7) ):
            try:
                self.filter = int(raw_input('Which Filter do you want to use?'))
            except:
                self.filter = -1
        self.inst_suite.mono.set_filter(self.filter)

        self.slit_width = -1
        while ( (self.slit_width < 50) or (self.slit_width > 2000) ):
            try:
                self.slit_width = int(raw_input("Enter the monochromator slit width :"))
            except:
                self.slit_width = -1

        self.inst_suite.mono.set_slit_width(self.slit_width)

    def get_Background(self):
        self.background = None

        # Set monochromator to correct wavelength
        self.inst_suite.mono.goto_wavelength(self.wl[0])

        # Open filter wheel
        self.inst_suite.mono.open_shutter()

        y = []

        test_angles = numpy.linspace(-0.5, 0.5, 11)
        # Find middle of PSF
        print "Finding Middle of beam"
        for angle in test_angles:
            txt_out = "Angle = "+str(angle)
            sys.stdout.write(txt_out)
            self.inst_suite.motor.goto(angle)
            time.sleep(self.parameters['SETTLE_TIME'])
            yval = self.inst_suite.srs.measure_const_SNR(10)
            y.append([yval[0], yval[1]])
            for i in numpy.arange(len(txt_out)):
                sys.stdout.write('\b')
            sys.stdout.flush()

        # Find the maximum value of y, go to that beta angle
        order = numpy.array(zip(*y)[0]).argsort()
        self.inst_suite.motor.goto(test_angles[order[-1]])

        # Close filter wheel
        self.inst_suite.mono.close_shutter()

        # Take background
        self.background = self.inst_suite.srs.measure_const_SNR(self.SNR)

    def get_Reference(self):
        self.reference = []
        
        # Open filter wheel
        self.inst_suite.mono.open_shutter()

        # Set monochromator to next wavelength, repeat (with step below)
        for wl in self.wl:
            self.inst_suite.mono.goto_wavelength(wl)
            time.sleep(self.parameters['SETTLE_TIME'])
            # Take reference measurement
            self.reference.append(self.inst_suite.srs.measure_const_SNR(self.SNR))
            

    def get_Transmission(self, out_file):
        self.measurement = []
        self.transmission = []

        # Open filter wheel
        self.inst_suite.mono.open_shutter()

        # Set monochromator to next wavelength, repeat (with step below)
        for wl in self.wl:
            self.inst_suite.mono.goto_wavelength(self.wl[0])
            time.sleep(self.parameters['SETTLE_TIME'])
            # Take measurement
            self.measurement.append(self.inst_suite.srs.measure_const_SNR(self.SNR))

        # Compute transmission
        for ref, meas in zip(self.reference, self.measurement):
            self.transmission.append([(meas[0]-self.background[0])/(ref[0]-self.background[0]),( (meas[1]/(ref[0]-self.background[0]))**2.0+(ref[1]*((meas[0]-self.background[0])/(ref[0]-self.background[0])**2))**2.0+(self.background[1]*(meas[0]-ref[0])/(ref[0]-self.background[0])**2.0)**2.0)**(0.5)])

        GT_util.save_figure([self.wl, self.wl], [zip(*self.reference)[0], zip(*self.measurement)[0]], [zip(*self.reference)[1], zip(*self.measurement)[1]], ['Reference Beam', self.material_name], 'Raw Transmission', 'Wavelength (nanometers)', 'Signal (V)', out_file+'_raw_.png')

        GT_util.save_figure([self.wl], [zip(*self.transmission)[0]], [zip(*self.transmission)[1]], [self.material_name+' Transmission'], 'Normalized Transmission', 'Wavelength (nanometers)', 'Transmission', out_file+'_.png')

        with open(out_file, 'a') as f:
            f.write('Background = '+str(self.background[0])+' +/- '+str(self.background[1])+'\n')
            for wl, ref, meas, trans in zip(self.wl, self.reference, self.measurement, self.transmission):
                f.write(str(wl)+', '+str(ref[0])+' +/- '+str(ref[1])+', '+str(meas[0])+' +/- '+str(meas[1])+', '+str(trans[0])+' +/- '+str(trans[1])+'\n')

    def run(self, files):
        ''' performs the transmission efficiency scan '''

        transmission_file = files[0]
        
        #background scan
        self.get_Background()

        #Reference scan
        self.get_Reference()

        #Transmission scan
        self.get_Transmission(transmission_file)

class BlazeEfficiency( Mode ):
    def __init__( self, inst_suite, parameters ):
        """ Sets up the blaze efficiency mode"""
        self.inst_suite = inst_suite
        self.wl = []
        self.beta = []
        self.alpha = None
        self.order = None
        self.disp_el = None
        self.parameters = parameters
        self.mailer = None
        self.twitterer = GT_util.Twitterer()

    def setup(self, dispersive_element):
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

        if ( dispersive_element == 'GRATING'):
            a = True
            while a == True:
                try:
                    blaze = float(raw_input('Enter the Blaze angle of the Grating : '))
                    gconst = float(raw_input('Enter the Grating constant (g/mm) : '))
                    n = float(raw_input('Enter the index of refraction : '))
                    name = raw_input('Enter the name/SN of the grating: ')
                    self.disp_el = Grating(n, blaze, name=name, d=gconst)
                    a = False
                except:
                    print "Error! Please enter a sensible number!"
        elif ( dispersive_element == 'GRISM'):
            delta = None
            blaze = None
            while( not delta ):
                try:
                    delta = float(raw_input('Enter wedge angle of the prism :'))
                except:
                    delta = None
                    
            while( not blaze ):
                try:
                    blaze = float(raw_input('Enter groove angle (blaze) :'))
                except:
                    blaze = None
                    
            a = True
            while a == True:
                try:
                    gperiod = float(raw_input('Enter the Groove Period (microns/groove) : '))
                    n = float(raw_input('Enter the index of refraction : '))
                    name = raw_input('Enter the name/SN of the grating: ')
                    self.disp_el = Grism(n, delta, blaze, sigma=gperiod, name=name)
                    a = False
                except:
                    print "Error! Please enter a sensible number!"


        while( not self.alpha ):
            try:
                self.alpha = float(raw_input('Enter Alpha angle (degrees) :'))
            except:
                self.alpha = None

        while (not self.order):
            try:
                self.order = float(raw_input('Enter the order you wish to test: '))
            except:
                self.order = None

        # Calculates the beta angles for all the wavelengths to be tested
        self.beta = self.disp_el.calc_beta(self.wl/1000.0, self.alpha, self.order)

	    #Sets up the email
        while (not self.mailer):
            try:
                emails = raw_input('Enter the email(s) you wish to notify when user intervention is required (separate with commas) :')
                self.mailer = GT_util.Mailer(self.parameters["MAILSERVER"], self.parameters["SENDER"], emails.split(','))
            except:
                print 'Error! Enter a valid email address!'

        # Gets the user's desired Signal-to-Noise ratio
        a = True
        while a == True:
            try:
                self.SNR = float(raw_input('Enter your desired Signal-to-Noise :'))
                a = False
            except:
                print "Error! Please enter a sensible number!"
                
        self.filter = -1
        while ( (self.filter <= 0) or (self.filter >= 7) ):
            try:
                self.filter = int(raw_input('Which Filter do you want to use?'))
            except:
                self.filter = -1
        self.inst_suite.mono.set_filter(self.filter)

        self.slit_width = -1
        while ( (self.slit_width < 50) or (self.slit_width > 2000) ):
            try:
                self.slit_width = int(raw_input("Enter the monochromator slit width :"))
            except:
                self.slit_width = -1

        self.inst_suite.mono.set_slit_width(self.slit_width)
        
        # Sets the self.background variable to an empty list
        self.background = []

    def get_Background(self, bg_file):
        x = numpy.linspace(self.parameters['REF_START_BETA'], self.parameters['REF_STOP_BETA'], self.parameters['N_REF_PTS'])
        y = []
        self.inst_suite.mono.close_shutter()
        junk = self.inst_suite.srs.measure_const_SNR(self.SNR)
        for b in x:
            txt_out = "Angle = "+str(b)
            sys.stdout.write(txt_out)
            self.inst_suite.motor.goto(b)
            time.sleep(self.parameters['SETTLE_TIME'])
            y.append(self.inst_suite.srs.measure_const_SNR(self.SNR))
            for i in numpy.arange(len(txt_out)):
                sys.stdout.write('\b')
            for i in numpy.arange(len(txt_out)):
                sys.stdout.write(' ')
            for i in numpy.arange(len(txt_out)):
                sys.stdout.write('\b')
            sys.stdout.flush()

        self.background = [x, y]

        GT_util.save_figure([x], [zip(*y)[0]], [zip(*y)[1]], ['Background Signal'], 'Background', 'Angle (degrees)', 'Signal (V)', bg_file+'.png')
        with open(bg_file, 'a') as f:
            for angle, reading, dreading in zip(x, zip(*y)[0], zip(*y)[1]):
                f.write(str(angle)+', '+str(reading)+', '+str(dreading)+'\n')
        self.inst_suite.mono.open_shutter()
        

    def get_PSF(self, psf_file):
        x = numpy.linspace(self.parameters['REF_START_BETA'], self.parameters['REF_STOP_BETA'], self.parameters['N_REF_PTS'])
        xpts = []
        ypts = []
        dypts = []
        names = []
        retval = []

        print self.wl[0]

        self.inst_suite.mono.goto_wavelength(self.wl[0])
        y = []

        print "Starting PSF Scan"
        for b, y_bkgnd in zip(self.background[0], self.background[1]):
            txt_out = "Angle = "+str(b)
            sys.stdout.write(txt_out)
            self.inst_suite.motor.goto(b)
            time.sleep(self.parameters['SETTLE_TIME'])
            yval = self.inst_suite.srs.measure_const_SNR(self.SNR)
            y.append([yval[0] - y_bkgnd[0], (yval[1]**2 + y_bkgnd[1]**2)**(0.5)])
            for i in numpy.arange(len(txt_out)):
                sys.stdout.write('\b')
            for i in numpy.arange(len(txt_out)):
                sys.stdout.write(' ')
            for i in numpy.arange(len(txt_out)):
                sys.stdout.write('\b')
            sys.stdout.flush()

        retval.append([x, y])

        # Find the maximum value of y, go to that beta angle
        #print y
        #print self.background[0]
        order = numpy.array(zip(*y)[0]).argsort()
        self.inst_suite.motor.goto(self.background[0][order[-1]])
        bkgnd = self.background[1][order[-1]]

        print "Sleeping for 10 seconds!"
        time.sleep(10)

        junk = self.inst_suite.srs.measure_const_SNR(self.SNR)
        print junk
        
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
        GT_util.save_figure([retval[0][0]], [zip(*retval[0][1])[0]], [zip(*retval[0][1])[1]], [str(self.wl[0])], 'Background Subtracted Reference PSF', 'Angle (degrees)', 'Signal (V)', psf_file+'.png')
        return retval

    def sweep(self, beta, psf, scale_factor, wl, raw_file):
        """ sweeps the detector through a range of angles centered on the predicted angle.
            Fits the observed data to the reference PSF """

        # creates array of beta angles to check
        x = numpy.linspace(beta+self.parameters['SWEEP_START_BETA'], beta+self.parameters['SWEEP_STOP_BETA'], self.parameters['N_SWEEP_PTS'])
        x_zoom = numpy.linspace(beta-0.15, beta+0.15, 11)
        below = scipy.where(x < beta-0.15)
        above = scipy.where(x > beta+0.15)
        new_x = []
        for bm in below[0]:
            new_x.append(x[bm])
        for i in x_zoom:
            new_x.append(i)
        for bm in above[0]:
            new_x.append(x[bm])

        x = numpy.array(new_x)
        y = []
        for b in x:
            txt_out = "Angle = "+str(b)
            sys.stdout.write(txt_out)
            self.inst_suite.motor.goto(b)
            time.sleep(self.parameters['SETTLE_TIME'])
            y.append(self.inst_suite.srs.measure_const_SNR(self.SNR))
            for i in numpy.arange(len(txt_out)):
                sys.stdout.write('\b')
            for i in numpy.arange(len(txt_out)):
                sys.stdout.write(' ')
            for i in numpy.arange(len(txt_out)):
                sys.stdout.write('\b')
            sys.stdout.flush()

        sys.stdout.write('\n')
        #print "PSF :", psf
        #print "PSF[1] :", psf[1]
        #print "Scale Factor :", scale_factor
        # Scales the PSF by the scale factor
        y_psf = numpy.array([s[0] for s in psf[1]])*scale_factor
        y_obs = [s[0] for s in y]

        #integrate under the beam
        beam = scipy.integrate.simps(y_psf, psf[0])
        print 'Beam signal = ', beam

        #integrate under the order
        observations = scipy.integrate.simps(y_obs, x)
        print 'Observed signal = ', observations
        
        retval = (1.0*observations)/(beam*1.0)

        print 'Efficiency = ', retval

        '''
        # beam and observations are interpolation objects.
        #print psf[0]
        #print y_psf
        beam = scipy.interpolate.interpolate.interp1d(psf[0], y_psf, bounds_error=False)
        observations = scipy.interpolate.interpolate.interp1d(x, y_obs, bounds_error=False)

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
            def errfunc(p, beta, y):
                newpsf = fitfunc(p, beta)
                #print "P = ", p, ", beta = ", beta
                bm = scipy.where( (numpy.isfinite(newpsf)) & (numpy.isfinite(y)) )
                #print bm
                return newpsf[bm] - y[bm]
            p1, success = scipy.optimize.leastsq(errfunc, pguess, args = (beta, y_new))
            #p1 = [0.3, 0.1]
            print "Fitting Coefficients are : ", p1
        except:
            print "ERROR!  Fit did not converge!"
            p1 = [0.0, 0.0]
        '''

        with open(raw_file, 'a') as f:
            f.write(str(wl)+'\n')
            for angle, reading, dreading in zip(x, zip(*y)[0], zip(*y)[1]):
                f.write(str(angle)+', '+str(reading)+', '+str(dreading)+'\n')
                    
        return retval, x, zip(*y)[0], zip(*y)[1]

    def run(self, files):
        mesg = 'Bully! I should be delighted to undertake this endeavor for you!  Haunting optics benches gets so dreadfully boring'
        self.mailer.send_message(mesg)
        status = self.twitterer.tweet(mesg)
        print status

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
            f.write('Dispersive Element Properties\n========================\n')
            if (type(self.disp_el) == Grism):
                f.write('Type: Grism\n')
            if (type(self.disp_el) == Grating):
                f.write('Type: Grating\n')
            f.write('Serial#/Identifier : '+str(self.disp_el.name)+'\n')
            f.write('Grating Constant : ' + str(self.disp_el.d)+' grooves/mm\n')
            f.write('                 = ' + str(self.disp_el.sigma)+' microns/groove\n')
            if (type(self.disp_el) == Grism):
                f.write('Wedge Angle : ' +str(self.disp_el.delta)+' degrees\n')
                f.write('index of refraction : '+str(self.disp_el.n)+'\n')
            f.write('Blaze Angle : ' + str(self.disp_el.blaze)+' degrees\n')
            f.write('Alpha Angle : ' + str(self.alpha)+' degrees\n')
            f.write('Diffraction Order : ' + str(self.order)+'\n')

        # Zeros out the readings variable
        self.readings = []
        self.d_readings = []

        # Ensure that all necessary variables are defined
        if (len(self.beta) <= 0) or (len(self.wl) <= 0):
            print "Error!  Set up the wavelength and alpha points!"
        else:

            # Measures the background
            self.get_Background(bg_file)
            mesg = 'Capital! The background scan has concluded!'
            status = self.twitterer.tweet(mesg)
            print status

            # Measures the reference PSF
            psfs = self.get_PSF(psf_file)

            # Gets ready to take grating data
            mesg = 'Non-corporeality has its drawbacks... I require assistance!  Would you be so kind as to place your grating in the beam!'
            self.mailer.send_message(mesg)
            status = self.twitterer.tweet(mesg)
            print status

            ch = raw_input('Place Grating in beam!')
            if len(self.wl) < 4:
                n_milestones = len(self.wl)
            else:
                n_milestones = 4
            milestones = numpy.linspace(self.wl[0], self.wl[-1], n_milestones+1)
            last_milestone = 0
            for pair in zip(self.wl, self.beta, psfs[1:]):
                wl = int(pair[0])
                beta = float(pair[1])
                psf_scale_factor = pair[2]

                if (len(milestones) > last_milestone+1):
                    if (wl > milestones[last_milestone+1]):
                        last_milestone += 1
                        fraction = str(((1.0*last_milestone)/(1.0*n_milestones)*100.0))
                        mesg = 'Most momentus! I am '+fraction+'% done with the task! Huzzah!'
                        status = self.twitterer.tweet(mesg)

                print "Monochromator, please go to ", wl, " nanometers"
                self.inst_suite.mono.goto_wavelength(wl)
                time.sleep(10.0)
                while (abs(self.inst_suite.mono.read_wavelength() - wl) > 1.0):
                    #print self.inst_suite.mono.read_wavelength()
                    time.sleep(2.0)
                print "Motor Controller, please go to ", beta+self.alpha, " degrees"
                coeff, x, y, dy = self.sweep(beta+self.alpha, psfs[0], (psf_scale_factor[1]/psfs[1][1]), wl, raw_file)
                self.readings.append(coeff)
                self.d_readings.append(coeff/10.0)
                xpts.append(x)
                ypts.append(y)
                dypts.append(dy)
                names.append(str(wl))
                with open(data_file, 'a') as f:
                    f.write(str(wl)+', '+str(beta)+', '+str(coeff)+'\n')
                print self.readings[-1]
            GT_util.save_figure(xpts, ypts, dypts, names, 'Raw Data', 'Angle (degrees)', 'Signal (V)', raw_file+'.png')
            print 'wl: ', self.wl
            print 'Readings : ', self.readings
            print 'unc : ', self.d_readings
            print 'name : ', self.disp_el.name
            GT_util.save_figure([self.wl], [self.readings], [self.d_readings], [self.disp_el.name], 'Blaze Efficiency', 'Wavelength (nm)', 'Efficiency', data_file+'.png')
        print "Done!"

        mesg = 'Great Scott! I have finshed testing your grating.  Now off to do some haunting!'
        self.mailer.send_message(mesg)
        status = self.twitterer.tweet(mesg)
        print status
        self.mailer.close()

class PSF_tester( Mode ):
    def __init__( self, inst_suite, parameters ):
        """ Sets up the blaze efficiency mode"""
        self.inst_suite = inst_suite
        self.wl = []
        self.parameters = parameters
        self.mailer = None
        self.twitterer = GT_util.Twitterer()

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

        #Sets up the email
        while (not self.mailer):
            try:
                emails = raw_input('Enter the email(s) you wish to notify when user intervention is required (separate with commas) :')
                self.mailer = GT_util.Mailer(self.parameters["MAILSERVER"], self.parameters["SENDER"], emails.split(','))
            except:
                print 'Error! Enter a valid email address!'

        # Gets the user's desired Signal-to-Noise ratio
        a = True
        while a == True:
            try:
                self.SNR = float(raw_input('Enter your desired Signal-to-Noise :'))
                a = False
            except:
                print "Error! Please enter a sensible number!"
                
        self.filter = -1
        while ( (self.filter <= 0) or (self.filter >= 7) ):
            try:
                self.filter = int(raw_input('Which Filter do you want to use?'))
            except:
                self.filter = -1
        self.inst_suite.mono.set_filter(self.filter)

        self.slit_width = -1
        while ( (self.slit_width < 50) or (self.slit_width > 2000) ):
            try:
                self.slit_width = int(raw_input("Enter the monochromator slit width :"))
            except:
                self.slit_width = -1

        self.inst_suite.mono.set_slit_width(self.slit_width)
        
        # Sets the self.background variable to an empty list
        self.background = []

    def get_Background(self, bg_file):
        x = numpy.linspace(self.parameters['REF_START_BETA'], self.parameters['REF_STOP_BETA'], self.parameters['N_REF_PTS'])
        y = []
        self.inst_suite.mono.close_shutter()
        for b in x:
            txt_out = "Angle = "+str(b)
            sys.stdout.write(txt_out)
            self.inst_suite.motor.goto(b)
            time.sleep(self.parameters['SETTLE_TIME'])
            y.append(self.inst_suite.srs.measure_const_SNR(self.SNR))
            for i in numpy.arange(len(txt_out)):
                sys.stdout.write('\b')
            for i in numpy.arange(len(txt_out)):
                sys.stdout.write(' ')
            for i in numpy.arange(len(txt_out)):
                sys.stdout.write('\b')
            sys.stdout.flush()

        self.background = [x, y]

        GT_util.save_figure([x], [zip(*y)[0]], [zip(*y)[1]], ['Background Signal'], 'Background', 'Angle (degrees)', 'Signal (V)', bg_file+'.png')
        with open(bg_file, 'a') as f:
            for angle, reading, dreading in zip(x, zip(*y)[0], zip(*y)[1]):
                f.write(str(angle)+', '+str(reading)+', '+str(dreading)+'\n')
        self.inst_suite.mono.open_shutter()
        

    def get_PSF(self, psf_file):
        x = numpy.linspace(self.parameters['REF_START_BETA'], self.parameters['REF_STOP_BETA'], self.parameters['N_REF_PTS'])
        xpts = []
        ypts = []
        dypts = []
        names = []
        retval = []

        print self.wl[0]

        self.inst_suite.mono.goto_wavelength(self.wl[0])
        y = []

        print "Starting PSF Scan"
        for b, y_bkgnd in zip(self.background[0], self.background[1]):
            txt_out = "Angle = "+str(b)
            sys.stdout.write(txt_out)
            self.inst_suite.motor.goto(b)
            time.sleep(self.parameters['SETTLE_TIME'])
            yval = self.inst_suite.srs.measure_const_SNR(self.SNR)
            y.append([yval[0] - y_bkgnd[0], (yval[1]**2 + y_bkgnd[1]**2)**(0.5)])
            for i in numpy.arange(len(txt_out)):
                sys.stdout.write('\b')
            for i in numpy.arange(len(txt_out)):
                sys.stdout.write(' ')
            for i in numpy.arange(len(txt_out)):
                sys.stdout.write('\b')
            sys.stdout.flush()

        retval.append([x, y])

        # Find the maximum value of y, go to that beta angle
        #print y
        #print self.background[0]
        order = numpy.array(zip(*y)[0]).argsort()
        self.inst_suite.motor.goto(self.background[0][order[-1]])
        bkgnd = self.background[1][order[-1]]

        print "Sleeping for 10 seconds!"
        time.sleep(10)

        junk = self.inst_suite.srs.measure_const_SNR(self.SNR)
        print junk
        
        for wl in self.wl:
            print "wavelength = "+str(wl)+"\n"
            self.inst_suite.mono.goto_wavelength(wl)
            time.sleep(self.parameters['SETTLE_TIME']*2.0)
            reading = self.inst_suite.srs.measure_const_SNR(self.SNR)
            print reading
            retval.append([wl, reading[0]-bkgnd[0], (reading[1]**2.0+bkgnd[1]**2.0)**(0.5)])

        with open(psf_file, 'a') as f:
            for angle, reading in zip(*retval[0]):
                f.write(str(angle)+', '+str(reading[0])+', '+str(reading[1])+'\n')
            for wavelength, reading, dreading in retval[1:]:
                f.write(str(wavelength)+ ', '+str(reading)+', '+str(dreading)+'\n')
        GT_util.save_figure([retval[0][0]], [zip(*retval[0][1])[0]], [zip(*retval[0][1])[1]], [str(self.wl[0])], 'Background Subtracted Reference PSF', 'Angle (degrees)', 'Signal (V)', psf_file+'.png')
        return retval

    def sweep(self, beta, psf, scale_factor, wl, raw_file):
        """ sweeps the detector through a range of angles centered on the predicted angle.
            Fits the observed data to the reference PSF """

        # creates array of beta angles to check
        x = numpy.linspace(beta+self.parameters['SWEEP_START_BETA'], beta+self.parameters['SWEEP_STOP_BETA'], self.parameters['N_SWEEP_PTS'])
        y = []
        for b in x:
            txt_out = "Angle = "+str(b)
            sys.stdout.write(txt_out)
            self.inst_suite.motor.goto(b)
            time.sleep(self.parameters['SETTLE_TIME'])
            y.append(self.inst_suite.srs.measure_const_SNR(self.SNR))
            for i in numpy.arange(len(txt_out)):
                sys.stdout.write('\b')
            for i in numpy.arange(len(txt_out)):
                sys.stdout.write(' ')
            for i in numpy.arange(len(txt_out)):
                sys.stdout.write('\b')
            sys.stdout.flush()

        sys.stdout.write('\n')
        #print "PSF :", psf
        #print "PSF[1] :", psf[1]
        #print "Scale Factor :", scale_factor
        # Scales the PSF by the scale factor
        y_psf = numpy.array([s[0] for s in psf[1]])*scale_factor
        y_obs = [s[0] for s in y]

        #integrate under the beam
        beam = scipy.integrate.simps(y_psf, psf[0])
        print 'Beam signal = ', beam

        #integrate under the order
        observations = scipy.integrate.simps(y_obs, x)
        print 'Observed signal = ', observations
        
        retval = (1.0*observations)/(beam*1.0)

        print 'Efficiency = ', retval

        '''
        # beam and observations are interpolation objects.
        #print psf[0]
        #print y_psf
        beam = scipy.interpolate.interpolate.interp1d(psf[0], y_psf, bounds_error=False)
        observations = scipy.interpolate.interpolate.interp1d(x, y_obs, bounds_error=False)

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
            def errfunc(p, beta, y):
                newpsf = fitfunc(p, beta)
                #print "P = ", p, ", beta = ", beta
                bm = scipy.where( (numpy.isfinite(newpsf)) & (numpy.isfinite(y)) )
                #print bm
                return newpsf[bm] - y[bm]
            p1, success = scipy.optimize.leastsq(errfunc, pguess, args = (beta, y_new))
            #p1 = [0.3, 0.1]
            print "Fitting Coefficients are : ", p1
        except:
            print "ERROR!  Fit did not converge!"
            p1 = [0.0, 0.0]
        '''

        with open(raw_file, 'a') as f:
            f.write(str(wl)+'\n')
            for angle, reading, dreading in zip(x, zip(*y)[0], zip(*y)[1]):
                f.write(str(angle)+', '+str(reading)+', '+str(dreading)+'\n')
                    
        return retval, x, zip(*y)[0], zip(*y)[1]

    def run(self, files):
        mesg = 'Bully! I should be delighted to undertake this endeavor for you!  Haunting optics benches gets so dreadfully boring'
        self.mailer.send_message(mesg)
        status = self.twitterer.tweet(mesg)
        print status

        # Sets up the data files
        data_file = files[0]
        raw_file = files[1]
        bg_file = files[2]
        psf_file = files[3]
        xpts = []
        ypts = []
        dypts = []
        names = []
        
        # Zeros out the readings variable
        self.readings = []
        self.d_readings = []
        self.beta = [0]*len(self.wl)

        # Ensure that all necessary variables are defined
        if (len(self.beta) <= 0) or (len(self.wl) <= 0):
            print "Error!  Set up the wavelength and alpha points!"
        else:

            # Measures the background
            self.get_Background(bg_file)
            mesg = 'Capital! The background scan has concluded!'
            status = self.twitterer.tweet(mesg)
            print status
            
            # Measures the reference PSF
            psfs = self.get_PSF(psf_file)

            # Gets ready to take grating data
            if len(self.wl) < 4:
                n_milestones = len(self.wl)
            else:
                n_milestones = 4
            milestones = numpy.linspace(self.wl[0], self.wl[-1], n_milestones+1)
            last_milestone = 0
            for pair in zip(self.wl, self.beta, psfs[1:]):
                wl = int(pair[0])
                beta = float(pair[1])
                psf_scale_factor = pair[2]

                if (len(milestones) > last_milestone+1):
                    if (wl > milestones[last_milestone+1]):
                        last_milestone += 1
                        fraction = str(((1.0*last_milestone)/(1.0*len(milestones))*100.0))
                        mesg = 'Most momentus! I am '+fraction+'% done with the task! Huzzah!'
                        status = self.twitterer.tweet(mesg)

                print "Monochromator, please go to ", wl, " nanometers"
                self.inst_suite.mono.goto_wavelength(wl)
                time.sleep(10.0)
                while (abs(self.inst_suite.mono.read_wavelength() - wl) > 1.0):
                    #print self.inst_suite.mono.read_wavelength()
                    time.sleep(2.0)
                print "Motor Controller, please go to ", beta, " degrees"
                print "psf scale factor : ", psf_scale_factor[1], ' ', psfs[1][1], ' ', psf_scale_factor[1]/psfs[1][1]
                coeff, x, y, dy = self.sweep(beta, psfs[0], (psf_scale_factor[1]/psfs[1][1]), wl, raw_file)
                self.readings.append(coeff)
                self.d_readings.append(coeff/10.0)
                xpts.append(x)
                ypts.append(y)
                dypts.append(dy)
                names.append(str(wl))
                with open(data_file, 'a') as f:
                    f.write(str(wl)+', '+str(beta)+', '+str(coeff)+'\n')
                print self.readings[-1]
            GT_util.save_figure(xpts, ypts, dypts, names, 'Raw Data', 'Angle (degrees)', 'Signal (V)', raw_file+'.png')
            print 'wl: ', self.wl
            print 'Readings : ', self.readings
            print 'unc : ', self.d_readings
            GT_util.save_figure([self.wl], [self.readings], [self.d_readings], ['PSF'], 'PSF Stability', 'Wavelength (nm)', 'Efficiency', data_file+'.png')
        print "Done!"

        mesg = 'Great Scott! I have finshed testing your grating.  Now off to do some haunting!'
        self.mailer.send_message(mesg)
        status = self.twitterer.tweet(mesg)
        print status
        self.mailer.close()


class ScatteredLight( Mode ):
    """ """
    def __init__( self, inst_suite ):
        """Sets up the scatter-scan mode"""
        self.inst_suite = inst_suite
        self.wl = []
        self.beta = []

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
