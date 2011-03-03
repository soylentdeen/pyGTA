import serial
import time
import numpy

class SRS510(object):
    '''
function: __init__(self, sim, port)
==========================================
Purpose: To initialize the serial port and the SRS510
Arguments:
  sim - (boolean) If True, we are simulating the data
  port - (integer) address of serial port

Returns:
  Nothing
    '''
    def __init__(self, sim, parameters):
        self.ser = serial.Serial(parameters['LOCKIN_COMM'])
        self.ser.baudrate = parameters['LOCKIN_BAUDRATE']
        self.ser.parity = parameters['LOCKIN_PARITY']
        self.ser.stopbits = parameters['LOCKIN_STOPBITS']
        self.ser.timeout = parameters['LOCKIN_TIMEOUT']

        self.t_min = parameters['LOCKIN_T_MIN']
        self.t_max = parameters['LOCKIN_T_MAX']

        self.sensitivity = {1:10e-9, 2:20e-9, 3:50e-9, 4:100e-9, 5:200e-9, 6:500e-9, 7:1e-6, 8:2e-6, 9:5e-6, 10:10e-6, 11:20e-6, 12:50e-6, 13:100e-6, 14:200e-6, 15:500e-6, 16:1e-3, 17:2e-3, 18:5e-3, 19:10e-3, 20:20e-3, 21:50e-3, 22:100e-3, 23:200e-3, 24:500e-3}
        self.gain = int(self.get_gain_setting())
        #print 'Gain setting is :' + str(self.gain)
        self.autorange = (parameters['LOCKIN_AUTORANGE']=='TRUE')
        self.previous_reading = 0.5               # stores the previous measurement.  To be used to determine if we need to switch ranges


    '''
measure_const_SNR(self, SNR_min)
Purpose: Starts a measurement, and dynamically adjusts the sensitivity, dynamic range, and time constants of the
lock-in amplifier to achieve optimum signal.
Arguments:
  SNR_min - (float) Minimum acceptable Signal-to-Noise ratio

Returns:
  measurement - (float) mean value of all the measurements
  noise - (float) noise associated with the measurement
    '''
    '''
    def measure_const_SNR(self, SNR_min):
     measurements = []
     snr = 1.0
     tstart = time.time()
     t = time.time()-tstart

     temp = self.ser.read(0)
     while ( ((snr <= SNR_min) or (t <= self.t_min)) and (t <= self.t_max)):
         self.ser.write('Q\r\n')
         time.sleep(0.35)
         nchar = self.ser.inWaiting()
         temp = self.ser.read(nchar)
         if ( temp[len(temp)-1] == '\r'):
             measurements.append(float(temp.split()[0]))
         else:
             time.sleep(0.25)
             nchar = self.ser.inWaiting()
             junk = self.ser.read(nchar)
             measurements.append(float(temp+junk))
         signal = numpy.mean(measurements)
         noise = numpy.std(measurements)/(len(measurements))**0.5
         if (len(measurements) >= 2):
             snr = numpy.abs(signal/noise)
         else:
             snr = 1.0
         t = time.time() - tstart

     avg = numpy.mean(measurements)
     std = numpy.std(measurements)

     return [avg, std/(len(measurements))**0.5]
    '''
    def get_gain_setting(self):
        self.ser.write('G\r\n')
        time.sleep(0.35)
        nchar = self.ser.inWaiting()
        temp = self.ser.read(nchar)
        if( temp[len(temp)-1] == '\r'):
            self.gain = float(temp)
        else:
            time.sleep(0.25)
            nchar = self.ser.inWaiting()
            junk = self.ser.read(nchar)
            self.gain = float(temp+junk)

        return self.gain

    '''
function: measure_const_SNR(self, SNR_min)
=================================
Purpose: Starts a measurement
Arguments:
  SNR_min - (float) Minimum acceptable Signal-to-Noise ratio

Returns:
  measurement - (float) mean value of all the measurements
  noise - (float) noise associated with the measurement
    '''
    def measure_const_SNR(self, SNR_min):
        measurements = []
        snr = 1.0
        tstart = time.time()
        t = time.time()-tstart

        temp = self.ser.read(0)      # clears out anything in the buffer

        # Need to see if previous measurement warrants increasing the gain
        full_scale = self.sensitivity[self.gain]
        next_full_scale = self.sensitivity[self.gain-1]
        if (self.previous_reading < 0.05*next_full_scale):
            self.gain -= 1
            self.ser.write('G%d\r\nD1\r\n' %int(self.gain))
            time.sleep(0.25)

        while ( ((snr <= SNR_min) or (t <= self.t_min)) and (t <= self.t_max)):
            # See if we are saturated
            if ( self.autorange ):
                t1 = time.time()
                self.ser.write('Y4\r\n')
                time.sleep(0.35)
                nchar = self.ser.inWaiting()
                temp = self.ser.read(nchar)
                if( temp[len(temp)-1] == '\r'):
                    ovrld = float(temp)
                else:
                    time.sleep(0.25)
                    nchar = self.ser.inWaiting()
                    junk = self.ser.read(nchar)
                    ovrld = float(temp+junk)
                if ( ovrld == 1 and self.gain < 24):
                    #print 'Overloaded!'
                    self.gain += 1
                    self.ser.write('G%d\r\n' % int(self.gain))
                    time.sleep(2.0)
                    t2 = time.time()
                    tstart += t2-t1

            # Go ahead with the measurement
            self.ser.write('Q\r\n')
            time.sleep(0.35)
            nchar = self.ser.inWaiting()
            temp = self.ser.read(nchar)
            if ( temp[len(temp)-1] == '\r'):
                measurements.append(float(temp.split()[0]))
            else:
                time.sleep(0.25)
                nchar = self.ser.inWaiting()
                junk = self.ser.read(nchar)
                measurements.append(float(temp+junk))
                
            signal = numpy.mean(measurements)

            # Calculate signal-to-noise ratio
            noise = numpy.std(measurements)/(len(measurements))**0.5
            if (len(measurements) >= 2):
                snr = numpy.abs(signal/noise)
            else:
                snr = 1.0
            t = time.time() - tstart

        avg = numpy.mean(measurements)
        std = numpy.std(measurements)
        self.previous_reading = avg
        return [avg, std/(len(measurements))**0.5]


    def close(self):
        self.ser.close()
