import serial
import numpy
import time

class Monochromator( object ):
    '''
    function: __init__(self, sim)
    =============================
    Purpose: To initialize the serial port and the monochromator object
    Arguments:
        sim - (boolean variable)  If True, we are simulating the data.  If false, we are taking live data.

    Returns:
        Nothing
    ''' 

    def __init__(self, sim, parameters):
        self.sim = sim
        self.filters = {1:"empty_1", 2:"F1", 3:"empty_2", 4:"AB3190", 5:"AB3300", 6:"Blank"}
        self.selected_filter = 1
        if (self.sim == False):
            NEWLINE_CONVERSION_MAP = ('\n', '\r', '\r\n')
            self.ser = serial.Serial(int(parameters['MONO_COMM']))
            self.ser.bytesize = parameters['MONO_BYTESIZE']
            self.ser.parity = parameters['MONO_PARITY']
            self.ser.timeout = parameters['MONO_TIMEOUT']
            self.ser.baudrate = parameters['MONO_BAUDRATE']
            self.ser.newline = NEWLINE_CONVERSION_MAP
        else:
            self.ser = 'junk'
        #self.move_filter_wheel(3)
            

    '''
    function: encode_wavelength(nm)
    ====================================
    Purpose: to encode the requested wavelength in the format recognized by the monochromator
    Arguments:
    nm- (float) - requested wavelength in nanometers

    Returns:
    ascii character array representation of wavelength in monochromator format.
    '''
    def encode_wavelength(self, nm):
        hundredths = nm*100.00
        h = hundredths // 65536.0
        m = (hundredths - 65536.0*h) // 256.0
        l = (hundredths - 65536.0*h) % 256.0
        return chr(int(h))+chr(int(m))+chr(int(l))

    def encode_slit_width(self, sw):
        h = sw // 256.0
        l = sw % 256.0
        return chr(int(h))+chr(int(l))
    
    '''
    function: decode_wavelength(word)
    ====================================
    Purpose: to decode the wavelength string received from the monochromator
    Arguments:
        word- (3 character string) read from monochromator

    Returns:
    nm-  Wavelength in nanometers.
    '''
    def decode_wavelength(self, word):
        return (65536*ord(word[0]) + 256*ord(word[1]) + ord(word[2]))/100.0


    def goto_wavelength(self, wl):
        self.send_command(16)
        time.sleep(0.5)
        junk = self.read_port()
        self.send_str(self.encode_wavelength(wl))
        time.sleep(0.5)
        junk = self.read_port()
        time.sleep(2.0)
  
    '''
    function: send_command(command)
    ===================================
    Purpose: writes a command to the monochromator
    Arguments:
        command- (int) - base 10 number to be converted to ascii character and sent to

    Returns:
        nothing
    '''
 
    def send_command(self, command):
        if (self.sim == False):
            self.ser.write(chr(command))
        else:
            print command

    '''
    function: send_str(command_string)
    ===================================
    Purpose: writes an ascii string to the monochromator.
    Arguments:
        command_string- (ascii sttring) - ascii string to be sent to the monochromator.  Must already be translated to ascii

    Returns:
        nothing
    '''

    def send_str(self, command_string):
        if (self.sim == False):
            self.ser.write(command_string)
        else:
            print "wrote "+command_string+" to the port (simulated)"

    '''
    function: read_port()
    ==================================
    Purpose: reads characters waiting at the port from the monochromator
    Arguments:
        nothing
    Returns:
        ascii string
    '''
    def read_port(self):
        if (self.sim == False):
            nchr = self.ser.inWaiting()
            counter = 0.0
            while ( (nchr == 0) and (counter < self.ser.timeout) ):
                time.sleep(1.0)
                #print 'waiting...'
                nchr = self.ser.inWaiting()
                counter += 1.0
            if counter < self.ser.timeout:
                return self.ser.read(nchr)
            else:
                return "junk"
        else:
            return "junk"

    def move_filter_wheel(self, position):
        """ Moves filter wheel to desired position """
        print "Moving to Filter Wheel detent : ", position
        self.send_command(15)
        time.sleep(0.5)
        echo = self.read_port()
        self.send_command(position)
        time.sleep(0.5)
        status = self.read_port()
        time.sleep(0.5)
        cancel_byte = self.read_port()
        self.send_command(24)     # Although not mentioned in the manual,
                                  # this <24>D cancel byte is necessary to
                                  # return the DK480 to a non-bricked state.
        time.sleep(2.0)
        self.send_command(24)
        #return echo, status, cancel_byte

    def set_filter(self, filter_position):
        self.selected_filter = filter_position

    def set_slit_width(self, slit_width):
        print "Setting Slit Width to : ", slit_width, " microns"
        self.send_command(14)
        time.sleep(0.5)
        echo = self.read_port()
        self.send_str(self.encode_slit_width(slit_width))
        status = self.read_port()
        time.sleep(3.0)
        cancel_byte = self.read_port()
        self.send_command(24)
        time.sleep(2.0)
        self.send_command(24)
        
    def open_shutter(self):
        self.move_filter_wheel(self.selected_filter)

    def close_shutter(self):
        self.move_filter_wheel(6)
        
    def read_wavelength(self):
        #print "reading wavelength now :"
        self.send_command(24)
        nchr = self.ser.inWaiting()
        if ( nchr > 0):
            #print nchr, 'Hi!'
            junk = self.read_port()
            #print junk
        b = ''
        while (len(b) <= 1):
            self.ser.write(chr(29))
            time.sleep(0.25)
            b = self.read_port()
            #print b, len(b)
        wavelength = self.decode_wavelength(b[1:4])
        time.sleep(2.0)
        return wavelength

    def close(self):
        self.ser.close()
