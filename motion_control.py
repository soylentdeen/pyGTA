import socket
import sys

class Motion_Controller( object ):

    def __init__(self, parameters):
        self.host = parameters['MOTOR_HOST']
        self.port = int(parameters['MOTOR_PORT'])
        self.sock = None
        self.gearing_ratio = parameters['GEARING_RATIO']
        self.microsteps = parameters['MICROSTEPS'] # number of steps per revolution of motor
        self.home_pos_angle = parameters['HOME_POS_ANGLE']   # Home position in degrees away from zero
        self.home_pos_steps = parameters['HOME_POS_STEPS']

        for res in socket.getaddrinfo(self.host, self.port, socket.AF_UNSPEC, socket.SOCK_STREAM):
            af, socktype, proto, canonname, sa = res
            try:
                self.sock = socket.socket(af, socktype, proto)
            except socket.error, msg:
                self.sock = None
                continue
            try:
                self.sock.connect(sa)
            except socket.error, msg:
                self.sock.close()
                self.sock = None
                continue
            break
        if self.sock is None:
            print 'Error!  Could not open the socket for the motor control!'
            sys.exit(1)

    def move_abs(self, position):
        self.sock.send('MOVE_ABS '+str(position)+'\r\n')
        data = self.sock.recv(1024)
        return 1

    def move_rel(self, position):
        self.sock.send('MOVE_REL '+str(position)+'\r\n')
        data = self.sock.recv(1024)
        return 1

    def goto(self, angle):
        delta_angle = self.home_pos_angle - angle
        position = self.home_pos_steps + self.microsteps*(self.gearing_ratio/360.0)*delta_angle
        retval = self.move_abs(position)
        return retval

    def read_pos(self):
        self.sock.send('READ_POS\r\n')
        return self.sock.recv(1024)

    def go_home(self):
        self.sock.send('HOME\r\n')
        retval = self.sock.recv(1024)
        self.home_pos_steps = float(self.read_pos())
        return 1

    def close(self):
        self.sock.close()
        
    def __del__(self):
        self.sock.close()
