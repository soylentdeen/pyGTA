import numpy
import smtplib
from email.mime.text import MIMEText
#import oauthtwitter
#import oauth
import twitter
import time
import matplotlib.pyplot as pyplot
from pylab import rand

def save_figure(xpts, ypts, dypts, names, title, xlabel, ylabel, filename):
    fig_width_pt = 246.0
    fig_width = 7.0
    inches_per_pt = 1.0/72.27
    pts_per_inch = 72.27
    golden_mean = (numpy.sqrt(5)-1.0)/2.0
    #fig_width = fig_width_pt*inches_per_pt
    fig_height = fig_width*golden_mean
    fig_size = [fig_width, fig_height]
    params = {'backend' : 'AGG',
              'axes.labelsize' : 12,
              'text.fontize' : 12,
              'legend.fontsize' : 12,
              'xtick.labelsize' : 10,
              'ytick.labelsize' : 10,
              'text.usetex' : True,
              'figure.figsize': fig_size}

    pyplot.rcParams.update(params)
    figure = pyplot.figure(0)
    figure.clf()
    axes = figure.add_subplot(1, 1, 1)
    #print "xpts = ", xpts
    #print "ypts = ", ypts
    #print "dypts = ", dypts
    #print "Names = ", names
    plots = []
    for x, y, dy, name in zip(xpts, ypts, dypts, names):
        a,b,c = axes.errorbar(x, y, dy, fmt='-', color = rand(100), label = name)
        plots.append(a)
    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)
    axes.set_title(title)
    figure.legend(plots, names, 'upper right')
    #pyplot.show()
    #print "Done!"
    pyplot.savefig(filename)

class Mailer( object ):
    def __init__(self, server, SENDER, RECIPIENTS):
        self.server = server
        self.SENDER = SENDER
        self.RECIPIENTS = RECIPIENTS

    def send_message(self, mesg):
        self.session = smtplib.SMTP(self.server)
        message = MIMEText(mesg)
        message['Subject'] = 'Message from Skynet'
        message['From'] = self.SENDER
        message['To'] = self.RECIPIENTS[0]
        smtpresult = self.session.sendmail(self.SENDER, self.RECIPIENTS, message.as_string())
        
    def close(self):
        self.session.quit()

class Personality( object ):
    def __init__(self, name):
        self.name = name
        self.exclamations = ['Zounds!', 'Bully!', 'I say!', 'Capital!', 'Huzzah!', 'Gadzooks!']
        self.oaths = ['Oh bother!', 'Rubbish!', 'Hogwash!', 'Horsefeathers!', 'What tomfoolery!']

    def exclamation(self):
        retval = ''


class Twitterer( object ):
    def __init__(self):
        self.consumer_key = 'W4Und2OJ4oDIlsgCTi5YGw'
        self.consumer_secret = 'VIcPs145UPgrjLQ90fUFON567S3vQWU3U5iqoy9vc'
        self.token_key = '182978048-nP508mwNfhj7tL6DJHzahTWyGUWQ5gxlEo69lviU'
        self.token_secret = 'STOYhRyXRp8FaVnwZYWlQtgH7N9iiuYGAwBeFo2iiE'
        #self.token = oauth.oauth.OAuthToken(self.token_key, self.token_secret)
        self.twit = twitter.Api(consumer_key=self.consumer_key, consumer_secret=self.consumer_secret, access_token_key=self.token_key, access_token_secret=self.token_secret)
        self.user = self.twit.GetUser('rowlands_ghost')
        self.last_update = time.time()
        

    def tweet(self, message):
        old_status = self.user.GetStatus()
        if (old_status == message):
            message = message + '_'
        self.user.SetStatus(message)
        try:
            new_time = time.time()
            while ((new_time - self.last_update) < 10):
                print 'Sleeping!'
                time.sleep(2)
                new_time = time.time()
            self.last_update = time.time()
            print self.last_update
            status = self.twit.PostUpdate(self.user.GetStatus())
            retval = status.text
        except:
            retval = 'Oh bother... the twittering mechanism seems to be down!'
        return retval
