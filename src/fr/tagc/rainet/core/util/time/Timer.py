import time

from fr.tagc.rainet.core.util.log.Logger import Logger


##
# This class permits to easy mark the various durations of workflow step
class Timer( object):
    
    ## The singleton instance
    __instance = None

    ##
    # The constructor
    def __init__(self):
        
        self.startTime = 0
        self.lastTime = 0
        self.currentStep = None

    ##
    # The singleton provider
    @staticmethod
    def get_instance():
        if Timer.__instance == None:
            Timer.__instance = Timer()
        return Timer.__instance

    ##
    # This method permits to initialize the chrono, displaying th
    # provided message
    #
    # @param message : string - The message to display
    def start_chrono(self):
        
        self.start_time = time.time()
        self.lastTime = self.start_time
        Logger.get_instance().info ( "\nSTART CHRONO\n")
        self.currentStep = None
    
    ##
    # This method permits to get the current duration from the last chrono start
    # without reinitializing the chrono
    def time_point(self):
        
        current_time = time.time()
        duration = current_time - self.lastTime
        Logger.get_instance().info ( "Current duration : " + Timer.format_duration( duration))
    
    ##
    # This method permits to get the stop the chrono, reinitializing it and
    # displaying the provided message
    #
    # @param message : string - The message to display
    def stop_chrono(self, message):
        
        current_time = time.time()
        step_duration = current_time - self.lastTime
        total_duration = current_time - self.start_time
        Logger.get_instance().info ( "Step duration : " + Timer.format_duration( step_duration))
        Logger.get_instance().info ( "\n\nSTOP CHRONO : " + message + ". Total duration " + Timer.format_duration(total_duration))
        self.lastTime = 0
        self.start_time = 0
    
    ##
    # This methods permit to indicate the duration from the last chrono start
    # and to start a new chrono at the same time with the provided message
    #
    # @param message : string - The message to display
    def step(self, message):
        
        current_time = time.time()
        duration = current_time - self.lastTime
        if self.currentStep != None:
            Logger.get_instance().info ( "Step duration : " + Timer.format_duration( duration) + "\n")
        self.lastTime = current_time
        Logger.get_instance().info ( "\n------------------------------------------------")
        Logger.get_instance().info ( "START STEP : '" + message+ "'")
        self.currentStep = message
        
    ##
    # This method provide a human readable version of the duration
    @staticmethod
    def format_duration(real_seconds):
        seconds = long(round(real_seconds))
        minutes, seconds = divmod(seconds, 60)
        hours, minutes = divmod(minutes, 60)
        days, hours = divmod(hours, 24)
        years, days = divmod(days, 365.242199)
     
        minutes = long(minutes)
        hours = long(hours)
        days = long(days)
        years = long(years)
     
        duration = []
        if years > 0:
            duration.append('%d year' % years + 's'*(years != 1))
        else:
            if days > 0:
                duration.append('%d day' % days + 's'*(days != 1))
            if hours > 0:
                duration.append('%d hour' % hours + 's'*(hours != 1))
            if minutes > 0:
                duration.append('%d minute' % minutes + 's'*(minutes != 1))
            if seconds > 0:
                duration.append('%d second' % seconds + 's'*(seconds != 1))
            if seconds == 0:
                duration.append( "less than 1 second")
        return ' '.join(duration)
        
        

