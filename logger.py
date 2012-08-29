##############################################################################
#
# Logger class:
# a class for writing a log file of the Lazarus run.
#
#
# Author: Victor Hanson-Smith
# Contact: victorhs@cs.uoregon.edu
#
###############################################################################

import os, time

class Logger:    
    def __init__(self, verbosityLevel, outputDirectory):
        self.verbosityLevel = verbosityLevel
        if self.verbosityLevel == False:
            self.verbosityLevel = 1
        
        if os.path.exists(outputDirectory) == False:
            os.system("mkdir " + outputDirectory)
            
        self.outputPath = outputDirectory + "/lazarus_job_status.log"
        f = open(self.outputPath, "w")
        f.close()
    
    #
    # t is the array of integers returned by time.localtime()
    #
    def localtimeToString(self, t):
        r = "(" + t[1].__str__() + "/" + t[2].__str__() + "/" + t[0].__str__() + " " + t[3].__str__() + ":" + t[4].__str__() + ":" + t[5].__str__() + ")"
        return r
    
    #
    # The status file contains a log of status updates.
    # Each update is two lines long: one line for the timestamp, and one line for the message.
    #
    def updateStatus(self, message):
        f = open(self.outputPath, "a")
        t = time.localtime()
        f.write( self.localtimeToString(t) + " " + message + "\n")
        f.close()
    #
    # display an error message to the screen and exit (i.e. kill) the process
    #
    def throwError(self, message):
        print "Hmmm, " + message
        self.updateStatus("An error occurred: " + message)
        raise AssertionError
    
    #
    # prints a message, assuming the verbosity level ("vl") of the message is equal to or below
    # the verbosity level specified by the user
    #    
    def printMessage(self, message, vl):
        if vl <= self.verbosityLevel:
            print message