#!/usr/bin/env python3
import logging

def log_init(rootLogger, logFile=None, verbose=False, quiet=False):
    rootLogger.setLevel(logging.DEBUG)
    logFormatter = logging.Formatter("%(asctime)s [%(levelname)s]  %(message)s")

    numeric_level = 20 # logging.INFOa
    if  verbose is True:
        numeric_level = 10  #logging.DEBUG
    if logFile is not None:
        fileHandler = logging.FileHandler("{0}".format(logFile))
        fileHandler.setLevel(numeric_level)
        fileHandler.setFormatter(logFormatter)
        rootLogger.addHandler(fileHandler)

    if quiet is True:
        numeric_level = 40 #logging.ERROR
    consoleHandler = logging.StreamHandler()
    consoleHandler.setLevel(numeric_level)
    consoleHandler.setFormatter(logFormatter)
    def decorate_emit(fn):
        # add methods we need to the class
        def new(*args):
            levelno = args[0].levelno
            if(levelno >= logging.CRITICAL):
                color = '\x1b[31;1m'
            elif(levelno >= logging.ERROR):
                color = '\x1b[31;1m'
            elif(levelno >= logging.WARNING):
                color = '\x1b[33;1m'
            elif(levelno >= logging.INFO):
                color = '\x1b[32;1m'
            elif(levelno >= logging.DEBUG):
                color = '\x1b[35;1m'
            else:
                color = '\x1b[0m'
            # add colored *** in the beginning of the message
            args[0].levelname = "{0}{1}\x1b[0m".format(color, args[0].levelname)
            #args[0].msg = "{0}{1}\x1b[0m".format(color, args[0].msg)
            return fn(*args)
        return new
    consoleHandler.emit = decorate_emit(consoleHandler.emit)

    rootLogger.addHandler(consoleHandler)
    return rootLogger

