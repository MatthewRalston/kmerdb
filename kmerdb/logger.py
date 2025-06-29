'''
   Copyright 2020 Matthew Ralston

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

'''




import sys
import os
import yaml, json


from collections import OrderedDict



from kmerdb import util


import logging
logger = logging.getLogger(__file__)


MAIN_LOG_FORMAT = "[%(levelname)s] %(asctime)s|%(filename)s %(funcName)s L%(lineno)s| %(message)s"
ABBREV_LOG_FORMAT = "[%(levelname)s] %(asctime)s| %(message)s"

LEVEL_STRS = ["DEBUG", "INFO", "WARNING", "WARN", "ERROR", "CRITICAL"]
LEVEL_NUMS = [10, 20, 30, 40, 50]

def get_root_logger(level:int, logfile:str=None):
    if level is not None and type(level) is not int:
        raise TypeError("get_root_logger requires the log-level to be an int")

    elif logfile is not None and type(logfile) is not str:
        raise TypeError("kmerdb.logger.get_root_logger requires the logfile to be a str")
    elif os.path.exists(logfile):
        sys.stderr.write("Existing log-file '{0}' found. WARNING! Overwriting\n".format(logfile))
    
    levels=[logging.WARNING, logging.INFO, logging.DEBUG]
    if level < 0 or level > 2:
        raise TypeError("{0}.get_root_logger expects a verbosity between 0-2".format(__file__))
    logging.basicConfig(level=levels[level], format=MAIN_LOG_FORMAT, datefmt="%Y/%m/%d %I:%M:%S")
    root_logger = logging.getLogger()

    
    for name in logging.Logger.manager.loggerDict.keys():
        if ('boto' in name) or ('urllib3' in name) or ('s3' in name) or ('findfont' in name):
            logging.getLogger(name).setLevel(logging.WARNING)

    
            
    return root_logger



class AppLogger:
    
    def __init__(self, logfile:str=None, level:int=None):

        if logfile is not None and type(logfile) is not str:
            raise TypeError("template_py.logger.AppLogger() requires a str as its argument")

        elif level is None or type(level) is not int:
            raise TypeError("template_py.logger.AppLogger() requires a verbosity level to control logging")

        self._levels = LEVEL_STRS
        self._levelnums = LEVEL_NUMS
        self._format = MAIN_LOG_FORMAT
        self._alt_format = ABBREV_LOG_FORMAT
        self.logs = []
        self.logfile = logfile
        self.level = level
        self._logger = get_root_logger(self.level, logfile=logfile)


        if os.path.exists(logfile):
            sys.stderr.write("Overwriting existing log file '{0}'...\n".format(logfile))
            pass
        elif not os.path.exists(os.path.dirname(logfile)):
            sys.stderr.write("Creating new log file '{0}'...\n".format(logfile))

        if logfile is not None:
            
            fh = logging.FileHandler(logfile, mode="w")
            fh.setFormatter(MAIN_LOG_FORMAT)
            fh.setLevel(logging.DEBUG)
            self._logger.addHandler(fh)
            
        if hasattr(sys, "_getframe"):
            self.currentframe = lambda: sys._getframe(2)
        
        
        
    def log_it(self, log_str:str, level:str="WARNING"):
        """
        NOTE: All formatting and type casting is out of scope of this function.


        This function *may* coerce simple dictionaries and lists to valid strings.


        The Python logging module outputs to stderr as well as to the log-file.
        """
        if (type(level) is str and level in LEVEL_STRS) or (type(level) is int and level in LEVEL_NUMS):
            pass
        else:
            raise ValueError("template_py.logger.AppLogger.log_it(): invalid log level '{0}'".format(level))
            
        if log_str is not None and (type(log_str) is not str and type(log_str) is not dict and type(log_str) is not list):
            raise TypeError("configurator.logger.AppLogger.log_it(): Unhashable type '{0}' supplied to logging function. Cannot coerce to a loggable string...")
        elif type(log_str) is not str:
            raise TypeError("template_py.logger.AppLogger.log_it(): Unknown log_str supplied to log_it")
        # Try to coerce dict/list to str
        if type(log_str) is dict or type(log_str) is list:
            try:
                yaml.add_representer(OrderedDict, util.represent_yaml_from_collections_dot_OrderedDict)
                log_str = yaml.dump(log_str)
            except Exception as e:

                self.logger.error("Could not convert input to a loggable ASCII representation")
                self.logger.error(log_str)
                raise e

        #for hndlr in 
        _frm = self.currentframe()
        _code = _frm.f_code
        _filename = os.path.basename(_code.co_filename)
        _funcname = _code.co_name
        _lineno = _frm.f_lineno

        # Append additional stack information to the front of the string.
        
        additl = f"|{_filename} {_funcname} L{_lineno}| "
        log_str = additl + log_str

        if level == 'DEBUG':
            self._logger.debug(log_str)
        elif level == 'INFO':
            self._logger.info(log_str)
        elif level == 'WARNING':
            self._logger.warning(log_str)
        elif level == 'ERROR':
            self._logger.error(log_str)
        elif level == 'CRITICAL':
            self._logger.critical(log_str)

        self.logs.append(log_str)

    def debug(self, log_str:str):
        self.log_it(log_str)
    def info(self, log_str:str):
        self._logger.info(log_str)
    def warn(self, log_str:str):
        self._logger.warning(log_str)
    def warning(self, log_str:str):
        self._logger.warning(log_str)
    def error(self, log_str:str):
        self._logger.error(log_str)
        
