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


my_log_format = "%(levelname)s: %(asctime)s %(funcName)s L%(lineno)s| %(message)s"



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
    logging.basicConfig(level=levels[level], format=my_log_format, datefmt="%Y/%m/%d %I:%M:%S")
    root_logger = logging.getLogger()

    
    for name in logging.Logger.manager.loggerDict.keys():
        if ('boto' in name) or ('urllib3' in name) or ('s3' in name) or ('findfont' in name):
            logging.getLogger(name).setLevel(logging.WARNING)

    
            
    return root_logger



class Loggah:
    
    def __init__(self, logfile:str=None, level:int=None):

        if logfile is None and type(logfile) is not str:
            raise TypeError("kmerdb.logger.Loggah() requires a str as its argument")

        elif level is None or type(level) is not int:
            raise TypeError("kmerdb.logger.Loggah() requires a int level")

        self._levels = ["DEBUG", "INFO", "WARNING", "ERROR"]
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
            fh.setFormatter(my_log_format)
            fh.setLevel(logging.DEBUG)
            self._logger.addHandler(fh)
            
        
        
        
    def log_it(self, log_str:str, level:str="WARNING"):
        """
        NOTE: All formatting and type casting is out of scope of this function.


        This function *may* coerce simple dictionaries and lists to valid strings.


        The Python logging module outputs to stderr as well as to the log-file.
        """
        if level is not None and (type(level) is not str or level not in self._levels):
            raise TypeError("kmerdb.logger.Loggah.log_it: invalid log level '{0}'".format(level))
        elif log_str is not None and (type(log_str) is not str and type(log_str) is not dict and type(log_str) is not list):
            raise TypeError("kmerdb.logger.Loggah.log_it: Unhashable type '{0}' supplied to logging function. Cannot coerce to a loggable string...")

        # Try to coerce dict/list to str
        if type(log_str) is dict or type(log_str) is list:
            try:
                yaml.add_representer(OrderedDict, util.represent_yaml_from_collections_dot_OrderedDict)
                log_str = yaml.dump(log_str)
            except Exception as e:

                self.logger.error("Could not convert input to a loggable ASCII representation")
                self.logger.error(log_str)
                raise e
            

        if level == "DEBUG":
            self._logger.debug(log_str)
        elif level == "INFO":
            self._logger.info(log_str)
        elif level == "WARNING":
            self._logger.warning(log_str)
        elif level == "ERROR":
            self._logger.error(log_str)

        self.logs.append(log_str)
