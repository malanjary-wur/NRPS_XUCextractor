#!/usr/bin/env python
# 2019 Mohammad Alanjary
# Wageningen University and Research (WUR)
# Bioinformatics department, Lab of Prof. Marnix Medema
#
# License: A copy of the GPLv3 can also be found at: <http://www.gnu.org/licenses/>. 

import logging

def init(logfile=None,toconsole=False,level="debug"):
    log=logging.getLogger("root")
    formatter = logging.Formatter(fmt='%(asctime)s - %(levelname)s - %(module)s - %(message)s')
    if logfile:
        log.handlers = [] #clear all handlers
        handler=logging.FileHandler(logfile)
        handler.setFormatter(formatter)
        log.addHandler(handler)
    if toconsole:
        log.handlers = [] #clear all handlers
        handler=logging.StreamHandler()
        handler.setFormatter(formatter)
        log.addHandler(handler)
    if level.upper()=="CRITICAL":
        log.setLevel(logging.CRITICAL)
    elif level.upper()=="ERROR":
        log.setLevel(logging.ERROR)
    elif level.upper()=="WARNING":
        log.setLevel(logging.WARNING)
    elif level.upper()=="INFO":
        log.setLevel(logging.INFO)
    elif level.upper()=="DEBUG":
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.NOTSET)
    return log