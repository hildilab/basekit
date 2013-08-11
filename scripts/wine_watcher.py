#! /usr/bin/env python
# initial idea by Dominic Theune


import os
import subprocess
from time import sleep, strftime


def get_active_windbg():
    output = subprocess.check_output([ "ps", "-ef" ]).split("\n")
    out2 = [ l for l in output if l.find("winedbg")!=-1 ]
    pids = []
    for l in map( lambda x: x.split(), out2 ):
        if l[7]=='winedbg':
            print l
            pids.append( l[1] )
    return pids

def kill_pids( pids ):
    for pid in pids:
        try:
            os.system( "kill -9 %s" % pid )
        except:
            print pid, "kill failed"


while True:
    print strftime('%X %x %Z')
    pids = get_active_windbg()
    kill_pids( pids )
    sleep(10)
