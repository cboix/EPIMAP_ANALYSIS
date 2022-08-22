#!/usr/bin/python
# ---------------------------------
# Parse the bash config for python
# Then read in exported environment
# ---------------------------------
import os
import pprint
import shlex
import subprocess
import socket
import re

# Check domain and set home and directories::
domain = socket.getfqdn()
if 'broadinstitute.org' in domain:
    homedir = '/broad/compbio/cboix'
else:
    homedir = os.getenv("HOME")

bindir = homedir + "/EPIMAP_ANALYSIS/bin/"

configfile = bindir + 'config_ChromImpute.sh'
command = shlex.split("env -i bash -c 'source " +
                      configfile + " && env'")
proc = subprocess.Popen(command, stdout = subprocess.PIPE)
for line in proc.stdout:
    line = line.decode("utf-8")
    line = line.split("\n")[0]
    if re.match('.*=.*', line):
        (key, _, value) = line.partition("=")
        os.environ[key] = value

# Check if anything left in buffer:
# proc.communicate()
# pprint.pprint(dict(os.environ))

print(os.environ['BINDIR'])
print("[STATUS] Loaded config correctly.")
