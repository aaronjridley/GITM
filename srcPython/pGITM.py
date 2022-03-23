#!/usr/bin/env python3
'''
A wrapper for PostGITM.exe to handle all unconcatenated output fragments 
produced during a GITM simulation.

This script is an alternative to pGITM (written in cshell) that works robustly
and securely on a range of unix-like environments with Python available
(i.e., any modern computer.)
'''

import os
from glob import glob
from subprocess import run
from argparse import ArgumentParser, RawDescriptionHelpFormatter

# Create object to handle command line arguments.
# Right now, only argument is "-h" (default help.)
parser = ArgumentParser(description=__doc__,
                        formatter_class=RawDescriptionHelpFormatter)
args = parser.parse_args()

# Move into the UA/data directory (required by PostGITM.exe)
os.chdir('UA/data')

# Get list of header files to process:
headers = glob('*.header')

# For each file, pipe file name to PostGITM and process away.
for head in headers:
    # Attempt to process file:
    print(f'Processing {head}...')
    process = run('../../PostGITM.exe', input=head, text=True)

    # Handle result robustly:
    if process.returncode:
        # On failure, warn user.
        print(f'WARNING: Error processing {head} (Code {process.returncode}')
    else:
        # On success, remove sub-files
        stem = head.replace('.header','')
        os.remove(head) # remove header
        for f in glob(stem+'.b????'):
            os.remove(f) # remove fragments
        if os.path.exists(stem+'.sat'):
            os.remove(stem+'.sat') # remove sat file.

    
    
