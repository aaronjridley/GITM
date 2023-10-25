#!/usr/bin/env python

import argparse
import os
from glob import glob
import time

IsVerbose = True
DoRm = True

# ----------------------------------------------------------------------
# Function to parse input arguments
# ----------------------------------------------------------------------

def parse_args_post():

    parser = argparse.ArgumentParser(description =
                                     'Post process and move model results')
    parser.add_argument('-file',
                        help = 'File that contains info. about remote system',
                        default = 'remote')

    parser.add_argument('-user',
                        help = 'remote user name (default none)',
                        default = 'none')
    
    parser.add_argument('-server',
                        help = 'remote system name (default none)',
                        default = 'none')
    
    parser.add_argument('-dir',
                        help = 'remote directory to use',
                        default = 'none')
    
    parser.add_argument('-sleep',
                        help = 'how long to sleep between loops',
                        default = 300, type = int)

    parser.add_argument('-q',
                        help = 'Run with verbose turned off',
                        action = 'store_true')
    
    parser.add_argument('-norm',
                        help = "don't remove any files",
                        action = 'store_true')
    
    parser.add_argument('-tgz',
                        help = "tar and zip raw GITM file instead of process",
                        action = 'store_true')
    
    args = parser.parse_args()

    return args

# ----------------------------------------------------------------------
# parse remote file
#   - remote file has the form:
# username
# remote server
# remote directory
# ----------------------------------------------------------------------

def parse_remote_file(file):

    if (IsVerbose):
        print('Reading file ', file)
    
    fpin = open(file, 'r')
    user = fpin.readline()
    server = fpin.readline()
    dir = fpin.readline()
    fpin.close()

    remote = {'user': user.strip(),
              'server': server.strip(),
              'dir': dir.strip()}
    return remote

# ----------------------------------------------------------------------
# do system command
# ----------------------------------------------------------------------

def run_command(command):
    if (IsVerbose):
        print("   -> Running Command : ")
        print("      ", command)
    os.system(command)
    return True

# ----------------------------------------------------------------------
# Check inputs:
# ----------------------------------------------------------------------

def check_inputs(user, server, dir):
    
    IsRemote = True
    if ((len(user) == 0) or (user == 'none')):
        print("Can't parse user information")
        IsRemote = False
    if ((len(server) == 0) or (server == 'none')):
        print("Can't parse server information")
        IsRemote = False
    if ((len(dir) == 0) or (dir == 'none')):
        print("Can't parse dir information")
        IsRemote = False

    return IsRemote

# ----------------------------------------------------------------------
# Check the test file to see if it is good
# ----------------------------------------------------------------------

def parse_test_file(file):

    if (os.path.exists(file)):
        IsGood = True
        if (IsVerbose):
            print('   --> Reading test file ', file)
        fpin = open(file, 'r')
        for line in fpin:
            if (line.find('No such') >= 0):
                IsGood = False
            if (line.find('Connection reset') >= 0):
                IsGood = False
        fpin.close()
    else:
        IsGood = False

    return IsGood
    
# ----------------------------------------------------------------------
# Checks to see if remote file or directory exists
# ----------------------------------------------------------------------

def test_if_remote_exists(user, server, dir):

    DidWork = True
    
    remote_command = user + "@" + server + " 'ls " + dir + "'"
    # check to see if the remote directory exists:

    if (IsVerbose):
        print('Checking to see if remote file or directory exists')
    command = 'ssh ' + remote_command + ' >& .test_file'
    DidWork = run_command(command)
    DidWork = parse_test_file('.test_file')

    if (DidWork):
        if (IsVerbose):
            print('   --> Remote directory (or file) exists!')
    else:
        print('--> Remote directory (or file) does NOT exist!')
        print('    Need to make this directory!')
        
    return DidWork

# ----------------------------------------------------------------------
# Make a remote directory
# ----------------------------------------------------------------------

def make_remote_dir(user, server, dir):

    DidWork = True
    
    remote_command = user + "@" + server + " 'mkdir " + dir + "'"
    # check to see if the remote directory exists:

    print('Making remote directory : ', dir)
    command = 'ssh ' + remote_command + ' >& .test_file'
    DidWork = run_command(command)
    DidWork = parse_test_file('.mkdir_command')

    return DidWork

# ----------------------------------------------------------------------
# Push log and data files to remote server
# ----------------------------------------------------------------------

def transfer_log_files(user, server, dir):

    remote = user + '@' + server + ':' + dir
    files = 'job* log.* UAM* *.txt *.dat*'
    outfile = '.output_rsync_log'
    command = 'rsync -vrae ssh ' + files + ' ' + remote + ' >& ' + outfile
    DidWork = run_command(command)

    return DidWork
    
# ----------------------------------------------------------------------
# Determine base filenames from header file
# ----------------------------------------------------------------------

def determine_base(headerFile, dir):

    start = len(dir)
    end = len(headerFile) - len('.header')
    
    baseFile = headerFile[start:end]

    return baseFile
    
# ----------------------------------------------------------------------
# tar and zip raw GITM files
# ----------------------------------------------------------------------

def tar_and_zip_gitm():

    print('-> Tar and Zipping files!')
    data_here = 'UA/data'
    
    # get header file list:
    headerFiles = sorted(glob(data_here + '/*.header'))

    for headerFile in headerFiles:
        baseFile = determine_base(headerFile, data_here + '/')
        print('--> Processing header file : ', baseFile)

        tarFile = baseFile + '.tgz'
        command = 'cd ' + data_here + ' ; rm -f ' + tarFile + ' ; cd ../..'
        DidWork = run_command(command)
        
        command = 'cd ' + data_here + ' ; tar -cvzf ' + tarFile + ' ' + baseFile + '.* ; cd ../..'
        DidWork = run_command(command)

        # Remove raw files:
        command = 'cd ' + data_here + ' ; rm -f ' + baseFile + '.b[0-9][0-9][0-9][0-9] ; cd ../..'
        DidWork = run_command(command)
        command = 'cd ' + data_here + ' ; rm -f ' + baseFile + '.header ; cd ../..'
        DidWork = run_command(command)
        command = 'cd ' + data_here + ' ; rm -f ' + baseFile + '.sat ; cd ../..'
        DidWork = run_command(command)

    DidWork = True
        
    return DidWork

# ----------------------------------------------------------------------
# post process GITM files
# ----------------------------------------------------------------------

def post_process_gitm():

    command = \
        'cd UA ; ' + \
        'chmod a+rx data ; '
    if (IsVerbose):
        command = command + './pGITM ; '
    else:
        command = command + './pGITM > .post_process ; '
    command = command + 'chmod a+r data/*.bin ; cd ..'
    DidWork = run_command(command)

    return DidWork

# ----------------------------------------------------------------------
# Transfer file, check if it made it, then delete local (if requested)
# ----------------------------------------------------------------------

def transfer(filelist, user, server, dir, DoRemove):

    DidWork = True
    remote = user + '@' + server + ':' + dir

    files = ''
    outfile = '.output_rsync_log'

    if (len(filelist) > 0):
        for file in filelist:
            chmod = 'chmod a+r ' + file
            DidWork = run_command(chmod)
            files = files + ' ' + file
            
        rsync = 'rsync -rav ' + files + ' ' + remote
        if (not IsVerbose):
            rsync = rsync + ' >> ' + outfile + ' 2>&1'
        DidWork = run_command(rsync)
        
    if (DoRemove):
        for file in filelist:
            sep = file.split('/')
            test_file = sep[-1]
            DidTransfer = test_if_remote_exists(user,
                                                server,
                                                dir + '/' + test_file)
            if (DidTransfer):
                if (IsVerbose):
                    print('   --> Remote file (' + test_file + ') exists!' +
                          '  Deleting local!')
                if (DoRm):
                    command = '/bin/rm -f ' + file
                    DidWork = run_command(command)
                    # Systems reject ssh commands if too many happen in
                    # too short of time, so sleep in between the commands. 
                    time.sleep(5)
            else:
                if (IsVerbose):
                    print('Remote file (' + test_file + ') does not exist!')

    return DidWork

# ----------------------------------------------------------------------
# Transfer a bunch of different file types to remote machine:
# ----------------------------------------------------------------------

def transfer_model_output_files(user, server, dir):

    data_here = 'UA/data'
    remote = user + '@' + server + ':' + dir

    print("Transfering files...")
    
    # run info file:
    file = 'run_information.txt'
    filelist = sorted(glob(data_here + '/' + file))
    DoRemove = True
    DidWork = transfer(filelist, user, server, dir, DoRemove)

    file = 'log*.dat'
    filelist = sorted(glob(data_here + '/' + file))
    DoRemove = False
    DidWork = transfer(filelist, user, server, dir, DoRemove)

    # data files - remove by default
    DoRemove = True
    file = '3D*.bin'
    filelist = sorted(glob(data_here + '/' + file))
    DidWork = transfer(filelist, user, server, dir, DoRemove)

    file = '2D*.bin'
    filelist = sorted(glob(data_here + '/' + file))
    DidWork = transfer(filelist, user, server, dir, DoRemove)

    file = '1D*.bin'
    filelist = sorted(glob(data_here + '/' + file))
    DidWork = transfer(filelist, user, server, dir, DoRemove)

    file = '0D*.bin'
    filelist = sorted(glob(data_here + '/' + file))
    DidWork = transfer(filelist, user, server, dir, DoRemove)

    # tar and zip - remove by default
    DoRemove = True
    file = '3D*.tgz'
    filelist = sorted(glob(data_here + '/' + file))
    DidWork = transfer(filelist, user, server, dir, DoRemove)

    file = '2D*.tgz'
    filelist = sorted(glob(data_here + '/' + file))
    DidWork = transfer(filelist, user, server, dir, DoRemove)

    file = '1D*.tgz'
    filelist = sorted(glob(data_here + '/' + file))
    DidWork = transfer(filelist, user, server, dir, DoRemove)

    file = '0D*.tgz'
    filelist = sorted(glob(data_here + '/' + file))
    DidWork = transfer(filelist, user, server, dir, DoRemove)

    return DidWork

# ----------------------------------------------------------------------
# Post process and then transfer files once:
# ----------------------------------------------------------------------

def do_loop(doTarZip, user, server, dir, IsRemote):

    DidWork = True
    
    if (IsRemote):
        DidWork = test_if_remote_exists(user, server, dir)
        if (not DidWork):
            return DidWork

    if (IsRemote and DidWork):
        DidWork = transfer_log_files(user, server, dir)
        if (not DidWork):
            return DidWork

    # Post process GITM files:
    print('Post Processing GITM files...')
    if (doTarZip):
        DidWork = tar_and_zip_gitm()
    else:
        DidWork = post_process_gitm()
        
    # Check if remote data directory exists, make it if it doesn't:
    data_remote = '/data'
    if (IsRemote and DidWork):
        dir = dir + data_remote
        DidWork = test_if_remote_exists(user, server, dir)
        if (not DidWork):
            DidWork = make_remote_dir(user, server, dir)
            
    DidWork = transfer_model_output_files(user, server, dir)
    
    return DidWork

# ----------------------------------------------------------------------
# Main Code
# ----------------------------------------------------------------------

if __name__ == '__main__':  # main code block

    args = parse_args_post()

    doTarZip = args.tgz

    file = args.file

    IsVerbose = not args.q
    if (args.norm):
        DoRm = False

    if (os.path.exists(file)):
        print('Found file: ', file)
        remote = parse_remote_file(file)
        user = remote['user']
        server = remote['server']
        dir = remote['dir']
    else:
        user = args.user
        server = args.server
        dir = args.dir

    # Check inputs:
    IsRemote = check_inputs(user, server, dir)

    DidWork = True
    
    while DidWork:
        DidWork = do_loop(doTarZip, user, server, dir, IsRemote)
        if (DidWork):
            print('Sleeping ... ', args.sleep, ' sec.')
            time.sleep(args.sleep)
