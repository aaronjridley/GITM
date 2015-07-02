#!/usr/bin/perl
# Usage: 
# Edit this file to set things and then run it as
#
# work.pl >& work.log &
#

my $Success1 = "SWMF.SUCCESS";
my $Success2 = "BATSRUS.SUCCESS";

print "running on ", `hostname`;
while(not -d "RESTART_t0003.0000h"){      # final restart file
    unlink $Success1, $Success2;          # remove success indicator
    `qsub job.devel`;                     # submit job
    print "job submitted on ", `date`;
    sleep 7000;                           # no need to check yet
    sleep 5 while not (-f $Success1 or -f $Success2); 
    print "job finished successfully on ", `date`;
    sleep 60;                             # wait for the job to quit
}
print "job reached final time on ", `date`;
exit 0;
