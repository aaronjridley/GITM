#!/usr/bin/perl -pi~
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

# Convert source code to use BATL variables
# The script does not do a perfect job. 
# The resulting code should be checked and edited as needed.
#
# Usage: 
#       FixBatl.pl file1.f90 file2.f90 ...

# UnusedBLK --> Unused_B
s/\bunusedblk\b/Unused_B/gi;

# east_ --> 1
s/\beast_\b/1/gi;

# west_ --> 2
s/\bwest_\b/2/gi;

# south_ --> 3
s/\bsouth_\b/3/gi;

# north_ --> 4
s/\bnorth_\b/4/gi;

# bot_ --> 5
s/\bot_\b/5/gi;

# top_ --> 6
s/\btop_\b/6/gi;

# Blkneighborlev --> DiLevelNei_IIIB
s/\bblkneighborlev\b/DiLevelNei_IIIB/gi;

# 1-gcn:nI+gcn --> MinI:MaxI
# 1-gcn,nK+gcn --> MinK,MaxK
s/1\s*\-\s*gcn\s*([,:])\s*n([ijk])\s*\+\s*gcn/Min$2$1Max$2/gi;

# use ModSize, ONLY: ...gcn... --> 
# use ModSize, ONLY: ...MinI, MaxI, MinJ, MaxJ, MinK, MaxK
s/^(\s*use.*)\bgcn\b/$1MinI, MaxI, MinJ, MaxJ, MinK, MaxK/i;

# -1:nI+2 --> MinI:MaxI
# -1,nK+2 --> MinK,MaxK
s/\-\s*1\s*([,:])\s*n([ijk])\s*\+\s*2/Min$2$1Max$2/gi 
    unless /n[IJK]\s*\-\s*1/;

# x_BLK( --> Xyz_DGB(x_,
# y_BLK( --> Xyz_DGB(y_,
# z_BLK( --> Xyz_DGB(z_,
s/\b([xyz]_)blk\(/Xyz_DGB\($1,/gi;

# use ... x_BLK, y_BLK, z_BLK --> use ... Xyz_DGB
s/^(\s*use.*)\bx_blk\b(\s*,\s*y_blk)?(\s*,\s*z_blk)?/$1Xyz_DGB/i;

# dx_BLK( --> CellSize_DB(x_,
# dy_BLK( --> CellSize_DB(x_,
# dz_BLK( --> CellSize_DB(x_,
s/\bd([xyz]_)blk\(/CellSize_DB\($1,/gi;

# use ... dx_BLK, dy_BLK, dz_BLK --> use ... CellSize_DB
s/^(\s*use.*)\bdx_blk\b(\s*,\s*dy_blk)?(\s*,\s*dz_blk)?/$1CellSize_DB/i;
