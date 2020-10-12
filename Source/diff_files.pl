#!/usr/local/bin/perl

# File: diff_files.pl
#
# For each source code file in the current directory, 
# diff with the corresponding file in a second directory.
#
# Usage: perl diff_files.pl
#
# Melannie Hartman
# September 5, 2018
#

$dir1 = "./";
#$dir2 = "/data/parton/melannie/GrassTree/ForCABBI_02.01.2018/Current_Files/DayCent_Photosyn_grasstree/";
$dir2 = "./SOURCE_Nov_2017_Elena/";

@srcfiles = glob("*.f *.c *.inc *.h");
foreach $srcfile (@srcfiles)
{
    print $srcfile, "\n";
    $srcfile1 = "${dir1}${srcfile}";
    $srcfile2 = "${dir2}${srcfile}";
    $cmdline = "diff ${srcfile1} ${srcfile2} > Diff/diff_${srcfile}.txt";
    print "$cmdline\n";
    system($cmdline);
}

