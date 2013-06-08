#!/usr/bin/perl

#---------------------------------------------------------------------
# *** autocollect.pl ***
# A small script to collect the latest restarting files.
# 06/08/2013  Kengo TOMIDA (tomida__AT__astro.princeton.edu)
# Usage:
# > perl autocollect.pl [problemid]
# > mpirun /home/tomidakn/athena/bin/athena -r restart.latest.rst
# Note:
# This script can be used to continue your simulations automatically.
# This generates restart-id*.latest.rst files in the current directory
# without any confirmation. So please be careful.
#---------------------------------------------------------------------

use File::Copy;

my $pre=$ARGV[0];

opendir DIR, "id0";
my @list=readdir DIR;
closedir DIR;
my %lstep=();


my $max=0;
my $step="";
foreach $file (@list)
{
	$file=~/\.([\d]+)\.rst$/;
	my $tstep=$1;
	$tstep=~/0*([\d]+)/;
	my $num=$1;
	if($num > $max) {$max=$num; $step=$tstep;}
}

if(-e "restart.latest.rst")
{
  print "Cleaning up previous restart files\n";
  system("rm -rf restart.*.rst");
}
print "Collecting restart files,  step = $step\n";

opendir DIR, ".";
my @ids=readdir DIR;
closedir DIR;

foreach $id (@ids)
{
	if($id eq "id0")
	{
		copy("$id/$pre.$step.rst", "restart.latest.rst");
	}
	elsif($id=~/^id([\d]+)$/)
	{
		copy("$id/$pre-$id.$step.rst", "restart-$id.latest.rst");
	}
}
print "Ready to go!\n";


exit;


