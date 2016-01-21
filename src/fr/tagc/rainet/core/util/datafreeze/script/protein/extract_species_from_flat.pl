#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
my %opts;
getopts('s:',\%opts);

my $species=$opts{s}||die("Need a species (-s)!\n");
# Change the line termination string so we read an entire entry at a time
local $/ = "\n//\n";


while(<>){
    next unless /OS\s+.*$species.*/i;
    print;
}
