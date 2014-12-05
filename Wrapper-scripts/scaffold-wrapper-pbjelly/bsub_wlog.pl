#!/usr/bin/perl -ws

my $logfile = shift;
#my $parentID = shift; #GET FROM LSF ENV
open(LOG,">>$logfile");

my $command = "bsub "+join(" ",@ARGV);
print STDERR $command;
my $output = `$command`;
print STDERR $output;
if($output =~ m/\<(\d+)\>.*submitted*/gi) {
  $jobID = $1;
  print LOG $parentID."\t".$jobID;
  close(LOG);
  exit(0);
}
else {
  print LOG $output;
  exit(1001);
}
