#!/usr/bin/perl -ws

my $logfile = shift;
#my $parentID = shift; #GET FROM LSF ENV
my $parentID = $ENV{"LSB_JOBID"};
my $parentName = $ENV{"LSB_JOBNAME"};
#my $logName = "bsub_log";
$logfile =~ s/.txt$//gi;
$logfile .= ".".$parentName if ($parentName ne "");
$logfile .= ".".$parentID.".txt";
open(LOG,">>$logfile");

my $command = "bsub -e ".$parentName."_".$parentID.".%J.err ".
    "-o ".$parentName."_".$parentID.".%J.out ".
    join(" ",@ARGV);

print STDERR $command;
my $output = `$command`;
print STDERR $output;
if($output =~ m/\<(\d+)\>.*submitted*/gi) {
  $jobID = $1;
  print LOG $parentID."\t".$jobID."\t".$command."\n";
  close(LOG);
  exit(0);
}
else {
  print LOG $output;
  exit(1001);
}
