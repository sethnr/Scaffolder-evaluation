#!/usr/bin/perl -ws

#GET VARIABLES
# no of cores, output files

#TEMPLATES




$template = <<<"END";



END;
my $template_lg = $template.format($rundist);
my $template_sm = $template.format($runlocal);

my $bsub_script = "...bsub.pl";

$tempfile_lg.format($reads,
                    $contigs,
                    $outfolder
                   );

$tempfile_sm.format($reads,
                    $contigs,
                    $outfolder
                   );
my $LASTID = -1;

#def submit single-job job
sub _submitLocalJob($command, $pre, ...) {
  # bsub, write to log file...
  # sub it, grab ID,
  $LASTID = $last_job_id;
  # make next job dependent on completion if ID
}


#def submit multi-job job
sub _submitLocalJob($command, $pre, ...) {
  # submit job, get Job ID,
  # wait, check, wait...
  # when complete
  # read out list of subJobs from Array
  # return array of job IDs to main runner
}


#SETUP
my $run_index = rand(...);
my $run_folder = "tmp.pbjelly.".$run_index;
system("mkdir $run_folder");
system("cd $run_folder");

open(TEMPLG,">".$run_folder."/template_large.xml");
  print TEMPLG $tempfile_lg;
close(TEMPLG)
  
open(TEMPSM,">".$run_folder."/template_small.xml")
  print TEMPSM $tempfile_sm;
close(TEMPSM)


# make suffixarray...
  
#MAPPING (multi)
# check exists m4 / bml?

#SUPPORT (single)
# check bml file?
  
#EXTRACT (single)
#ASSEMBLY (multi)
#OUTPUT (single)


#clean up (template file, log files, etc)
