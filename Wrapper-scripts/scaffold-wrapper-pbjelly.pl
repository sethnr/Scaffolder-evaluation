#!/usr/bin/perl -ws

use strict;
use warnings;
use Getopt::Long;
use File::Spec;
use File::Spec::Link;
use File::Which;
use File::Basename;
use FindBin;

my %options;
my $no_clean = -1;
my $debug = -1;
my $extras = "";
my($blasr_opts, $NJOBS, $no_cores, $dataset_name,
    $contigs, @reads, $outfile, $outdir);

#PBJELLY OPTIONS
$NJOBS = 10;


my $localmem = 8000;
my $distmem = 2000;

my $preexec = "~/Pacbio_tests/setPaths.sh";

#BLASR_OPTIONS
my %blasr = (
"minMatch", 8,
"sdpTupleSize", 8,
"minPctIdentity", 75,
"bestn", 1,
"nCandidates", 10,
"maxScore", 500, 
"nrpoc", 12
);


my $options_ok = GetOptions(\%blasr,
'name=s'    => \$dataset_name,
'noclean'   => \$no_clean,
'blasr=s'   => \$blasr_opts,
'jobs=i'    => \$NJOBS,
'cores=i'   => \$no_cores,
'debug'     => \$debug,

'reads|r=s'   => \@reads,
'contigs|c=s' => \$contigs,
'o=s'       => \$outfile

);

print join("\n",%blasr);

if (!($options_ok)) {
print STDERR "usage: $0 [options] -c <contigs.fa> -o <output dir> -r <reads.fq> [-r <reads.fq>]

Options:
-blasr_ops \"STRING\"
Put blasr options in quotes. Unless this option is used, defaults will be used.
-noclean
Don't clean up files. Default is to delete most files when finished
";
exit(1);
}


$contigs = File::Spec->rel2abs($contigs);

print STDERR "ARGV: ".$#ARGV."\n";
print STDERR "@ARGV";
print STDERR "\n";
print STDERR "SCRIPT: ".$0."\nSCRIPTDIR: ".dirname($0);

my $READSDIR = dirname($reads[0]);
$reads[0] = basename($reads[0]);

my $READS = "";
foreach my $read (@reads) {
    $READS .= "<job>".basename($read)."</job>\n";
}

#my $outdir = $ARGV[2];

#######################
# SETUP
#######################
my $run_index = int(rand(100000));
my $run_folder = File::Spec->rel2abs("tmp.pbjelly.".$run_index);

#print STDERR $run_folder."\n";
system("mkdir $run_folder"); # or die "could not make rundir";
chdir($run_folder) or die "could not cd into $run_folder";

print STDERR $run_folder."\n";

#######################
#TEMPLATES
#######################
$contigs = "/lustre/scratch108/parasites/snr/refs/Pf3D7.version3.uniquely-tagged.contigs.fasta";
$outdir = "/lustre/scratch108/parasites/snr/Pacbio_tests/out";
$dataset_name = "dataset_name";

my $bsub_script = "$FindBin::Bin/scaffold-wrapper_pbjelly/bsub_wlog.pl";
print "BSUB_SCRIPT: ".$bsub_script."\n";
# my $bsub_script =  File::Spec->rel2abs('./scaffold-wrapper_pbjelly/bsub_wlog.pl');

my $LSF_SUBLINE = $bsub_script.' ${JOB_ID} -q normal -R \'select[mem>48000] rusage[mem=48000] span[ptile=12]\' -n12  -M48000 ${CMD}';
my $LOCAL_SUBLINE='${CMD} 2> ${STDERR} 1> ${STDOUT} &amp';

my $BLASR_LINE = "";
foreach my $key (keys(%blasr)) {
    $BLASR_LINE .= " -$key $blasr{$key}"
}

my $template = <<'END';
<jellyProtocol>
    <reference>$contigs</reference>  
    <outputDir>$outdir/$dataset_name</outputDir>
    <cluster>
       <command>$COMMAND_LINE</command>
        <nJobs>$NJOBS</nJobs>
    </cluster>
    <blasr>$BLASR_LINE</blasr>
    <input baseDir="$READSDIR">
      $READS
    </input>
</jellyProtocol>

END
;

####################
# print large and small templates
####################

my $COMMAND_LINE = "";

$COMMAND_LINE = $LSF_SUBLINE;
my $tempfile_lg = File::Spec->rel2abs("./template_large.xml");
my $template_lg = $template;
$template_lg =~ s/(\$\w+)/$1/eeg;
#print $template_lg;
open(TEMPLG,">./template_large.xml");
print TEMPLG $template_lg;
close(TEMPLG);

$COMMAND_LINE = $LOCAL_SUBLINE;
my $tempfile_sm = File::Spec->rel2abs("./template_small.xml");
my $template_sm = $template;
$template_sm =~ s/(\$\w+)/$1/eeg;
#print $template_sm;
open(TEMPSM,">$tempfile_sm");
print TEMPSM $template_sm;
close(TEMPSM);
$COMMAND_LINE = "";


my $LASTID = -1;
my $LASTDEPENDENCY;
#def submit single-job job

sub _submitLocalJob {
    my $stage = shift;
    my $queue = shift;
    my $dependency = shift;
    
    $dependency = "-w \"".$dependency."\"" if $dependency;
    #print $dependency."\n" if $dependency;
    $dependency = "" if !$dependency;
    $preexec = "-E ".$preexec if $preexec;

  # bsub, write to log file...
    my $job_name = "pbj_".$dataset_name."_".$stage;
    #print $job_name."\n";
    my $command = <<"END";
 bsub -J $job_name -o $job_name.%J.out -e $job_name.%J.err \\
 -q $queue -R \'select[mem>2000] rusage[mem=2000] span[ptile=4]\' \\
 $dependency $preexec  -M2000 -n4  \\
 Jelly.py $stage $tempfile_sm
END
$command .= " -x $extras" if $extras;

#print STDERR $command;
_bsub($command)

  # sub it, grab ID,
#  my $LASTID = $1;
  # make next job dependent on completion if ID
}


#def submit multi-job job
sub _submitDistJob {
    my $stage = shift;
    my $queue = shift;
    my $dependency = shift;
    
    $dependency = "-w \"".$dependency."\"" if $dependency;
    #print $dependency."\n" if $dependency;
    $dependency = "" if !$dependency;
    $preexec = "-E ".$preexec if $preexec;

  # bsub, write to log file...
    my $job_name = "pbj_".$dataset_name."_".$stage;
    #print $job_name."\n";
    my $command = <<"END";
 bsub -J $job_name -o $job_name.%J.out -e $job_name.%J.err \\
 -q $queue -R \'select[mem>2000] rusage[mem=2000] span[ptile=4]\' \\
 $dependency $preexec  -M2000 -n4  \\
 Jelly.py $stage $tempfile_lg
END

$command .= " -x $extras" if $extras;

#print STDERR $command;
    _bsub($command);
#wait for return, 
    my $result = _wait_done();
  # make next job dependent on completion if ID
  # submit job, get Job ID,
  # wait, check, wait...
  # when complete
  # read out list of subJobs from Array
  # return array of job IDs to main runner
}

sub _wait_done {
    my $bsub = `bjobs -noheader $LASTID`;
    my @F = split($bsub);
    while $F[2] ne "EXIT" && $F[2] ne "DONE" {
	wait(60);
	my $bsub = `bjobs -noheader $LASTID`;
	my @F = split($bsub);
    }
    return $F[2];
}

sub _bsub {
    our $LASTID;
    my $command = shift;
#    my $output = `$command`;
    my $output = "thingy <123> submitted";
    if($output =~ m/\<(\d+)\>.*submitted*/gi) {
	$LASTID =  $1;
	print STDERR 'submitted: '.$command."\tjob_id: ".$LASTID."\n";
	return $output;
    }
    else{
	print STDERR "could not submit command\n".$command;
	$LASTID = "-1";}
}

# setup: make suffixarray...
_submitLocalJob("setup","normal");

# MAPPING
_submitDistJob("mapping","normal","done($LASTID)");
#check mapping worked:
#check exists 'm4' file

# SUPPORT
_submitLocal("support","normal",$LASTDEPENDENCY);
# check bml file?
exit(1);

#EXTRACT (single)
_submitLocal("extract","normal","done($LASTID)");

#ASSEMBLY (multi)
_submitDistJob("assembly","normal","done($LASTID)");

#OUTPUT (single)
_submitLocal("support","normal",$LASTDEPENDENCY);

#clean up (template file, log files, etc)
#check exists out.fasta...
# if so, delete tmp
