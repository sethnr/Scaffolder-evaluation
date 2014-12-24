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


my $local_mem = 8000;
my $dist_mem_head = 2000;
my $dist_mem_node = 48000;
my $DIST_CORES = 12;
my $preexec = "~/Pacbio_tests/setPaths.sh";

my $jelly = "Jelly.py"

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
'jelly=s'   => \$jelly,
'blasr=s'   => \$blasr_opts,
'jobs=i'    => \$NJOBS,
'cores=i'   => \$no_cores,
'debug'     => \$debug,
'distmem=i' => \$dist_mem_node,
'localmem=i' => \$local_mem,
'reads|r=s'   => \@reads,
'contigs|c=s' => \$contigs,
'o=s'       => \$outdir

);

#print join("\n",%blasr);

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
#$reads = File::Spec->rel2abs($reads);

print STDERR "SCRIPT: ".$0."\nSCRIPTDIR: ".dirname($0)."\n";

my $READSDIR = dirname($reads[0]);
$READSDIR = File::Spec->rel2abs($READSDIR);
$reads[0] = basename($reads[0]);

unless($dataset_name) {
    $dataset_name = $reads[0];
    $dataset_name =~ s/.fast[aq]$//gi;
    $dataset_name =~ s/.fa$//gi;
}
if($outdir) {
    $outdir = File::Spec->rel2abs($outdir);
}
else {
    $outdir = File::Spec->rel2abs(".");
}

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

print STDERR "DEBUG: mkdir $run_folder \n";
system("mkdir $run_folder"); # or die "could not make rundir";
chdir($run_folder) or die "could not cd into $run_folder";
my $LOGFILE = $run_folder."/logfile.txt";
open(LOGFILE,">>$LOGFILE");

print STDERR "DEBUG: mkdir $outdir/$dataset_name";
system("mkdir $outdir/$dataset_name\n"); # or die "could not make rundir";
system("ls $outdir");


#######################
#TEMPLATES
#######################
$contigs = "/lustre/scratch108/parasites/snr/refs/Pf3D7.version3.uniquely-tagged.contigs.fasta";
#$outdir = "/lustre/scratch108/parasites/snr/Pacbio_tests/out";
#$dataset_name = "dataset_name";

my $bsub_script = "$FindBin::Bin/scaffold-wrapper-pbjelly/bsub_wlog.pl";
print "BSUB_SCRIPT: ".$bsub_script."\n";
# my $bsub_script =  File::Spec->rel2abs('./scaffold-wrapper_pbjelly/bsub_wlog.pl');

my $LSF_SUBLINE = $bsub_script.' '.$LOGFILE.' -q normal -R \"select[mem\>'.$dist_mem_node.'] rusage[mem='.$dist_mem_node.'] span[ptile='.$DIST_CORES.']\" -n'.$DIST_CORES.'  -M'.$dist_mem_node.' ${CMD}';
my $LOCAL_SUBLINE='${CMD} 2> ${STDERR} 1> ${STDOUT} &amp;';

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
    $READS</input>
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
    my $preexec = "-E ".$preexec if $preexec;

  # bsub, write to log file...
    my $job_name = "pbj_".$dataset_name."_".$stage;
    #print $job_name."\n";
    my $command = <<"END";
 bsub -J $job_name -o $job_name.%J.out -e $job_name.%J.err \\
 -q $queue -R \'select[mem>$dist_mem_head] rusage[mem=$dist_mem_head] span[ptile=4]\' \\
 $dependency $preexec  -M$dist_mem_head -n4  \\
 $jelly $stage $tempfile_sm
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
    my $preexec = "-E ".$preexec if $preexec;

  # bsub, write to log file...
    my $job_name = "pbj_".$dataset_name."_".$stage;
    #print $job_name."\n";
    my $command = <<"END";
 bsub -J $job_name -o $job_name.%J.out -e $job_name.%J.err \\
 -q $queue -R \'select[mem>2000] rusage[mem=2000] span[ptile=4]\' \\
 $dependency $preexec  -M2000 -n4  \\
 $jelly $stage $tempfile_lg
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
    print STDERR "BSUB: :".$bsub.":\n";
    my $status = "";
    if($bsub) {
	my @F = split($bsub);
	$status = $F[2];
    }
    my $waits = 0;
    while ($status ne "EXIT" && $status ne "DONE") {
	print STDERR "waiting... ".$status." ".$waits."\n";
	$waits++;
	sleep(30);
	$bsub = `bjobs -noheader $LASTID`;
	print "BSUB: ".$bsub;
	if($bsub) {
	    my @F = split(/\s+/,$bsub);
	    $status = $F[2];
	}
	#if still waiting and nothing on job list - assume failed to submit
	if ($waits > 2 && $status ne "PEND" 
	    && $status ne "RUN" && $status ne "DONE") {
	    print STDERR "job $LASTID not returning status [$status] - not submitted?\n";
	    exit(1);
	}
	# check through log file for failed jobs;
	# will fail if finds one with EXIT status
	_get_all_jobs(); 
	
	
    }
    return $status;
}

sub _bsub {
    #our $LASTID;
    my $command = shift;
    #print STDERR $command."\n";
    my $output = `$command`;
#    my $output = "thingy <123> submitted";
    $command =~ s/\n/ /gi;
    if($output =~ m/\<(\d+)\>.*submitted*/gi) {
	$LASTID = $1;
	print LOGFILE $LASTID."\t".$LASTID."\t".$command."\n";
	return $output;
    }
    else{
	print LOGFILE "-1\t-1\tcould not submit command\n".$command;
	$LASTID = "-1";}
}

sub _get_child_jobs {
    my $parent_id = shift;
    my @children;
    open(READLOG,$LOGFILE);
    while(<READLOG>) {
	print $_;
	my @F = split("\t",$_);
	my ($parent, $child) = @F[0..1];
	if ($parent == -1) {print STDERR "job not submitted, exiting\n";
			    exit(1);
	}
	if ($parent == $parent_id) {
	    push(@children, $child);
	}
    }
    close(READLOG);
    return \@children;
}
sub _get_all_jobs {
    my @jobs;
    open(READLOG,$LOGFILE);
    while(<READLOG>) {
	my @F = split("\t",$_);
	my ($parent, $job) = @F[0..1];
	if ($job == -1) {print STDERR "job not submitted, exiting\n";
			 exit(1);}
	my $status = _get_job_status($job);
	if ($status eq "EXIT") {
	    print STDERR "job ".$job." exited with status ".$status.", exiting\n";
	    exit(1);}
	
	push(@jobs, $job);
    }
    close(READLOG);
    return \@jobs;
}



sub _get_job_status {
    my $job_id = shift;
    my $bjob = `bjobs -noheader $job_id`;
    my $status = "NULL";
    if($bjob && $bjob ne "") {
	my @F = split(/\s/,$bjob);
	$status = $F[2];
    }
    return $status;
}

sub _make_depends_children {
    my @depends = ();
    foreach my $child (@_) {
	my $depend = "done(".$child.")";
	push(@depends,$depend);
    }
    return join(" && ",@depends);
}

# setup: make suffixarray...
_submitLocalJob("setup","normal");

# MAPPING
_submitDistJob("mapping","normal","done($LASTID)");
#check mapping worked:

my $children = _get_child_jobs($LASTID);
foreach my $child (@{$children}) {
    print $child."\t";
    my $status =  _get_job_status($child);
    print $status."\n";
}
$LASTDEPENDENCY = _make_depends_children(@{$children});
print STDERR $LASTDEPENDENCY."\n";
#check exists 'm4' file?

# SUPPORT
_submitLocalJob("support","normal",$LASTDEPENDENCY);
# check bml file?

#EXTRACT (single)
_submitLocalJob("extract","normal","done($LASTID)");

#ASSEMBLY (multi)
_submitDistJob("assembly","normal","done($LASTID)");

my $children = _get_child_jobs($LASTID);
$LASTDEPENDENCY = _make_depends_children(@{$children});

#OUTPUT (single)
_submitLocalJob("support","normal",$LASTDEPENDENCY);

#clean up (template file, log files, etc)
#check exists out.fasta...
# if so, delete tmp
