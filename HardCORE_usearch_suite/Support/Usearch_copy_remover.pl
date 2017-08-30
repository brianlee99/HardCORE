#!/usr/bin/perl

use strict;
use warnings;
use Parallel::ForkManager;
use Sys::Info;
use Sys::Info::Constants qw(:device_cpu);


my $threads = $ARGV[0];
my $info;
my $cpu;
my $cpus;
my $pm;

if($threads == 0){
    $info=Sys::Info->new;
    $cpu  = $info->device( CPU => my %options );
    $cpus=scalar($cpu->identify);
    $pm=Parallel::ForkManager->new($cpus);
}

else{
    $pm=Parallel::ForkManager->new($threads);
}

my @files=glob "*.faa";

RUN:

foreach my $file (@files)
{
$pm->start and next RUN;
my $name=$file;
$name =~ s{\.[^.]+$}{};
my $name_sort=$name."_sort.faa";

USEARCH_SORT($file,$name_sort);
USEARCH_CLUSTER($name_sort,$name);
REMOVE_COPIES($file,$name);

sub REMOVE_COPIES{

my $file=shift @_;
my $name=shift @_;
my (%fasta,$new_fasta,%remove);
open(CLUSTER,"<usearch_cluster_$name.txt");
while (<CLUSTER>)
	{
		chomp $_;
		my @line=split("\t",$_);
		$remove{$line[0]}=();
	}
close CLUSTER;
open(FASTA,"<$file");
while (<FASTA>)
	{
		if ($_=~/$name/)
			{
				$new_fasta=$new_fasta.$_
			}
		else
			{
				$_=~s/\n//g;
				$new_fasta=$new_fasta.$_;
			}
	}
close FASTA;
my @new_fasta=split('>',$new_fasta);
shift @new_fasta;

foreach (@new_fasta)
	{
		my @fasta=split("\n",$_);
		$fasta[0]=~s/\s.*//;
		$fasta{$fasta[0]}=$fasta[1];
	}

my $name_refine=$name.".refine.faa";
open(OUT,">$name.refine.faa");
 
foreach my $key (keys %fasta)
	{
		unless (exists $remove{$key})
			{
				print OUT ">$key\n$fasta{$key}\n"
			}
	}
close OUT;
}#end REMOVE_COPIES


sub USEARCH_CLUSTER{

my $name_sort=shift @_;
my $name=shift @_;
system("usearch -cluster_smallmem $name_sort -id 0.7 -trunclabels --blast6out usearch_cluster_$name.txt");
}#end USEARCH_CLUSTER


sub USEARCH_SORT{
my $file=shift @_;
my $name_sort=shift @_;
system("usearch -sortbylength $file -fastaout $name_sort");
return $name_sort;
}#end USEARCH_SORT

$pm->finish;
}
$pm->wait_all_children;

system("mkdir Clusters");
system ("rm *sort*");
system ("mv usearch_cluster*txt ./Clusters");

exit;
