#! /usr/bin/perl

	
use Algorithm::Combinatorics qw( combinations );
use Parallel::ForkManager;
use Sys::Info;
use Sys::Info::Constants qw(:device_cpu);

my $id = $ARGV[0];
my $threads = $ARGV[1];
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
    
    
    
    

my @strings = glob "*.refine.faa";
my $iter = combinations(\@strings, 2);

RUN:

while (my $c = $iter->next) {
	$pm->start and next RUN;
	my @files = @$c;	
	$a= $files[0]; #each of the files in this combination
	$b= $files[1];
	$a =~ s/.faa//; #remove the file designation
	$b =~ s/.faa//;

	Usearch($a,$b);

	$pm->finish;
}
$pm->wait_all_children;
exit;

sub Usearch{
	$a= $_[0]; #each of the files in this combination
	$b= $_[1];
	system("cp $b.faa ./$a$b.comb");
	system("cat $a.faa >> $a$b.comb"); #concatenate the two files together
	#Threshold Adjustment!!
	system("usearch -usearch_global $a$b.comb -threads $threads -db $a$b.comb -id $id -maxaccepts 2 -userout $a$b.usearch -userfields query+target+qs+tl+alnlen+ids");
}#end Usearch

