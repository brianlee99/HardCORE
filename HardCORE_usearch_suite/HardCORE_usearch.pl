#! /usr/bin/perl

use Cwd;
use List::Util qw(first);
use Getopt::Long;
use Term::ANSIColor;
use FindBin;

sub ProkkaPrep_faa
{
	my $home = $_[0];
	chdir($home);
	$home = getcwd; #if the user enters './' as a dir this will solve issues
	chdir($home);
	my @prokka_files = glob "*.faa";
	system("mkdir PROKKA_backups");
	foreach(@prokka_files){
		my $Genome = $_;
		$Genome =~ s{\.[^.]+$}{}; #strip .ffn
		open(CHECK,"<$Genome.faa") or die print "\n\nERROR\n\nCannot open $Genome.faa";
		my $firstLine = <CHECK>;
		my @split = split(/\|/,$firstLine);
		if($split[0] ne ">gi" && $split[2] ne "ref" && $split2[1] ne "\[$Genome\]")
		{
			open(IN, "<$_") or die;
			open(OUT,">$Genome.mod.faa") or die;

			while(<IN>){
				chomp();
				my $line = $_;
				if($_ =~ m/>/){
					my @split = split(/\_/,$line);
					my @split2 = split(" ",$split[1]);
					my $prokka_id = $Genome . "_" . $split2[0];
					my $new = ">gi|$prokka_id|ref|1234| $prokka_id "."[$Genome]";
					print OUT $new . "\n";
				}#end if
				else{print OUT $line . "\n";}
			}#end while IN
			close(IN);
			close(OUT);
			close(CHECK);
			system("mv $Genome.faa ./PROKKA_backups/$Genome.orig.faa");
			system("mv $Genome.mod.faa $Genome.faa");
		}#end if
	}#end foreach
}#end ProkkaPrep_faa

sub AllvsAll_uSearch
{
	my $home = $_[0];
	chdir($home);
	$home = getcwd; #if the user enters './' as a dir this will solve issues
	chdir($home);
	my $id = $_[1];
	my $plen = $_[2];
	my $threads = $_[3];
	$Program_dir = "$FindBin::RealBin/Support";
	#Get Protein List
	#system("vsearch --allpairs_global AllProteinList.all --self --userout AllProteinList.vsearch --id 0.80 --minqt 0.9 --strand both --userfields query+target+qs+tl+alnlen+ids+qstrand+tstrand");
	#system("usearch -sortbylength AllProteinList.all -fastaout AllProteinList.all.uready");
	#system("usearch -cluster_smallmem AllProteinList.all.uready -id 0.70 -userout AllProteinList.all.cluster -userfields query+target+qs+tl+alnlen+ids");
	#system("usearch -allpairs_global AllProteinList.all.uready -id 0.70 -userout AllProteinList.all.cluster -userfields query+target+qs+tl+alnlen+ids");
	system("perl $Program_dir/Usearch_copy_remover.pl $threads");
	system("perl $Program_dir/Parallel_usearch.pl $id $threads");
	system("rm *.comb");
	system("touch AllProteinList.all.cluster");
	system("cat *.usearch >> AllProteinList.all.cluster");
	#system("rm *.usearch");
	open(IN, "<AllProteinList.all.cluster") or die print "Could not open AllProteinList.all.cluster\n";
	open(OUT,">AllProteinList.all.cluster.fixed");
	#Need to replace the second column with proper tags (i.e. gi|SI01_00001|ref|1234| SI01_00001 [SI01])
	#Also will interchange the target and query columns!
	while(<IN>){
		my $line = $_;
		my @split = split(/\t/,$line);
		my @tag = split(/\|/,$split[0]);
		my @strain = split(/\_/,$tag[1]);
		my $query = "$split[0] $tag[1] \[$strain[0]\]";
		my $target = $split[1];
		my $qlen = $split[3];
		my $tlen = $split[2];
		my $aln = $split[4];
		my $ident = $split[5];
		my $lenfilt;
		if($qlen >= $tlen){$lenfilt = $tlen/$qlen;}
		if($tlen > $qlen){$lenfilt = $qlen/$tlen;}
		################################################################################
		#                     Threshold Adjustment For Length Score!!                  #
		if($lenfilt >= $plen ){print OUT "$target\t$query\t$qlen\t$tlen\t$aln\t$ident";}#
          ################################################################################
	}#end while(<IN>)
	close(IN);
	close(OUT);
	system("rm AllProteinList.all.cluster");
	system("mv AllProteinList.all.cluster.fixed AllProteinList.all.cluster");
}#end AllvsAll_VSEARCH

sub SortUSearch
{
	my $home = $_[0];
	chdir($home);
	$home = getcwd; #if the user enters './' as a dir this will solve issues
	chdir($home);
	open(US,"<AllProteinList.all.cluster") or die print "Could not open AllProteinList.all.cluster\n";
	open(OUT,">AllProteinList.all.cluster.sorted");
	
	#Going to take the first two columns ($Query and $Target) and use them as key.
	#The value is going to be the entire line itslef.
	my %SORT;
	while(<US>){
		my $line = $_;
		my @split = split(/\t/,$line);
		my $key = $split[0] . $split[1];
		my $value = $line;
		$SORT{$key}=$value;
	}#end while
	foreach(sort keys %SORT){
		print OUT $SORT{$_};
	}#end foreach
	close(US);
	close(OUT);
}


sub ParseBlast_USearch
{
	my $ffn_dir = $_[2];
	my $home = $_[0];
	my $threads = $_[3];
	chdir($home);
	$home = getcwd; #if the user enters './' as a dir this will solve issues
	chdir($home);
	open (ALL,"<AllProteinList.all") or die;
	my %PAN; #This will hold the entire PAN GENOME
	my %AllProteins; #Make a hash where the key is the GI SeqID and the value is the "NOT_YET"
	my %FULL_QUERY; #need a seperate hash for keys to hold the info for metadata
	my $key; #this will store the genome title
	my $value;
	my $count=0; #need a count to make sure I don't push nothing onto the hash
	while(<ALL>){
		my $line = $_;
		if($line =~ m/>/){
			my @split = split(/\|/,$line);
			my $gi = $split[1];
			$AllProteins{$gi}="NOT_YET";
		}#end count if
	}#end While ALL
	close(ALL);
	#Creates a hash of all the strain names
	open(STRAINS,"<$_[1]") or die;
	my %STRAINS;
	my $Hash_strain_size = keys %STRAINS;
	while(<STRAINS>){
		chomp();
		$STRAINS{$_}="NOT_DONE";	
	}
	close(STRAINS);
	
	#Start working on the PAN
	system("mkdir Pan_Genome");
	open (ALL,"<AllProteinList.all.cluster.sorted") or die print "Cannot open AllProteinList.all.cluster.sorted\n";
	my $previous;
	my $burned_GI;
	open(PAN,">$home/Pan_Genome/PanGenome.usearch");
	while(<ALL>){
		chomp();
		my $line = $_;
		my @Tab_split = split(/\t/,$line);
		my @Query_split = split(/\|/,$Tab_split[0]);
		my $Query_GI = $Query_split[1]; #This will yield QueryGenome_#####
		my @Query_split = split(/\_/,$Query_GI);
		my $Query_strain = "\[$Query_split[0]\]";
		my @Hit_split= split(/\|/,$Tab_split[1]);
		my $Hit_GI = $Hit_split[1]; #This will yield HitGenome_#####
		my $Hit_title = $Hit_split[4];
		$Hit_title=~ s/\s//g;
		my @Hit_strain = split(/\[/,$Hit_title);
		my $Hit_strain = $Hit_strain[1]; #This will yield HitGenome
		#my $Qstrand = $Tab_split[6];
		#my $Tstrand = $Tab_split[7];
		chop($Hit_strain);
		my $Query_out = "\>$Query_GI\$$Query_GI$Query_strain";
		my $Hit_out = ">$Hit_GI\$$Hit_title"; # >Genome_#####$Genome
		#print "Query_GI:$Query_GI\;Hit_GI:$Hit_GI\;Hit_title:$Hit_title\;Hit_out:$Hit_out";
		#This first 'if' will only create a new file for the hits if its the first instance of the GI
		if($AllProteins{$Query_GI} ne "KEY" && $Hit_GI ne $Query){ 
			$PAN{$Query_GI} .= "$Hit_out\t"; #concatenating multiple records onto one key
			$AllProteins{$Hit_GI} = "KEY";
			#$FULL_QUERY{$Query_GI}="$Query_out\t";
		}#end if
		if($Hit_GI eq $Query){
			$PAN{$Query_GI} .= "$Hit_out\t"; #concatenating multiple records onto one key
			$AllProteins{$Hit_GI} = "KEY";
		}#end if
	}#end while

	#open(OUT, ">HashTest") or die;
	#while (my ($key, $value) = each %AllProteins) {
	#			print OUT "$key $value\n";
	#}#end hash print While
	#close(OUT);

	foreach(sort keys %PAN)
	{
		print PAN $PAN{$_} . "\n";
	}
	close(PAN);

	open(PAN,"<$home/Pan_Genome/PanGenome.usearch");
	my %PAN;
	while(<PAN>){
		chomp();
		my $line = $_;
		my @split = split(/\t/,$line);
		my $length = length $PAN{$split[0]};
		my $new;
		for(my $x=0; $x<= scalar @split; $x++){
			if($PAN{$split[0]} && $length < scalar @split && $split[0] ne $split[$x]){				
				$new .= "$split[$x]\t";
			}
			elsif($split[0] ne $split[$x]){
				$new .= "$split[$x]\t";
			}
		}#end foreach
		$PAN{$split[0]} = $new;
		$new="";
	}#end while
	close(PAN);

	open(OUT, ">$home/Pan_Genome/PanGenome.usearch2") or die;
	while (my ($key, $value) = each %PAN) {
				print OUT "$key\t$value\n";
	}#end hash print While
	close(OUT);

	#Need to get rid of extra \t and get rid of the doubling up of seq names and numbers between the $

	open(PAN, "<$home/Pan_Genome/PanGenome.usearch2") or die;
	open(OUT, ">$home/Pan_Genome/PanGenome.usearch3") or die;
	while(<PAN>){
		my $line = $_;
		my @split = split(/\t/,$line);
		foreach(@split){
			my @split2 = split(/\$/,$_);
			my $length = length($split2[1]);
			if($length > 5){print OUT ">$split2[1]\t";} #this is to ensure I don't print extra '>'
		}#end foreach
		print OUT "\n";
	}#end while

	close(OUT);
	close(PAN);
	
	system("rm $home/Pan_Genome/PanGenome.usearch");
	system("rm $home/Pan_Genome/PanGenome.usearch2");
	system("mv $home/Pan_Genome/PanGenome.usearch3 $home/Pan_Genome/PanGenome.usearch");

	chdir("$home/Pan_Genome");
	$Program_dir = "$FindBin::RealBin/Support";
	system("perl $Program_dir/reFFNer_v2.pl ./PanGenome.usearch $ffn_dir $threads");
	chdir($home);

}#end sub ParseBlast


sub Get_Core_Singletons
{
	my $home = $_[0];
	chdir($home);
	$home = getcwd; #if the user enters './' as a dir this will solve issues
	chdir($home);
	#Need to find out how many STRAINS there are:
	open(STRAINS,"<$_[1]") or die;
	$strain_count = 0;
	while(<STRAINS>){
		$strain_count++;
	}#end while
	#print $strain_count . "\n";
	close(STRAINS);
	system("mkdir Core_Genome");
	system("mkdir Singletons");
	#open the Pan Genome file
	open(PAN,"<$home/Pan_Genome/PanGenome.usearch") or die;
	open(CORE,">$home/Core_Genome/CoreGenome.usearch");
	open(SINGLE,">$home/Singletons/Singletons.usearch") or die;

	while(<PAN>){
		chomp();	#Need to chomp() twice to get rid of the \n 
		chomp();	#and the extra \t
		my $line = $_;
		my @split_tab = split(/\t/,$line);
		#print scalar @split_tab . "\n";		
		if(scalar @split_tab == 1){print SINGLE "$line\n";}
		if(scalar @split_tab == $strain_count){print CORE "$line\n";}
	}#end while(<PAN>)

	close(CORE);
	close(PAN);
	close(SINGLE);
}#end Get_Core

#The user will input if they want the fastas of the Pan, Core or Singletons
sub GetFaa{
	my $home = $_[0];
	chdir($home);
	$home = getcwd; #if the user enters './' as a dir this will solve issues
	chdir($home);
	my $choices_ref = shift;
	my $strains = shift;
	my @choices = @{$choices_ref};
	my $strain_count; #This will help for getting the core genome
	open(STRAINS,"<$strains") or die;
	$strain_count = 0;
	while(<STRAINS>){
		$strain_count++;
	}#end while
	#Going to store AllProteins.all in a hash. WIll make sure the key matches the 
	#output from the PanGenome.usearch file
	open(AllP,"<AllProteinList.all") or die "\n\tCannot open AllProteinList.all\n";
	my %AllProteins;
	my $title_key;
	my $seq_value;
	my $first = 0;
	while(<AllP>){
		chomp();
		my $line = $_;
		if($line =~ m/\>/){
			if($first != 0){
				$AllProteins{$title_key} = $seq_value;
				$seq_value="";
				$title_key="";
			}
			my @split = split(" ",$line);
			$title_key = ">$split[1]$split[2]";
			$first=1;
		}#end if
		else{
			$seq_value .= $line . "\n";
		}
		#Need to hit the last one
		$AllProteins{$title_key} = $seq_value;
	}#end while(AllP)
	foreach(@choices){
		#print $_ . "\n";
		if($_ eq "Pan"){
			chdir("Pan_Genome")or die "\n\tCan't find the Pan_Genome directory.\n\n";
			system("mkdir Pan_Fastas");
			open(PAN,"<PanGenome.usearch")or die "\n\tCan't open PanGenome.usearch\n\n";
			while(<PAN>){
				chomp();
				my $line = $_;
				my @tab_split = split(/\t/,$line);
				my @split2 = split(/\[/,$tab_split[0]);
				my $out_file = $split2[0];
				$out_file =~ s/\>//;
				open(OUT,">>./Pan_Fastas/$out_file.fasta");
					foreach(@tab_split){
						print OUT "$_\n$AllProteins{$_}";#Use the hash call to print to a file
					}#end foreach
				close(OUT);
			}#end while	
			close(PAN);
			chdir($home);
		}#end if(PAN)

		elsif($_ eq "Core"){
			chdir("Core_Genome")or die "\n\tCan't find the Core_Genome directory.\n\n";
			system("mkdir Core_Fastas");
			open(CORE,"<CoreGenome.usearch")or die "\n\tCan't open CoreGenome.usearch\n\n";
			while(<CORE>){
				chomp();
				my $line = $_;
				my @tab_split = split(/\t/,$line);
				if(scalar @tab_split == $strain_count){
					my @split2 = split(/\[/,$tab_split[0]);
					my $out_file = $split2[0];
					$out_file =~ s/\>//;
					open(OUT,">>./Core_Fastas/$out_file.fasta");
						foreach(@tab_split){
							print OUT "$_\n$AllProteins{$_}";#Use the hash call to print to a file
						}#end foreach
					close(OUT);
				}#end if
			}#end while	
			close(CORE);
			chdir($home);
		}#end elsif(CORE)

		elsif($_ eq "Singletons"){
			chdir("Singletons")or die "\n\tCan't find the Pan_Genome directory.\n\n";
			system("mkdir Singleton_Fastas");
			open(SING,"<Singletons.usearch")or die "\n\tCan't open Singletons.usearch\n\n";
			while(<SING>){
				chomp();
				my $line = $_;
				my @tab_split = split(/\t/,$line);
				if(scalar @tab_split == 1){
					my @split2 = split(/\[/,$tab_split[0]);
					my $out_file = $split2[0];
					$out_file =~ s/\>//;
					open(OUT,">>./Singleton_Fastas/$out_file.fasta");
						foreach(@tab_split){
							print OUT "$_\n$AllProteins{$_}";#Use the hash call to print to a file
						}#end foreach
					close(OUT);
				}#end if
			}#end while	
			close(SING);
			chdir($home);
		}#end elsif(Singletons)
	}#end foreach
}#sub GetFaa


sub Core_Subset{
	my $home = $_[0];
	chdir($home);
	$home = getcwd; #if the user enters './' as a dir this will solve issues
	chdir($home);
	my $file = $_[1];
	my @Subset;
	open(OUT,">$file.core.subset");
	open(SUBSET,"<$file");
	#Load the subset strains into an array
	while(<SUBSET>){
		chomp();
		push(@Subset,$_);
	}#end while(<SUBSET>)
	#foreach(@Subset){print $_ . "\n";}
	close(SUBSET);
	#Will be comparing arrays, always easier if sorted first
	my @sorted_subset = sort @Subset;
	my $subset_size = scalar @Subset;
	#Now need to open pan genome file and search each line for the subset
	open(PAN,"$home/Pan_Genome/PanGenome.usearch") or die "\n\n\tCould not open PanGenome.usearch file!\n";
	while(<PAN>){
		chomp();
		my $line = $_;
		my @tab_split = split(/\t/,$line);
		my $line_size = scalar @tab_split;
		my $count = 0;
		if($subset_size == $line_size){
			my @strains;
			foreach(@tab_split){
				my @strain_split= split(/\[/,$_);
				my $strain = $strain_split[1];
				$strain =~ s/\]//; #I now have the strain that is the []
				#print "$strain\n";
				push(@strains,$strain);
			}#end foreach
			my @sorted_strains = sort @strains;
			$count = 0;
			for(my $i=0; $i< $line_size; $i++){
				if($sorted_strains[$i] eq $sorted_subset[$i]){$count++;}
				else{last;}
			}#end for
			#Here is where we have a subset found and now need to export the necessary info
			if($count == $subset_size){
				print OUT "$line\n";
			}#if($count == my $subset_size)
		}#end if($subset_size == $line_size)
	}#end while(<PAN>)
	close(PAN);
	close(OUT);
}#end sub Core_Subset


########################################################################################################################
#			MAKE SURE THERE ARE NO UNDERSCORES IN STRAIN NAMES
########################################################################################################################
########################################################################################################################
sub MAIN{

	#Read options into a hash
	my %options;
	GetOptions (\%options, "dir=s", "strains=s", "fasta=s", "core_subset=s", "ffn=s", "id=s", "plen=s", "threads=s", "only_get_fasta","help", "h");

	if($options{'help'} || $options{'h'}){
		print color ('bold');
		print colored(['bright_red on_bright_yellow'], "\n\t", 'USAGE:');
		print color ('reset');
		print colored(['bright_red'], "\n\t", 'Please note: For the strain file, this is just a list of the strains.', "\n\t", 'Each strain is followed by a line break and no ".faa" file designation.'); 
		print color ('bold');	
		print colored(['bright_red on_black'], "\n\n\t", 'Most basic run:');	
		print color ('reset');
		print "\n\tHardCORE_usearch.pl ";
		print colored(['red on_bright_yellow'], '-dir');  
		print " /path/to/*.faa(s) ";
		print colored(['red on_bright_yellow'], '-strain');
		print " user_strain_list_file\n\n";
		exit;
	}#end if
	if($options{'dir'} && $options{'strains'} && $options{'ffn'} && $options{'id'} && $options{'plen'} && !$options{'only_get_fasta'} && !$options{'core_subset'}){
        #Specifies the number of threads to be used throughout the analysis
        my $threads;
   		if($options{'threads'}){$threads = $options{'threads'};}
   		else{$threads=0;}
   		
		my $dir = $options{'dir'};
		my $strains = $options{'strains'};
		my $ffn_dir = $options{'ffn'};
		my $id = $options{'id'};
		my $plen = $options{'plen'};
		chdir($dir)or die "\n\tCan't find $dir\n\tIs it a valid directory with all your .faa files and strains file?\n\n";
		ProkkaPrep_faa($dir);
		system("cat *.faa > AllProteinList.all"); #Making a global list of proteins
		print colored(['red on_bright_yellow'], "\n", 'At Usearch stage', "\n");		
		AllvsAll_uSearch($dir, $id, $plen, $threads);
		SortUSearch($dir);
		print color ('bold red');
		print "\n" . "At parse stage" . "\n";	
		print color ('reset');	
		ParseBlast_USearch($dir,$strains,$ffn_dir,$threads);
		print color ('bold red');
		print "\n" . "Generating Pan - Core - Singleton genome files" . "\n";
		print color ('reset');	
		Get_Core_Singletons($dir,$strains);
		print color ('bold red');
		chdir($dir);
		system("rm *.refine*");
		print "\n" . "All done basic run!" . "\n\n";
		print color ('reset');	
	}#end if
	if($options{'fasta'} || $options{'only_get_fasta'}){
			my $dir = $options{'dir'};
			chdir($dir)or die "\n\tCan't find $dir\n\tIs it a valid directory with all your .faa files?\n\n";
		if(lc $options{'fasta'} eq lc "Core"){
			print color ('bold red');
			print "\n\t" . "Gathering the fastas of all the CORE sequences.\n\n";
			print color ('reset');
			my @choice = "Core";
			GetFaa(\@choice,$options{'strains'});
		}#end if(Core)
		elsif(lc $options{'fasta'} eq lc "Pan"){		
			print color ('bold red');
			print "\n\t" . "Gathering the fastas of all the PAN sequences.\n\n";
			print color ('reset');
			my @choice = "Pan";
			GetFaa(\@choice,$options{'strains'});
		}#end if(Pan)
		elsif(lc $options{'fasta'} eq lc "Singletons"){		
			print color ('bold red');
			print "\n\t" . "Gathering the fastas of all the SINGLETON sequences.\n\n";
			print color ('reset');
			my @choice = "Singletons";
			GetFaa(\@choice,$options{'strains'});
		}#end if(Singletons)
		elsif(lc $options{'fasta'} eq lc "All"){
			print color ('bold red');
			print "\n\t" . "Gathering the fastas of all the CORE PAN and SINGLETON sequences.\n\n";
			print color ('reset');
			my @choice = ("Core", "Pan", "Singletons");
			GetFaa(\@choice,$options{'strains'});
		}#end elsif(ALL)
		else{
			print color ('bold red');
			print "\n\t" . "The ";
			print colored(['red on_bright_yellow'], '-fasta');  
			print color ('bold red');
			print " option is required and must designate either Pan Core Singletons or All\n\n";
			print color ('reset');	
			exit;
		}#end else
	}#end if($options{'fasta'})

	if($options{'core_subset'} && $options{'dir'} && $options{'strains'}){
		my $dir = $options{'dir'};
		open(SUBSET,"<$options{'core_subset'}") or die "\n\n\tCannot open core_subset file. Must provide a file similar to STRAINS file that contains the desired subset.\n\n";
		close(SUBSET);
		Core_Subset($dir,$options{'core_subset'});
	}#end if($options{'core_subset'})



}#end MAIN

MAIN();

exit;


