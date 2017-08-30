#! /usr/bin/perl
use Algorithm::Combinatorics qw( combinations );
use Parallel::ForkManager;
use Sys::Info;
use Sys::Info::Constants qw(:device_cpu);
use Cwd;
#reFFNer
#This will go through the PanGenome.usearch file and double check every single entry

#First read the entire PanGenome.usearch into a hash
sub PanGenomeHash{
	my $path = $_[0];
	open(PAN,"<$path") or die "\n\nCould not open PanGenome.usearch\n\n";
	my %PAN;	
	while(<PAN>){
		chomp();
		my $line = $_;
		my @split = split(/\t/,$line);
		my $value;
		#Want to get rig of the [xxxx] so that the values match the keys
		for(my $i=0; $i < scalar @split; $i++){
			my @value_split = split(/\[/,$split[$i]);
			$value .= $value_split[0] . "\t";
		}#end for()
		my @split2 = split(/\[/,$_); #removing the [XXXX] from each key
		$PAN{$split2[0]} = $value;
	}#end while()	
	return \%PAN; #This will return a reference to the HASH
	close(PAN);
}

sub FFN_Hash{	#Need to concatenate all FFNs and read all entries into a hash
	my $path = $_[0];
	my %FFN;
	system("cat $path/*.ffn >> ALL.ffn");
	open(FFN,"<ALL.ffn") or die "\n\nCould not open ALL.ffn \n\n";
	my $first = 0;
	my $key="";
	my $value="";
	while(<FFN>){
		my $line = $_;
		if($line =~ m/\>/){
			if($first != 0){
				$FFN{$key}="$key\n$value";
				$key="";				
				$value="";
			}
			my @split = split(" ",$line);
			$key = $split[0];
			$first =1; 	
		}#end if
		else{
			$value .= $line;
		}#end else
	}
	$FFN{$key}="$key\n$value";
	close(FFN);
	return \%FFN;
}

sub Master_Hash{	#Need to concatenate all FFNs and read all entries into a hash
	my %MASTER;
	open(FFN,"<ALL.ffn") or die "\n\nCould not open ALL.ffn \n\n";
	my $first = 0;
	my $key="";
	my $value="";
	while(<FFN>){
		my $line = $_;
		if($line =~ m/\>/){
			if($first != 0){
				$FFN{$key}="Not_Yet";
				$key="";				
			}
			my @split = split(" ",$line);
			$key = $split[0];
			$first =1; 	
		}#end if
	}
	$FFN{$key}="Not_Yet";
	close(FFN);
	return \%FFN;
}

sub Checker{
	#Need a New Temp Directory
	system("mkdir Self_Blast_Temp");
	my $home_dir = getcwd;
	#Get all the hashes you need	
	my $Path_to_PanGenome = $_[0];
	my $Path_to_ffns = $_[1];
	my $threads = $_[2];
	my $FFN_HASH = FFN_Hash($Path_to_ffns);
	my $PAN_HASH = PanGenomeHash($Path_to_PanGenome);
	#my $MASTER_HASH = Master_Hash(); Don't think I need MASTER_HASH anymore
	my @To_Delete :shared;
	chdir("$home_dir/Self_Blast_Temp");

    ############################################################################################
    # PART 1: CHECK IF WE ARE USING THE BEST REPRESENTITIVE SEQUENCE FOR THE PAN_GENOME FAMILY
    ############################################################################################

	foreach my $key (keys %{ $PAN_HASH }) {
		my $file_name = $key; #Don't want the > in the file name
		substr($file_name, 0, 1) = "";
		my @Value_Split = split(/\t/,$PAN_HASH->{$key});
		$num_values = scalar @Value_Split;
		#This is looking at the Pan_Genome one line at a time
		#If there are more than 3 genomes on one line (i.e. more than 3 genomes in a Pan_genome family)
		#A new file is created ######.SELF_BLAST.
		#This file contains each gene from each genome that is part of that Pan_Genome family.
		if($num_values >= 3){ 
			open(SELF_BLAST,">$file_name.SELF_BLAST");
			foreach(@Value_Split){
				print SELF_BLAST "$FFN_HASH->{$_}";
			}#end foreach()
			close(SELF_BLAST);
		}#end if
	}#end foreach()
	#Now we are going to go through each Pan_Genome family with 3 or more members using the ####.SELF_BLAST files generated
	#Each one is going to be BLASTed against itself
	my @files = glob "*.SELF_BLAST";
	my $info;
    my $cpu;
    my $cpus;
    my $pm;
	#Get things ready to run in parallel
	if($threads == 0){
        $info=Sys::Info->new;
        $cpu  = $info->device( CPU => my %options );
        $cpus=scalar($cpu->identify);
        $pm=Parallel::ForkManager->new($cpus);
    }

    else{
        $pm=Parallel::ForkManager->new($threads);
    }
		
	RUN:
		foreach(@files){
			$pm->start and next RUN;
			my $name = $_;
			$name =~ s{\.[^.]+$}{};
			system("formatdb -p F -i $_");
			#Here is the self-blast
			system("blastn -query $_ -db $_ -outfmt '6 qseqid stitle score' -word_size 15 > $name.BLAST");
			system("rm $_.nhr");
			system("rm $_.nin");
			system("rm $_.nsq");
			system("rm $_"); #This removes the FFN file with all the genes from the Pan_Genome family
			#The ####.BLAST is full of self-hits, need to be rid of them
			open(BLAST,"<$name.BLAST");
			open(OUT,">$name.sort.BLAST");
			my %TEMP;
		   #Here is where we get rid of self-hits and save them in #####.sort.BLAST
		   # #####.sort.BLAST is the file we will be using from here on out
		   # For --DEBUG-- purposes look at the #####.sort.BLAST files
			while(<BLAST>){
				chomp();
				my @split = split(/\t/,$_);
				if($split[0] ne $split[1]){
					my $key = "$split[2]\t$split[0]";
					my $value = "$split[0]\t$split[1]\t$split[2]";
					$TEMP{$key}=$value;
				}
			}
			foreach my $key (sort {$b <=> $a} keys %TEMP){
				print OUT $TEMP{$key} . "\n";
			}
			close(OUT);
			close(BLAST);
			system("rm $name.BLAST"); #Just need to keep the sorted blast file
			$pm->finish;
		}#end foreach()
		$pm->wait_all_children;
		    
	##BIG STEP## 
	#Now need to check first line to see if we are currently using the Ref with the best score
			my @files = glob "*.sort.BLAST";		
			foreach(@files){
				my $name = $_;
				$name =~ s/.sort.BLAST//;
				
				open my $file, '<', "$name.sort.BLAST";
				my $Actual_Ref = "";
				while(<$file>){
				    my $line = $_;
				    my @split = split(/\t/,$line);
				    #Check if current "top hit" is already a key
				    #In this case, break if the top score genome is already the ref
				    if($split[0] eq $name){
				        #print"$split[0] First IF\n" . $PAN_HASH->{"\>$name"} . "\n\n";
				        $Actual_Ref = $split[0];
				        last;
				    }
				    #If the best score genome doesn't exist as a key, that is the actual ref
				    if(!exists $PAN_HASH->{"\>$split[0]"}){
  				        #print"$split[0] Key does not exist, split[0]\n" . $PAN_HASH->{"\>$name"} . "\n\n";
                        $Actual_Ref = $split[0];
                        #Give the value of the current ref to the actual one
                        $PAN_HASH->{"\>$Actual_Ref"} = $PAN_HASH->{"\>$name"};
                        #Delete old ref as a key
                        delete $PAN_HASH->{"\>$name"};
                        last;
 			        }
			        #If top hit left column is a key, check the right column
			       if(exists $PAN_HASH->{"\>$split[0]"} && !exists $PAN_HASH->{"\>$split[1]"}){
   			           #print"$split[0] Key does exist, but key split[1] does not\n" . $PAN_HASH->{"\>$name"} . "\n\n";
			           $Actual_Ref = $split[1];
                       $PAN_HASH->{"\>$Actual_Ref"} = $PAN_HASH->{"\>$name"};
                       delete $PAN_HASH->{"\>$name"};
                       last;
			       }
			       #If neither condition is met, go to the next line to find Actual Ref
                   else {next;}
   			    }#end while
				close($file);
			}#end foreach(@files)		
			
		#system("rm formatdb.log");
		chdir($home_dir);
		system("rm -r Self_Blast_Temp"); #Can now get rid of the folder
	
	
	#Ths hash print would allow you to see the KEY for each Pan_Genome family as well as its values on the next line	
=c
		open(HASH,">PAN_HASH_TEST");
		while (my ($key) = each(%$PAN_HASH)) {
			print HASH $key . "\n" . $PAN_HASH->{$key} . "\n";
		}
		close(HASH);
=cut
	
		
		
	#########################################################
    #PART 2:
    ###########################################################
        
	chdir($home_dir);
	system("mkdir RepBlast_Temp");
	chdir("$home_dir/RepBlast_Temp");
	#Copy all full genome ffn files to RepBlast_Temp folder
	system("cp $Path_to_ffns/*.ffn $home_dir/RepBlast_Temp");
	my @ffn = glob("*.ffn");
	foreach(@ffn){system("formatdb -p F -i $_");}

#This goes into PAN_HASH{} and grabs the representive sequence
	foreach my $key (keys %{ $PAN_HASH }) {
		my $file_name = $key; #Don't want the > in the file name
		substr($file_name, 0, 1) = "";
		my @Value_Split = split(/\t/,$PAN_HASH->{$key});
		$num_values = scalar @Value_Split;
		if($num_values >= 2){
			open(SELF_BLAST,">$file_name.FFN_of_Rep_to_BLAST");
				print SELF_BLAST "$FFN_HASH->{$key}";
			close(SELF_BLAST);
		}#end if
	}#end foreach()
	
#Doing the same thing as Checker() except this time need to take the rep gene and blast it against each Prokka .ffn
	my @files = glob "*.FFN_of_Rep_to_BLAST";
	
	my $info;
    my $cpu;
    my $cpus;
    my $pm;
	#Get things ready to run in parallel
	if($threads == 0){
        $info=Sys::Info->new;
        $cpu  = $info->device( CPU => my %options );
        $cpus=scalar($cpu->identify);
        $pm=Parallel::ForkManager->new($cpus);
    }

    else{
        $pm=Parallel::ForkManager->new($threads);
    }

	RUN:
		foreach(@files){
			$pm->start and next RUN;
			my $name = $_;
			$name =~ s{\.[^.]+$}{}; #remove FFN_of_Rep_to_BLAST
			#Get list of genomes to run against:
			my $To_Run = $PAN_HASH->{"\>$name"};
			my @To_Run = split(/\t/,$To_Run);
			my @FFNs;
			#The array @To_Run has all the Pan_Genome family members for that given representative
			foreach(@To_Run){
				my $gen = $_;
				my @split = split(/\_/,$gen);
				$gen = $split[0];
				substr($gen, 0, 1) = "";
				push(@FFNs,$gen);
			}
			#This array was created from the @To_Do and it stripped the PROKKA gene number.
			#This way we can make the associated ##.ffn file the database for the blast against the representative of this Pan_Genome family
			my $db;
			foreach(@FFNs){
				$db = $_;
				#A low word size must be used because the THRESHOLDS are base on Amino Acid blasts.
				#Therefore, if the THRESHOLD is met with the Amino, when translated, there may be more mismatches with FFN
				system("blastn -query $name.FFN_of_Rep_to_BLAST -db $db.ffn -outfmt '6 qseqid stitle score' -word_size 5 > $name.$db.BLAST.old");
			}
			system("rm $name.FFN_of_Rep_to_BLAST");	#Can now delete the Pan_Genome family FFN file
			                                        #Can't delete $db.ffn since that is the full .ffn for each genome in the analysis. It is used for every blast			
			$pm->finish;
	}#end foreach(@files)
	$pm->wait_all_children;

	my @OLD = glob("*.BLAST.old");
	foreach(@OLD){
		my $file = $_;
		my $name = $file;
		#Must isolate the large ####.BLAST.old down to its PROKKA ID name that can be found in PAN_HASH and FNN_HASH and MASTER_HASH
		#Stip the .BLAST.old from the filename
		$name =~ s/\.BLAST\.old//;
		#Now need to get rid of everything after PROKKA number
		my @split = split(/\_/,$name);
		my @split2 = split(/\./,$split[1]);
		$name = $split[0] . "_" . "$split2[0]"; 
		my $db = $split[1];
		substr($db, 0, 6) = "";
		#print "NAME: $name \n DB: $db\n";
		open(BLAST,"<$name.$db.BLAST.old");
		open(OUT,">$name.$db.sort.BLAST");
		my %TEMP;
		#Here is where the BLAST.old files are sorted so the best hit is on top
		while(<BLAST>){
			chomp();
			my @split = split(/\t/,$_);
			my $key = "$split[2]\t$split[0]";
			my $value = "$split[0]\t$split[1]\t$split[2]";
			$TEMP{$key}=$value;
		}
		foreach my $key (sort {$b <=> $a} keys %TEMP){
			print OUT $TEMP{$key} . "\n";
		}
		close(OUT);
		close(BLAST);
	}
	system("rm *.BLAST.old"); #Only need to keep the sorted blast at this point
	system("rm *.nhr");
	system("rm *.nin");
	system("rm *.nsq");
	
		#Go through all elements (keys and values) or %PAN_HASH and set the values in %Master_Hash to done for all that appear.
		#foreach my $key (keys %{ $PAN_HASH }) {
		#	my @tab_split = split(/\t/,$PAN_HASH->{$key});
		#	push(@tab_split,$key);
		#	foreach(@tab_split){
		#		$MASTER_HASH->{$_} = "Done";
		#		
		#	}
		#}#end foreach
		
		#-- FRO DEBUG PURPOSES -- 
		#This will print MASTER_HASH{} and allow you to see the status on the use of each FFN gene
		#foreach my $key (keys %{ $MASTER_HASH }) {
        #    print "$key\t$MASTER_HASH->{$key}\n";
        #}
		
		#Want to have a file that has all the genes that have been replaced
		open(REPLACED, ">$home_dir/replacements.txt");
		
		#Now need to check first line to see if we are currently using the Ref with the best score
		my @files = glob "*.sort.BLAST";
		#Going through each sorted BLAST file seperately		
		foreach(@files){
			my $file_checker = $_; #Complete name of the file
			my $name = $_; #Want name to store the REP genome
			$name =~ s/\.sort\.BLAST//;
			my @split = split(/\_/,$name);
    		my @split2 = split(/\./,$split[1]);
	    	$name = $split[0] . "_" . "$split2[0]"; 
	    	my $db = $split[1];
	    	substr($db, 0, 6) = "";
	    	#print "NAME: $name \n DB: $db\n";
			open my $file, '<', "$name.$db.sort.BLAST";
			my $check_ref = <$file>;
			close($file);
			#This will retrieve the blast hit in the 2nd column of the report. 
			#It will cut out the predicted function
			my @check_ref = split(/\t/,$check_ref);
			my $Actual_Ref = $check_ref[1];
			my @space_split = split(" ",$Actual_Ref);
			#Here is the final version of $Actual_Ref : >#######_00000
			$Actual_Ref = "\>$space_split[0]";
			#print $name . "\t" . "$db" . "\t" . $Actual_Ref . "\n";
			my @split = split('_',$Actual_Ref);
			my $Ref_genome = $split[0];#This strips the prokka number
			#If the $Actual_Ref is not the key in the PAN_HASH then it is not the current Ref file
			#Check %Master_Hash as well to make sure the genome hasn't been included yet!
			
			#$final_value is outside the loop because it will store the new replacement
			my $final_value = "";
			#$count is used to keep track of where we are in the Pan_Genome family
			my $count =0;
			
			#If the values for the current file-name (i.e. the representitive of the current Pan_Genome family) already contain $Actual_Ref, we move on
			if($PAN_HASH->{"\>$name"} =~ m/$Actual_Ref/){next;}
			
			else{ #Use to have: elsif(!exists $PAN_HASH->{$Actual_Ref} && $MASTER_HASH->{$Actual_Ref} ne "Done" but I don't think it is necessarry
				#go through the values to find the one that needs to be replaced
				#$MASTER_HASH->{$Actual_Ref} = "Done";
				#This split will give us all the genes in the Pan_Genome family of the current Ref/Key $name
				my @tab_split = split(/\t/,$PAN_HASH->{"\>$name"}); 
				foreach(@tab_split){
					#print $_ . "\n\n";
					my $genome = $_;#The current genome in the loop in the Pan_Genome family
					substr($genome, 0, 1) = ""; #Get rid of the '>'
					my $Ref_genome_pure = $Ref_genome;
					substr($Ref_genome_pure, 0, 1) = "";
					my @undscore_split = split(/\_/,$genome);

                    #print "UNDERSCORE:$undscore_split[0]:\tREF_GENOME_PURE:$Ref_genome_pure:ACTUAL:$Actual_Ref:\tNAME:$name:\tMASTER_HASH:$MASTER_HASH->{$Actual_Ref}:\n";

                    #Here we are trying to make the replacement:
                    #In the first case, if the first genome in the values of PAN_HASH{$name} is the genome we are after
					if($undscore_split[0] eq $Ref_genome_pure && $count == 0){
						$final_value = "$Actual_Ref" . "\t";
						print REPLACED "$genome\treplaced by\t$Actual_Ref\n";
					}
					#Here the genome we are after is 
					elsif($undscore_split[0] eq $Ref_genome_pure && $count != 0){
						$final_value = $final_value . "$Actual_Ref" . "\t";
						print REPLACED "$genome\treplaced by\t$Actual_Ref\n";
					}
					elsif($count == 0){
						$final_value = $_ . "\t";
					}
					else{
						$final_value .= $_ . "\t";
					}
					$count++;
				}#end foreach
				$PAN_HASH->{"\>$name"}=$final_value;
			}#end if(exists)
		}#end foreach(@files)
		close(REPLACED);	
		chdir($home_dir);
		open(HASH,">PanGenome.usearch.reFFNed");
			while (my ($key) = each(%$PAN_HASH)) {
				print HASH $PAN_HASH->{$key} . "\n";
			}
		close(HASH);
		system("rm -r RepBlast_Temp");
}#end sub checker()

sub Rename_Fix{
	my $home_dir = getcwd;
	system("mv PanGenome.usearch PanGenome.usearch.backup");

	open(IN,"<PanGenome.usearch.reFFNed") or die "Cannot open PanGenome.usearch.reFFNed";
	open(OUT,">PanGenome.usearch");

	while(<IN>){
		chomp();
		my $line = $_;
		my @tab_split = split(/\t/,$line);
		my $final="";
		foreach(@tab_split){
			my $genome = $_;
			my @underscore_split = split(/\_/,$genome);
			substr($underscore_split[0], 0, 1) = "";
			$final .= "$genome\[$underscore_split[0]\]\t";
		}
		print OUT $final . "\n";
	}#end while(<IN>)	


	close(OUT);
	close(IN);
	system("rm ALL.ffn");
	#system("rm PAN_HASH_TEST");
}#sub Rename_Fix


my $Path_to_PanGenome = $ARGV[0];
my $Path_to_ffns = $ARGV[1];

my $threads = $ARGV[2];

Checker($Path_to_PanGenome, $Path_to_ffns, $threads);
Rename_Fix();


#my $FFN_HASH = FFN_Hash($Path_to_ffns);
#my $PAN_HASH = PanGenomeHash($Path_to_PanGenome);
#my $MASTER_HASH = Master_Hash();

#A quick check to see if PanGenome Hash Reference is populated
#foreach my $key (keys %{ $PAN_HASH }) {
#    print $key . "\n" . $PAN_HASH->{$key} . "\n";
#}
#print "\n\n";
#A quick check to see if FFN Hash Reference is populated
#foreach my $key (keys %{ $FFN_HASH }) {
#    print $FFN_HASH->{$key};
#}
#print "\n\n";
#A quick check to see if MASTER Hash Reference is populated
#foreach my $key (keys %{ $MASTER_HASH }) {
#    print "$key\t$MASTER_HASH->{$key}\n";
#}
