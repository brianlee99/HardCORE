#! /usr/bin/perl
use Algorithm::Combinatorics qw( combinations );
use Parallel::ForkManager;
use MCE::Shared;
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
	my @ALL_check = glob "*.ffn";
	if($ALL_check[0] ne "ALL.ffn"){system("cat $path/*.ffn >> ALL.ffn");}
	else{print "\nNo need to generate ALL.ffn\n";}
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

    tie my %SELF_BLAST, "MCE::Shared";

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

    #Going to place all BLAST outputs into a hash using the backticks (``)
    #The key of the hash is the name of the file without the ".SELF_BLAST" extension
  
		
	RUN:
		foreach(@files){
			$pm->start and next RUN;
			my $name = $_;
			$name =~ s{\.[^.]+$}{};
			system("formatdb -p F -i $_");
			#Here is the self-blast
			#Use a variable to store the STDOUT of BLAST
			my $blast_out = `blastn -query $_ -db $_ -outfmt '6 qseqid stitle score' -word_size 15`;
			#Put the output into the %SELF_BLAST
			$SELF_BLAST{$name} = $blast_out;
			system("rm $_.nhr");
			system("rm $_.nin");
			system("rm $_.nsq");
			system("rm $_"); #This removes the FFN file with all the genes from the Pan_Genome family
			#The ####.BLAST is full of self-hits, need to be rid of them
		   #Here is where we get rid of self-hits and save them in #####.sort.BLAST
		   	my %TEMP;
			my @line_split = split(/\n/,$SELF_BLAST{$name});
			foreach(@line_split){
				my @tab_split = split(/\t/,$_);
					if($tab_split[0] ne $tab_split[1]){
						my $key = "$tab_split[2]\t$tab_split[0]";
						my $value = "$tab_split[0]\t$tab_split[1]\t$tab_split[2]";
						$TEMP{$key}=$value;
					}
			}
			my $new_value;
			foreach my $key (sort {$b <=> $a} keys %TEMP){
		 		$new_value .= $TEMP{$key} . "\n";
			}
			undef %TEMP;
			$SELF_BLAST{$name}= $new_value;
			$pm->finish;
		}#end foreach(@files)
		$pm->wait_all_children;
		    
	##BIG STEP## 
	#Now need to check first line to see if we are currently using the Ref with the best score
			foreach my $key (keys %SELF_BLAST){
				my $Actual_Ref = "";
				my $name = $key;
				my @line_split = split(/\n/,$SELF_BLAST{$key}); #split that value on a per hit basis
				foreach(@line_split){
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
   			    }#end foreach(@line_split)
			}#end foreach(keys %SELF_BLAST)		
		chdir($home_dir);
		system("rm -r Self_Blast_Temp"); #Can now get rid of the folder
	#Ths hash print would allow you to see the KEY for each Pan_Genome family as well as its values on the next line	

		#Remove hash from memory
		undef %SELF_BLAST;

		#open(HASH,">PAN_HASH_TEST");
		#while (my ($key) = each(%$PAN_HASH)) {
		#	print HASH $key . "\n" . $PAN_HASH->{$key} . "\n";
		#}
		#close(HASH);


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

    tie my %REP_BLAST, "MCE::Shared"; #Like above this is a shared hash that will store all the blast records

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
			my $To_Run = $PAN_HASH->{"\>$name"}; #Get the whole pan-genome family for this representitive
			my @To_Run = split(/\t/,$To_Run); #This is an array of the enitre pan genome family
			my @FFNs;
			#The array @To_Run has all the Pan_Genome family members for that given representative. 
			#@FFN only has the genome name (no info on the gene, this portion is taken out in this loop)
			#@FFN is going to be used to in the foreach(@FFNs) loop. The representative is going to be BLASTed against each genome family member
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
				my $db_blastout = "__" . $db; #There are 2 underscores here to help with splitting!
				#A low word size must be used because the THRESHOLDS are base on Amino Acid blasts.
				#Therefore, if the THRESHOLD is met with the Amino, when translated, there may be more mismatches with FFN
				#OLD CALL: system("blastn -query $name.FFN_of_Rep_to_BLAST -db $db.ffn -outfmt '6 qseqid stitle score' -word_size 5 > $name.$db.BLAST.old");
				my $out = `blastn -query $name.FFN_of_Rep_to_BLAST -db $db.ffn -outfmt '6 qseqid stitle score' -word_size 5`;
				#Sort the output of $out first:
				my %TEMP; #temporary hash used to sort
				my @line_split = split(/\n/,$out);
				my $REP_BLAST_out = "$name" . "$db_blastout"; #This looks like: Newport.2797.BMH_04460__Montevideo.2005.BMH IT IS THE KEY for %REP_BLAST
				#my @split_names = split('__',$REP_BLAST_out);
				#my $rep = $split_names[0];
				#my $db = $split_names[1];
				foreach(@line_split){
					chomp($_);
					my @split = split(/\t/,$_);
					my $key = "$split[2]\t$split[0]";
					my $value = "$split[0]\t$split[1]\t$split[2]";
					$TEMP{$key}=$value;
				}
				my $new_value;
				foreach my $key (sort {$b <=> $a} keys %TEMP){
					$new_value .= $TEMP{$key} . "\n";
				}
				$REP_BLAST{$REP_BLAST_out}= $new_value;	#Placing the BLAST results in order of highest score back into %REP_BLAST		
				undef %TEMP;
			}
			system("rm $name.FFN_of_Rep_to_BLAST");	#Can now delete the Pan_Genome family FFN file
			                                            #Can't delete $db.ffn since that is the full .ffn for each genome in the analysis. It is used for every blast			
			$pm->finish;
	}#end foreach(@files)
	$pm->wait_all_children;

	#Get rid of db files
	system("rm *.nhr");
	system("rm *.nin");
	system("rm *.nsq");
	system("rm formatdb.log");

	#For debug purposes; Print out the %REP_BLAST hash which is all the blast results from the previous loop
	#open(OUT,">REP_BLAST.hash");
	#foreach my $key (keys %REP_BLAST){
	#	print OUT "KEY:$key\n$REP_BLAST{$key}\n";
	#}
	#close(OUT);

	
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
		foreach my $key (keys %REP_BLAST){
			my @key_split = split('__',$key);
			my $Rep_genome = @key_split[0];
			my $name = $Rep_genome; #Needed to re-purpose old code which used this name variable
			my $Blasted_against = @key_split[1];
			#Need to retrieve the best blast hit for the rep and cut out the predicted function
			my @line_split = split(/\n/,$REP_BLAST{$key});
			my $Top_Blast_Line = $line_split[0];
			my @tab_split = split(/\t/,$Top_Blast_Line);
			my $Top_hit = $tab_split[1];
			my @space_split = split(' ',$Top_hit);
			my $Actual_Ref = "\>$space_split[0]"; #The predicted function is now gone, need the ">" for hash purposes
			my @underscore_split = split('_',$Actual_Ref);
			my $Ref_genome = $underscore_split[0];#This strips the prokka number
			#If the $Actual_Ref is not the key in the PAN_HASH then it is not the current Ref
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
					my $Ref_genome_pure = $Ref_genome; #pure signifies no ">"
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
		}#end foreach my $key (keys %REP_BLAST)

		close(REPLACED);	
		chdir($home_dir);
		open(HASH,">PanGenome.usearch.reFFNed");
			while (my ($key) = each(%$PAN_HASH)) {
				print HASH $PAN_HASH->{$key} . "\n";
			}
		close(HASH);
		undef %REP_BLAST;
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
