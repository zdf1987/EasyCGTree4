#!/usr/bin/perl
use warnings;
use strict; 
no warnings 'experimental::smartmatch';
use File::Copy qw(move mv);
use Getopt::Long;


my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $rmon=$mon +1;
my $ryear=$year+1900;

print "\n##############################\nReading command line... checking options and input data...\n\n";##########################注意更新版本和日期
my @usage=qq(
====== EasyCGTree ======
     Version 4.0 by Dao-Feng Zhang (Hohai University, China)
	 Update 2022-09-08

Usage: perl EasyCGTree.pl [Options]

Options:
-proteome <String>	(Essential)
	Input data (proteomes) directory
	
-task <String, 'all', 'hmmsearch', 'refine', 'alignment', 'tree_infer'> (Optional)
	Set run mode. [default: all]
-tree <String, 'sm', 'st'>(Optional)				
	Approach (sm, supermatrices; st, supertree) used for tree inference. [default: sm]
-tree_app <String, 'fasttree', 'iqtree'> (Optional)
	Applications used for tree inference. [default: fasttree]
-hmm <String, 'bac120', 'rp1', et al.> (Optional)
	Profile HMM used for gene calling. [default: bac120]
-thread <Int> (Optional)		
	Number of threads to be used by hmmer, clustalo (Linux) and iqtree. [default: 2]
-trim <String, 'gappyout','strict', 'strictplus'> (Optional)		
	Standard for trimming alignment used by trimAl. [default: strict]
-evalue <Real, 10..0> (Optional)
	Expect value for saving hits during . [default: 1e-10]
-gene_cutoff <Decimal, 0.5..1> (Optional)
	Cutoff for omitting low-prevalence gene. [default: 0.8]
-genome_cutoff <Decimal, 0.5..1> (Optional)
	Cutoff for omitting low-quality genomes. [default: 0.8]
	
-help (Optional)				
	Display this message.

);
my ($inputDir,$user_genesDir,$check1,$check2,$check3,@query,@queryR);


########## Default parameters ################
my $task="all";
my $evalue=1e-10;
my $geneCutoff=0.8;
my $hmm = "bac120";
my $trim="strict";
my $genomeCutoff=0.8;
my $thread=2;
my $tree_type ="sm";
my $tree_app = "fasttree";
##########################################

my %opt=qw();
GetOptions(\%opt,"proteome:s","task:s","tree:s","tree_app:s","hmm:s","trim:s","thread:i","evalue:f","gene_cutoff:f","genome_cutoff:f","help!");

if (scalar(keys %opt )==0 || exists($opt{"help"})) {
	print join("\n",@usage)."\n\n";
	exit;
}

unless (exists($opt{"proteome"})) {
	print join("\n",@usage)."\n\nERROR: '-proteome' must be set to start.\n\n";
	exit;	
}	

my $fn1;
{
	$inputDir = $opt{"proteome"};
	unless ($inputDir=~/\w/ || $inputDir=~/\d/) {
		print join("\n",@usage)."\n\nERROR: Argument '-proteome $inputDir': it should be a directory name, e.g., myGenome.\n\n";
		exit;
	}
	&rename($inputDir); 
	$check1= &checkDataType($inputDir);
	if ($check1 eq "nucl") {
		print join("\n",@usage)."\n\nERROR: Argument '-proteome $inputDir': directory '$inputDir' should contains file(s) of protein sequence, while some file(s) of nucleotide sequence were found.\n\n";
		exit;
	} elsif ($check1 eq "prot") {
	} elsif ($check1=~ /\.fas/ || $check1=~ /open input/) {
		print join("\n",@usage)."\n\nERROR: Argument '-proteome $inputDir': $check1\n\n";
		exit;
	} else {
		print join("\n",@usage)."\n\nERROR: Argument '-proteome $inputDir': directory '$inputDir' should only contains fasta-formated file(s) of protein sequence, while $check1\n\n";
		exit;		
	}
	opendir(DIR, $inputDir);
	@query = readdir(DIR);
	closedir(DIR);
	splice @query, 1, 1 unless $query[1]=~ /\.fas/;
	splice @query, 0, 1 unless $query[0]=~ /\.fas/;
	$fn1=$#query+1;
	unless ($fn1 >=5) {
		
		print join("\n",@usage)."\n\nERROR: Number of files (taxa on a tree) in Input directory '$inputDir' ($fn1) must be >=5 to start a run.\n\n";
		exit;
	}
}
	
my @tasks= qw(all hmmsearch refine alignment tree_infer);

if (exists($opt{"task"})) {
	$task = $opt{"task"};
	unless ($task ~~ @tasks) {
			
		print join("\n",@usage)."\nERROR: Argument '-task $task': the value '$task' should be one of the following: all, hmmsearch, refine, alignment, or tree_infer.\n\n";
		exit;
	}
}

$evalue=$opt{"evalue"} if exists($opt{"evalue"});

my $treesig;
if (exists($opt{"tree"})) {
	$tree_type = $opt{"tree"};
	if ($tree_type eq "st") {
		$treesig=1;
		
	} elsif ($tree_type eq "sm") {
	
	} else {
		print join("\n",@usage)."\nERROR: Argument '-tree $tree_type': the value '$tree_type' cannot be recognized (requested to be 'st' or 'sm').\n\n";
		exit;
	}
} 

if (exists($opt{"trim"})) {
	$trim = $opt{"trim"};
	if ($trim eq "gappyout") {
	} elsif ($trim eq "strict") {
	} elsif ($trim eq "strictplus") {
	} else {
		print join("\n",@usage)."\nERROR: Argument '-trim $trim': the value '$trim' cannot be recognized (requested to be 'gappyout','strict', or 'strictplus').\n\n";
		exit;
	}
} 

if (exists($opt{"tree_app"})) {
	$tree_app = $opt{"tree_app"};
	unless ($tree_app eq "fasttree" || $tree_app eq "iqtree") {
		print join("\n",@usage)."\nERROR: Argument '-tree_app $tree_app': the value '$tree_app' cannot be recognized (requested to be 'fasttree' or 'iqtree').\n\n";
		exit;
	}
} 

if (exists($opt{"hmm"})) {
	$hmm = $opt{"hmm"};
	my $ex="./HMM/$opt{'hmm'}.hmm";
	unless (-e $ex) {
		print join("\n",@usage)."\nERROR: $ex Argument '-hmm $hmm': the profile HMM '$hmm' cannot be found in the directory './HMM'.\n\n";
		exit;
	}
} 



if (exists($opt{"genome_cutoff"})) {
	my $qwe=$opt{"genome_cutoff"};
	if ($opt{"genome_cutoff"} =~ /^0?\.?\d+$/) {
		if ($opt{"genome_cutoff"}>=0.5 && 1>=$opt{"genome_cutoff"}) {
			$genomeCutoff= $opt{"genome_cutoff"};
		} else {
			print join("\n",@usage)."\nERROR: Argument '-ggenome_cutoff $qwe': the value '$qwe' does not fall in recommended interval 0.5-1.\n\n";
			exit;
		}
	} else {
		print join("\n",@usage)."\nERROR: Argument '-genome_cutoff $qwe': '$qwe' is not a decimal.\n\n";
		exit;
	}
}	
	
my $BNM;
if (exists($opt{"thread"})) {	
	$BNM=$opt{"thread"};
	if ($opt{"thread"} =~ /^\d+$/) {
			$thread=$opt{"thread"};
	} else {
		print join("\n",@usage)."\nERROR: Argument '-thread $BNM': '$BNM' is not an integer.\n\n";
		exit;
	}
}



open (LOG, ">$inputDir\_$ryear-$rmon-$mday-$hour-$min.log") or die "Can't open '$inputDir\_$ryear-$rmon-$mday-$hour-$min.log': $!\n";

my $order;
my @nouseop;
if ($task eq "all") { 
	$order= "The user ordered a complete analysis (-task all) of the pipeline!\n\n";
	
} elsif ($task eq "hmmsearch") {
	$order=  "The user ordered a HMM search task (-task hmmsearch) of the pipeline!\n\n";
	my @need= qw(hmm thread evalue);
	my @noneed= qw(tree tree_app gene_cutoff genome_cutoff trim);
	foreach my $op (keys %opt) {
		next if $op ~~ @need;
		push @nouseop, $op if $op ~~ @noneed;
	}
} elsif ($task eq "refine") {
	$order=  "The user ordered a task (-task refine) of screening significant homologs from the pre-existing HMM searching result and generating gene clusters!\n\n";
	my @need= qw(gene_cutoff genome_cutoff);
	my @noneed= qw(hmm thread evalue tree tree_app trim);
	foreach my $op (keys %opt) {
		next if $op ~~ @need;
		push @nouseop, $op if $op ~~ @noneed;
	}
	
} elsif ($task eq "alignment") {
	$order= "The user ordered a task (-task alignment) of aligning sequences of each gene cluster!\n\n";
	my @need= qw(thread trim);
	my @noneed= qw(hmm evalue tree tree_app gene_cutoff genome_cutoff);
	foreach my $op (keys %opt) {
		next if $op ~~ @need;
		push @nouseop, $op if $op ~~ @noneed;
	}	
} elsif ($task eq "tree_infer") {
	$order= "The user ordered a task (-task tree_infer) of infering phylogenetic tree from pre-existing alignment(s)!\n\n";
	my @need= qw(thread tree tree_app);
	my @noneed= qw(hmm evalue gene_cutoff genome_cutoff trim);
	foreach my $op (keys %opt) {
		next if $op ~~ @need;
		push @nouseop, $op if $op ~~ @noneed;
	}		
}



if (@nouseop) {
	print "\nWARNING: '-task $task' dose not require the following settings, which will be ignored in current analysis:";
	print LOG "\nWARNING: '-task $task' dose not require the following settings, which will be ignored in current analysis:";
	foreach my $zz (@nouseop) {
		print " -$zz $opt{$zz}";
		print LOG " -$zz $opt{$zz}";
	}
	print ".";
	print LOG ".";
}

my $rewq;
if (exists($opt{"gene_cutoff"}))	 {
	$rewq=$opt{"gene_cutoff"};
	if ($opt{"gene_cutoff"} =~ /^0?\.?\d+$/) {
		if ($opt{"gene_cutoff"}>=0.5 || 1>=$opt{"gene_cutoff"}) {
			$geneCutoff= $opt{"gene_cutoff"};
		} else {
			print join("\n",@usage)."\nERROR: Argument '-gene_cutoff $rewq': the value '$rewq' does not fall in recommended interval 0.5-1.\n\n";
			exit;
		}
	} else {
		print join("\n",@usage)."\n\nERROR: Argument '-gene_cutoff $rewq': '$rewq' is not a decimal.\n\n";
		exit;
	}
}
if ($treesig) {
	print "\nWARNING: '-tree $tree_type' requires '-gene_cutoff' to be set as '1'.\nThe value '1' will be used instead.\n" if $geneCutoff <1;
	print LOG "\nWARNING: '-tree $tree_type' requires '-gene_cutoff' to be set as '1'.\nThe value '1' will be used instead.\n" if $geneCutoff <1;	
	$geneCutoff=1;
} 



print "\n\n###############Starting Informations ###############\n\n";
print LOG "\n\n###############Starting Informations ###############\n\n";

print LOG '====== EasyCGTree ======
	Version 4.0
by Dao-Feng Zhang (Hohai University, China)

';
print '====== EasyCGTree ======
	Version 4.0
by Dao-Feng Zhang (Hohai University, China)

';

print LOG "#########################\n\nOptions: \n";  
print "#########################\n\nOptions: \n";

print "-proteome $inputDir\n-task $task\n-tree $tree_type\n-tree_app $tree_app\n-trim $trim\n-hmm $hmm\n-thread $thread\n-evalue $evalue\n-gene_cutoff $geneCutoff\n-genome_cutoff $genomeCutoff\n";
print LOG "-proteome $inputDir\n-task $task\n-tree $tree_type\n-tree_app $tree_app\n-trim $trim\n-hmm $hmm\n-thread $thread\n-evalue $evalue\n-gene_cutoff $geneCutoff\n-genome_cutoff $genomeCutoff\n";

print LOG "\n#########################\n\nJob Started at: $hour:$min:$sec,$ryear-$rmon-$mday\n\n*****$order\n\n";  
print "\n#########################\n\nJob Started at: $hour:$min:$sec,$ryear-$rmon-$mday\n\n*****$order\n\n";

my $TEMdir= $inputDir."_TEM";

############### Task 1: HMM search #################
#============= Step 1: HMM search ==============# 
if ($task eq "all" || $task eq "hmmsearch") {
	print "\n############### Task 1: HMM search ###############\n#============= Step 1: HMM search ==============# \n";
	print LOG "\n############### Task 1: HMM search ###############\n#============= Step 1: HMM search ==============# \n";
	if (-e $TEMdir) {
		my $del1="rm -rf $TEMdir";                ##############Different from Windows version
		system ($del1); 
		die "Step 1: Directory '$TEMdir' already exist, but permission denied when delete it.\n" if (-e "$TEMdir"); 
	}
	
	mkdir "$TEMdir" || die " Permission denied to create directories. Please contact the Administrator of your System\n"; 
	
	mkdir "$TEMdir/TEM1_HMMsearch_out" || die " Permission denied to create directories. Please contact the Administrator of your System\n";
	my $genomeNum = 0;
	
	open(HMM, "<HMM/$hmm.hmm") or die "Step 1: Can't open 'HMM/$hmm.hmm': $!\n";
	open(MMIN, ">$TEMdir/HMMinfo.txt") or die "Step 1: Can't open '$TEMdir/HMMinfo.txt': $!\n";
	foreach my $hmmin (<HMM>) {
		$hmmin=~ s/\n// unless $hmmin=~s/\r\n//;
		next unless $hmmin;
		if ($hmmin =~ s/NAME\s+//) {
			print MMIN "$hmmin	";
		} elsif ($hmmin =~ s/LENG\s+//) {
			print MMIN "$hmmin\n";
		}
	}
	close HMM;
	close MMIN;
	
	print "Step 1: Do HMM search against $fn1 files in directory '$inputDir': \n\n";
	foreach my $file (@query) {

		$genomeNum++;
		my $cmd;
		
		$cmd = "./bin/hmmsearch --tblout ./$TEMdir/TEM1_HMMsearch_out/$file -E $evalue --cpu $thread ./HMM/$hmm.hmm ./$inputDir/$file";    ###############Different from Windows version
		
		my $response1= `$cmd` || die "Step 1: Can't execute 'hmmsearch': $!\nPlease make sure: 'hmmsearch' is present in directory 'bin'.\n";  

		my $rat= &givePercentage($genomeNum,$fn1);

		
		print ("\r","	Searching against './$inputDir/$file': $genomeNum/$fn1 ($rat).           ");
	}
	print "\n\n";
	print "=====Results of $genomeNum HMM search were wroten in directory './$TEMdir/TEM1_HMMsearch_out/'!=====\n";
	print LOG "=====Results of $genomeNum HMM search were wroten in directory './$TEMdir/TEM1_HMMsearch_out/'!=====\n";
}
	

############### Task 2: Filtrate HMM Search Results #################
#============= Step 2: Filtrate HMM Search Results ==============# 
if ($task eq "all" || $task eq "refine") {

	print "\n\n############### Task 2: Filtrate HMM Search Results ###############\n#============= Step 2: Filtrate HMM Search Results ==============# \n";
	print LOG "\n\n############### Task 2: Filtrate HMM Search Results ###############\n#============= Step 2: Filtrate HMM Search Results ==============# \n";
	               
	if (-e "$TEMdir/TEM2_HMMsearch_outS") {########Different from Windows version
		my $del2="rm -rf $TEMdir/TEM2_HMMsearch_outS"; ###########Different from Windows version
		system ($del2);
		die "Step 2: Directory '$TEMdir/TEM2_HMMsearch_outS' already exist, but permission denied when delete it.\n" if (-e "$TEMdir/TEM2_HMMsearch_outS"); 
	}
	mkdir "$TEMdir/TEM2_HMMsearch_outS";
	
	opendir(DIR3, "$TEMdir/TEM1_HMMsearch_out") || die "Step 1: Can't open input directory '$TEMdir/TEM1_HMMsearch_out': $!\n";
	my @query2 = readdir(DIR3);
	closedir(DIR3);
	splice @query2, 1, 1 unless $query2[1]=~ /\.fas/;
	splice @query2, 0, 1 unless $query2[0]=~ /\.fas/;
	die "Step 3: Input  directory $TEMdir/TEM1_HMMsearch_out does not contain any files.\n" unless scalar(@query2);
	
	open(MMIN, "<$TEMdir/HMMinfo.txt") or die "Step 1: Can't open '$TEMdir/HMMinfo.txt': $!\n";
	my %score;
	foreach my $hmmin (<MMIN>) {
		$hmmin=~ s/\n// unless $hmmin=~s/\r\n//;
		next unless $hmmin;
		my @hmmi =split /\t/, $hmmin;
		$score{$hmmi[0]} = $hmmi[1];
	}
	close MMIN;
	
	
	foreach my $file3 (@query2) {
		my (%out3,%score3);
		next unless $file3 =~ /\.fas$/;                    
		open (FIL3, "<$TEMdir/TEM1_HMMsearch_out/$file3") or die "Can't open '$TEMdir/TEM1_HMMsearch_out/$file3': $!\n";
		my (@key3,@result3);
		my $geneC=0;
		my $head;
		foreach my $in3 (<FIL3>) {
			$in3=~ s/\n// unless $in3=~s/\r\n//;
			next unless $in3;
			if ($in3 =~ /\#/) {
				$head .= "$in3\n" if $head;
				$head = "$in3\n" unless $head;
				next;
			}
			my @in3 =split /\s+/, $in3;                           
			
			next if $in3[4] > $evalue;
			
			my $sco = $score{$in3[2]}/4;      #####Drop it? If not, how to use it?
			next if $in3[5] < $sco;
			my $keyr;
			if ($in3[3] =~ /-/) {
				$keyr =$in3[2];
			} else {
				$keyr =$in3[3];
			}
			push @key3, $keyr;									
			unless ($out3{$keyr}) {
			
				$out3{$keyr}=$in3;
				$geneC++;
			} 
		}
		close FIL3;
		my $outfile3="$geneC"."__$file3";
		open(OU3, ">$TEMdir/TEM2_HMMsearch_outS/$outfile3") or die "Step 3: Can't open '$TEMdir/TEM2_HMMsearch_outS/$outfile3': $!\n";
		print OU3 "$head";
		foreach my $in (keys %out3) {
			print OU3 "$out3{$in}\n";
		}
		close OU3;
	}
	print "=====Screened HMM search results has been written to '$TEMdir/TEM2_HMMsearch_outS'!=====\n";
	print LOG "=====Screened HMM search results has been written to '$TEMdir/TEM2_HMMsearch_outS'!=====\n";

###################################################
	opendir(DIR, "$TEMdir/TEM2_HMMsearch_outS") || die "Step 3: Can't open input directory '$TEMdir/TEM2_HMMsearch_outS': $!\n";
	my @query3 = readdir(DIR);
	closedir(DIR);
	die "Task 2: Input directory TEM2_HMMsearch_outS does not contain any files.\n" unless scalar(@query3);

	my (@querylist0,@querylistN);
	my (@finalgenelist,@outgenelist,@finalgenome,@outgenome);
	my (@finalfile);
	my $geneNum=0;
	my $genomeNum;

	foreach my $kk (@query3) {
		next unless $kk =~ /.\.fas$/;
		$genomeNum++;
		my @yy=split /__/,$kk;
		$geneNum = $yy[0] if $geneNum < $yy[0]
	}
	my $cutoff2=int($geneNum*$genomeCutoff);     ####################
	print "\n These genomes harboring more than $cutoff2 ($geneNum * $genomeCutoff) genes, which will be used in following analysis. This is the list (genome: gene_number):\n";
	print LOG "\n These genomes harboring more than $cutoff2 ($geneNum * $genomeCutoff) genes, which will be used in following analysis. This is the list (genome: gene_number):\n";

	foreach my $kk (@query3) {
		next unless $kk =~ /.\.fas$/;
		my @yy=split /__/,$kk;
		if ($yy[0]>=  $cutoff2) {
			push @finalfile, $kk ;
			push @finalgenome, "$yy[1]:$yy[0]";
			print "$yy[1]: $yy[0]\n";
			print LOG "$yy[1]: $yy[0]\n";
		} else {
			push @outgenome, "$yy[1]: $yy[0]\n";
		}
	}

	my $dis=$#outgenome+1;
	my $ingenomeN=$#finalgenome+1;	

	print ".\n=====Totally, $ingenomeN of $genomeNum genomes were selected!=====\n\n";
	print LOG ".\n=====Totally, $ingenomeN of $genomeNum genomes were selected!=====\n\n";
	if (@outgenome) {
	print "The following $dis genomes were excluded(genome/gene_number): @outgenome.\n\n";
	print LOG "The following $dis genomes were excluded (genome/gene_number): @outgenome.\n\n";
	}

	my $num=0;
	foreach my $file3 (@finalfile) {
		open (FIL3, "<$TEMdir/TEM2_HMMsearch_outS/$file3") or die "Step 3: Can't open '$TEMdir/TEM2_HMMsearch_outS/$file3': $!\n";
		foreach my $in3 (<FIL3>) {
			$in3=~ s/\n// unless $in3=~s/\r\n//;
			next unless $in3;
			next if $in3 =~ /\#/;
			my @in3 =split /\s+/, $in3;
			my $keyr;
			if ($in3[3] =~ /-/) {
				$keyr =$in3[2];
			} else {
				$keyr =$in3[3];
			}
			unless (@querylist0) {
				push @querylist0,$keyr;
				$querylistN[$num]++;
				$num++;	
			} else {
				unless ($keyr ~~ @querylist0) {
					push @querylist0,$keyr;
					$querylistN[$num]++;
					$num++;
				} else {
					foreach my $inn3 (0..$#querylist0) {
						if ($keyr eq $querylist0[$inn3]) {
					
							$querylistN[$inn3]++;	
						}
					}
				}
			}
		}
		close FIL;
	}
	my $cutoff=int($ingenomeN*$geneCutoff);           ##############

	print "\nThese genes were present in >= $cutoff ($ingenomeN * $geneCutoff) of the $ingenomeN selected genomes, which will be used for tree inference (gene/prevalence):\n";
	print LOG "\nThese genes were present in >= $cutoff ($ingenomeN * $geneCutoff) of the $ingenomeN selected genomes, which will be used for tree inference (gene/prevalence):\n";

	foreach my $kk (0..$#querylistN) {
		if ($querylistN[$kk]<  $cutoff) {
			push @outgenelist, "$querylist0[$kk]: $querylistN[$kk]\n";
		} else {
			print "$querylist0[$kk]: $querylistN[$kk]\n";
			print LOG "$querylist0[$kk]: $querylistN[$kk]\n";
			push @finalgenelist, $querylist0[$kk];
		}
	}
	my $ingene =$#finalgenelist+1;
	my $outgene =$#outgenelist+1;

	print ".\n=====Totally, $ingene of $geneNum genes were selected!=====\n";
	print LOG ".\n=====Totally, $ingene of $geneNum genes were selected!=====\n";
	if (@outgenelist) {
	print "\nThe following $outgene/$geneNum genes were excluded (gene/prevalence):\n@outgenelist\n";
	print LOG "\nThe following $outgene/$geneNum genes were excluded (gene/prevalence):\n@outgenelist\n";
	}
	open (TXT3, ">$TEMdir/GenomeGeneScreened.txt") or die "Step 3: Can't open '$TEMdir/GenomeGeneScreened.txt': $!\n";
	print TXT3 "Genome_List=@finalgenome\n\nGene_List=@finalgenelist";
	print LOG "The Genome and genes used in following analysis were also listed in '$TEMdir/GenomeGeneScreened.txt'\n";
	print "The Genome and genes used in following analysis were also listed in '$TEMdir/GenomeGeneScreened.txt'\n";
	close TXT3;
}

############### Task 3: Multiple Sequences Alignment #################  
#============= Step 3: Retrieve Sequences from Each proteome ==============# 
if ($task eq "all" || $task eq "alignment") {

	print "\n\n############### Task 3: Multiple Sequences Alignment #################\n#============= Step 3: Retrieve Sequences from Each proteome ==============#\n\n";
	print LOG "\n\n############### Task 3: Multiple Sequences Alignment #################\n#============= Step 3: Retrieve Sequences from Each proteome ==============#\n\n";

	my $numg=0;
	
	if (-e "$TEMdir/TEM3_GeneSeqs") {################Different from Windows version
		my $del2="rm -rf $TEMdir/TEM3_GeneSeqs";  ################Different from Windows version
		system ($del2);
		die "Step 3: Directory '$TEMdir/TEM3_GeneSeqs' already exist, but permission denied when delete it.\n" if (-e "$TEMdir/TEM3_GeneSeqs"); 
	} 
	

	mkdir "$TEMdir/TEM3_GeneSeqs";
	open (TXT4, "<$TEMdir/GenomeGeneScreened.txt") or die "Step 3: Can't open '$TEMdir/GenomeGeneScreened.txt': $! If it is absent, you should run 'task refine' again.\n";
	my (@fgenome,@fgene);
	foreach my $k3 (<TXT4>) {
		$k3=~ s/\n// unless $k3=~ s/\r\n//;	
		next unless $k3;		
		if ( $k3 =~ s/Genome_List=//) {
			@fgenome = split /\s/, $k3;
		} elsif( $k3 =~ s/Gene_List=//) {
			@fgene = split /\s/, $k3;
		}
	}
	close TXT4;
	my $fg=$#fgenome +1;
	my $fgn=$#fgene +1;
	print "According to the record in '$TEMdir/GenomeGeneScreened.txt', $fgn common genes will be extracted from $fg genomes.\n";
	print LOG "According to the record in '$TEMdir/GenomeGeneScreened.txt', $fgn common genes would be extracted from $fg genomes.\n";

	foreach my $i1 (0..$#fgenome) {
		$numg++;
		my @hh=split /\:/, $fgenome[$i1];
		my $fgenome =$hh[0];
		my %geno;
		my ($seq,$seqID);

		open (FILE, "<$inputDir/$fgenome") or die "Step 4: Can't open '$inputDir/$fgenome': $!\n";
		
		foreach my $i3 (<FILE>) {
			$i3=~ s/\n// unless $i3=~ s/\r\n//;
			next unless $i3;		
			if ( $i3 =~ s/\>//) {
				if ($seq) {
					$geno{$seqID} = $seq;
					undef $seq;
				}
				my @ii = split /\s/, $i3;
				$seqID = $ii[0];			
			} else { 
				$i3 =~ s/\*/-/g;
				$seq .= $i3 if $seq;
				$seq = $i3 unless $seq;
			}
		}
		$geno{$seqID} = $seq;
		close FILE;
			
		my $seqnum=0;
		open(OUT, ">$TEMdir/TEM3_GeneSeqs/$fgenome") or die "Step 4: Can't open '$TEMdir/TEM3_GeneSeqs/$fgenome': $!\n";
		
		my $hh = "$hh[1]__$hh[0]";
		open (FIL, "<$TEMdir/TEM2_HMMsearch_outS/$hh") or die "Step 4: Can't open '$TEMdir/TEM2_HMMsearch_outS/$hh': $!\n";
		my $seqTag="$hh[0]";		
		$seqTag=~ s/\.fas/_/;		
		foreach my $i2 (<FIL>) {
			$i2=~ s/\n// unless $i2=~s/\r\n//;
			next unless $i2;
			next if $i2=~ /\#/;
			my @in0 = split /\s+/, $i2;
			my $keyr;
			if ($in0[3] =~ /-/) {
				$keyr =$in0[2];
			} else {
				$keyr =$in0[3];
			}	
			my $geneSeq = $geno{$in0[0]};			
			print OUT ">$seqTag"."$keyr\n$geneSeq\n";
			$seqnum++;
		}
		close FIL;	
			
		close OUT;
		print "	$numg: Got $seqnum/$fgn sequences from genome: $fgenome\n";
	
	}
	print "=====\nTotally, $numg sequence files has been written to '$TEMdir/TEM3_GeneSeqs'!=====\n";
	print LOG "=====\nTotally, $numg sequence files has been written to '$TEMdir/TEM3_GeneSeqs'!=====\n";


#============= Step 4: Gather the Orthologs in Single File ==============# 
	print "\n\n\n#============= Step 4: Gather the Orthologs in Single File ==============# \n\n";
	print LOG "\n\n\n#============= Step 4: Gather the Orthologs in Single File ==============# \n\n";

	opendir(DIR1, "$TEMdir/TEM3_GeneSeqs") || die "Step 4: Can't open input directory '$TEMdir/TEM3_GeneSeqs': $!\n";
	my @in1 = readdir(DIR1);
	closedir(DIR1);
	die "Step 4: Input directory '$TEMdir/TEM3_GeneSeqs' does not contain any files.\n" unless scalar(@in1);

	if (-e "$TEMdir/TEM4_GeneCluster") {################Different from Windows version
		my $del2="rm -rf $TEMdir/TEM4_GeneCluster";  ################Different from Windows version
		system ($del2);
		die "Step 4: Directory '$TEMdir/TEM4_GeneCluster' already exist, but permission denied when delete it.\n" if (-e "$TEMdir/TEM4_GeneCluster"); 
	} 
	
	
	mkdir "$TEMdir/TEM4_GeneCluster";
	my %outClu;
	my $num5=0;
	open (TRE, ">$TEMdir/intreeInfo.txt") or die "Step 4: Can't open '$TEMdir/intreeInfo.txt': $!\n";
	my $newID=10000;
	foreach my $in5 (@fgenome) {	
		my $tt =$in5;
		$tt =~ s/\:\d+$//;
		open (FILE, "<$TEMdir/TEM3_GeneSeqs/$tt") or die "Step 4: Can't open '$TEMdir/TEM3_GeneSeqs/$tt': $!\n";
		$tt =~ s/\.fas$//;
		my ($id,$sig);
		$newID++;
		print TRE "SPE$newID	$tt\n";
		foreach my $i1 (<FILE>) {
			$i1=~ s/\n// unless $i1=~s/\r\n//;
			next unless $i1;
			if ( $i1 =~ /\>/ ) {
				my @seqids=split /_/, $i1;
				$id =$seqids[1];
				$sig=1 if $id ~~ @fgene; 
			
			} elsif ($sig) {
				unless ($outClu{$id}) {
					$outClu{$id}= ">$tt\n$i1\n";
				} else {
					$outClu{$id} .= ">$tt\n$i1\n"; 
				}
				$sig=0;
			}
		}
		close FILE;
		$num5++;
		my $ggggg=$#fgenome+1;
		my $rat =&givePercentage($num5, $ggggg);
		print ("\r", "	Step 4: Reading the gene set ($rat) '$TEMdir/TEM3_GeneSeqs/$tt'. 		");

	}
	close TRE;
	my @rrr =keys %outClu;
	print "\n\n";
	my $n=0;
	my $cc=$#rrr+1;	
	foreach my $in (keys %outClu) {
		my $cnt;
		$n++;
		open (OUT, ">./$TEMdir/TEM4_GeneCluster/$in.fas") or die "Step 4: Can't open '$TEMdir/TEM4_GeneCluster/$in.fas': $!\n";
		print OUT "$outClu{$in}";
		$outClu{$in} =~ s/>/"X".$cnt++/ge;
		my $rat= &givePercentage($n, $cc);
		print ("\r", "	Step 4: Writing each gene clusters in single file. ($rat)		");
		print LOG "	Step 4: Writing $cnt sequences of gene cluster '$in' in file './$TEMdir/TEM4_GeneCluster/$in.fas'\n";
		close OUT;
	}
	print ".\n";	
	print LOG "In total, $n gene clusters were written to '$TEMdir/TEM4_GeneCluster'.\n";

#============= Step 5: Do Alignment ==============# 
	print "\n\n\n#============= Step 5: Do Alignment ==============# \n\n"; 
	print LOG "\n\n\n#============= Step 5: Do Alignment ==============# \n\n";
	
	if (-e "$TEMdir/TEM5_Alignment") {
		my $del2="rm -rf $TEMdir/TEM5_Alignment";  #############Different from Windows version
		system ($del2);
		die "Step 5: Directory '$TEMdir/TEM5_Alignment' already exist, but permission denied when delete it.\n" if (-e "$TEMdir/TEM5_Alignment"); 
	} 
	
	mkdir "$TEMdir/TEM5_Alignment";
	my $num7 = 0;

	foreach (@fgene) {
		my $program = "./bin/clustalo";############Different from Windows version
		my $outname = $_ . ".fas.fasta";
		my @para =("-i", "./$TEMdir/TEM4_GeneCluster/$_.fas", "-o", "./$TEMdir/TEM5_Alignment/$outname", "--outfmt=fasta", "--output-order=tree-order","--threads=$thread","-v","--force");###Different from Windows version
		
		unless (system "$program", @para) {
			$num7++;
			print "\nStep 5: $_ has been completed! It is the $num7/$fgn file.\n\n";
		} else {
			die "Step 5: Can't execute './bin/clustalo': $!.\n";   ("-i", "./$TEMdir/TEM6_DedupGeneCluster/$_.fas", "-o", "./$TEMdir/TEM7_Alignment/$outname", "--outfmt=fasta", "--output-order=tree-order","--threads=$thread","-v","--force");###Different from Windows version
		}
	}
	print LOG "\nTotally, $num7 alignments has been created in '$TEMdir/TEM5_Alignment'!\n\n";

#============= Step 6: Trim Algnments ==============# 
	print "\n\n\n#============= Step 6: Trim Algnments ==============#\n\n";
	print LOG "\n\n\n#============= Step 6: Trim Algnments ==============#\n\n";
	
	if (-e "$TEMdir/TEM6_AlnTrimmed") {####################
		my $del2="rm -rf $TEMdir/TEM6_AlnTrimmed";  #####Different from Windows version
		system ($del2);
		die "Step 6: Directory '$TEMdir/TEM6_AlnTrimmed' already exist, but permission denied when delete it.\n" if (-e "$TEMdir/TEM6_AlnTrimmed"); 
	} 
	
	mkdir "$TEMdir/TEM6_AlnTrimmed";

	opendir(DIR, "$TEMdir/TEM5_Alignment") || die "Step 6: Can't open input directory '$TEMdir/TEM5_Alignment':$!\n";
	my @in8 = readdir(DIR);
	closedir(DIR);
	die "Step 6: Input directory $TEMdir/TEM5_Alignment does not contain any files.\n" unless scalar(@in8);

	my (@problem,@error);
	my $num8 = 0;

	my $uuu=$#fgene+1;
	foreach my $filex (@fgene) {
		my $file = $filex;
		$file .= "\.fas\.fasta";
		(my $outname = $file) =~ s/\.fas//;

		my $cmd76 = "./bin/trimal -in $TEMdir/TEM5_Alignment/$file -out $TEMdir/TEM6_AlnTrimmed/$outname -$trim";####Different from Windows version
		system ($cmd76);
		$num8++;
		my $rat =&givePercentage($num8, $uuu);
		print ("\r", "	Step 6: Trimming gene cluster ($rat): $filex.");

	}
	print LOG "\nTotally, $num8 alignments has been trimmed successfully and written to '$TEMdir/TEM6_AlnTrimmed'!\n\n";
	print "\nTotally, $num8 alignments has been trimmed successfully and written to '$TEMdir/TEM6_AlnTrimmed'!\n\n";
	
}


############### Task 4: Tree Inference #################  
my ($fast,$iq);
if ($task eq "all" || $task eq "tree_infer") {
	open (OPT, "<tree_app-options.txt") or die "Can't open 'tree_app-options.txt': $!\n";
	
	foreach my $i3 (<OPT>) {
		$i3=~ s/\n// unless $i3=~s/\r\n//;
			next unless $i3;
		if ( $i3 =~ s/IQ-TREE Command-line=// ) {
				$iq= $i3;
		} elsif ($i3 =~ s/FastTree Command-line=//) {
			$fast= $i3;
		}
	}
	close OPT;
	if ($tree_type eq "sm") {
#============= Step 7: Concatenate Sequences (If Necessary) ==============# 

		print "\n\n############### Task 4: Tree Inference #################\n#============= Step 7: Concatenate Sequences ==============# \n\n";
		print LOG "\n\n############### Task 4: Tree Inference #################\n#============= Step 7: Concatenate Sequences ==============# \n\n";

		opendir(DIR, "$TEMdir/TEM6_AlnTrimmed") || die "Step 7: Can't open input directory '$TEMdir/TEM6_AlnTrimmed': $!\n";
		my @in10 = readdir(DIR);
		closedir(DIR);
		splice @in10, 1, 1 unless $in10[1]=~ /\.fasta/;
		splice @in10, 0, 1 unless $in10[0]=~ /\.fasta/;
		die "Step 7: Input directory $TEMdir/TEM6_AlnTrimmed does not contain any files. Maybe you need to run the former task(s).\n" unless $in10[2];

		my $iii=$#in10+1;
		my $num10 = 0;
		my (%out,@id,@seq,$id,$seq);
		my $gap0;   ###############for gaps "-";
		foreach my $file ( @in10) {
			next unless $file =~ /\.fasta$/;
			my $gap1; 
			unless (@id) {
				open (FILE, "$TEMdir/TEM6_AlnTrimmed/$file") or die "Step 7: Can't open '$TEMdir/TEM6_AlnTrimmed/$file': $!\n";
				foreach my $i1 (<FILE>) {
					$i1=~ s/\n// unless $i1=~s/\r\n//;
					next unless $i1;
					if ( $i1 =~ /\>/ ) {
						if ($seq) {
							push @seq, $seq; 
							$seq=undef;
						}
						my @tri =split /\s/,$i1;
						push @id, $tri[0];
					} else {
						$i1 =~ s/\*/-/g;
						$seq .="$i1";
					}
				}
				close FILE;
				push @seq,$seq;
				$gap1=$seq;
				$gap1=~ s/[A-Z]/-/g;
				if ($gap0) {
					$gap0 .="$gap1";
				} else {
					$gap0 ="$gap1";
				}
				foreach my $i2 (0..$#id) {
					$out{$id[$i2]}=$seq[$i2];
				}
			} else {
				open (FILE, "<$TEMdir/TEM6_AlnTrimmed/$file") or die "Step 7: Can't open '$TEMdir/TEM6_AlnTrimmed/$file': $!\n";
				my ($gapid,@keys7);
				foreach my $i1 (<FILE>) {
					$i1=~ s/\n// unless $i1=~s/\r\n//;
					next unless $i1;
					if ( $i1 =~ /\>/ ) {
						my @tri =split /\s/,$i1;
						
						$id= $tri[0];
						$gapid =$tri[0] unless $gapid;
						push @keys7,$id;
						unless ($id ~~ @id) {
							push @id, $id;
							$out{$id} ="$gap0";
						}	
					} else {
						$i1 =~ s/\*/-/g;
						$out{$id} .="$i1";
						if ($id eq $gapid) {
							$gap1 .="$i1" if $gap1;
							$gap1 ="$i1" unless $gap1;
						}
					}
				}
				close FILE;
				$gap1=~ s/[A-Z]/-/g;
				foreach my $vv (@id) {
					$out{$vv} .="$gap1" unless $vv ~~ @keys7;
				}
			}
			$num10++;
			my $rat = &givePercentage($num10, $iii);
			print ("\r", "	Step 7: Reading the files: $num10/$iii ($rat)  .");

		}
		print ".\n\n";

		open (OUT, ">$inputDir.concatenation.fas") or die "Can't open '$inputDir.concatenation.fas': $!\n";
		my $nn=0;
		foreach (sort keys %out) {
			$nn++;
			print OUT "$_\n";
			my $i1=$out{$_};
			while (1) {
				if (length($i1)>60) {
					my $sub = substr($i1,0,60);
					print OUT "$sub\n";	
					$i1 =~ s/$sub//;	
				} else {
					print OUT "$i1\n";
					last;
				}
			}
		}
		close OUT;
		print "Got $nn concatenated sequences written to '$inputDir.concatenation.fas'!\n";
		print LOG "Got $nn concatenated sequences written to '$inputDir.concatenation.fas'!\n";

#============= Step 8: Phylogeny Inference ==============# 
		print "\n\n\n#============= Step 8: Phylogeny Inference ==============#\n\n";
		print LOG "\n\n\n#============= Step 8: Phylogeny Inference ==============#\n\n";

		
		print "Infering phylogeny of the supermatrices approach......\n\n";
		
		
		
		my $cmd11;
		if ($tree_app eq "fasttree") {
			$cmd11="./bin/$fast $inputDir.concatenation.fas > $inputDir.supermatrix.$tree_app.tree";	##Different from Windows version
			&FastreeResponse ($cmd11);
			print LOG "The supermatrix tree has been written successfully in '$inputDir.supermatrix.$tree_app.tree'.\n\n";
			print "The supermatrix tree has been written successfully in '$inputDir.supermatrix.$tree_app.tree'.\n\n";
		} elsif ($tree_app eq "iqtree") {
			$cmd11="./bin/$iq -s $inputDir.concatenation.fas -T $thread";##Different from Windows version
			system ($cmd11); 
			open (FIL77, "<$inputDir.concatenation.fas.contree") or die "Step 8: Can't open '$inputDir.concatenation.fas.contree': $!\n";
			open (FIL88, ">$inputDir.supermatrix.$tree_app.tree") or die "Step 8: Can't open '$inputDir.supermatrix.$tree_app.tree': $!\n";
			foreach my $i0 (<FIL77>) {
				print FIL88 "$i0";
						
			}
			close FIL77;
			close FIL88;
			print LOG "The supermatrix tree has been written successfully in '$inputDir.supermatrix.$tree_app.tree'.\n\n";
			
			print "The supermatrix tree has been written successfully in '$inputDir.supermatrix.$tree_app.tree'.\n\n";
		}

	} elsif ($tree_type eq "st") {
	#============= Step 7: Construct gene Trees ==============#
		opendir(DIR, "$TEMdir/TEM6_AlnTrimmed") || die "Step 7: Can't open input directory '$TEMdir/TEM6_AlnTrimmed': $!\n";
		my @in11 = readdir(DIR);
		closedir(DIR);
		die "Step 7: Input directory $TEMdir/TEM6_AlnTrimmed does not contain any files.\n" unless scalar(@in11);
		splice @in11, 1, 1 unless $in11[1]=~ /\.fasta/;
		splice @in11, 0, 1 unless $in11[0]=~ /\.fasta/;
		my $seqnum=1;
		my $report="**Sequence numbers in some files:\n";
		foreach my $file11 ( @in11) {
			next unless $file11 =~ /\.fasta$/;
			my $nnum=&countSeq("$TEMdir/TEM6_AlnTrimmed/$file11");
			$report .="$TEMdir/TEM6_AlnTrimmed/$file11	$nnum\n";
			if ($seqnum==1) {
				$seqnum=$nnum;
			} elsif ($seqnum>$nnum || $seqnum<$nnum) {
				die "ERROR:\nTask 4 (-task tree_infer): Options '-tree st' requires the alignment files in directory '$TEMdir/TEM6_AlnTrimmed/' contains the same number of sequences.\n$report\nPlease run EasyCGTree from 'Task 2 (-task refine)' or run the whole pepline.\n\n****Tips: if this ERROR still occures while you have run 'Task 2' or the whole pepline, the reason must be that some proteomes have no sequence retained during trimming of some gene clusters by trimAl. Please remove related cluster(s) from '$TEMdir/TEM6_AlnTrimmed/' manually and run 'Task 4 (-task tree_infer)' again.\n\n";
			}
		}
		print "\n\n############### Task 4: Tree Inference #################\n#============= Step 7: Construct gene Trees ==============#\n\n";
		print LOG "\n\n############### Task 4: Tree Inference #################\n#============= Step 7: Construct gene Trees ==============#\n\n";
		if (-e "$TEMdir/TEM7_GeneTrees") {####################
			my $del2="rm -rf $TEMdir/TEM7_GeneTrees";  ######Different from Windows version
			system ($del2);
			die "Step 7: Directory '$TEMdir/TEM7_GeneTrees' already exist, but permission denied when delete it.\n" if (-e "$TEMdir/TEM7_GeneTrees"); 
		} 
		
		mkdir "$TEMdir/TEM7_GeneTrees";

		my $num10=$#in11+1;
		my $num11 = 0;
		foreach my $file11 ( @in11) {
			next unless $file11 =~ /\.fasta$/;
			$num11++;
			
			print "Step 7: Infering phylogeny from '$TEMdir/TEM6_AlnTrimmed/$file11': $num11/$num10.\n";
			my $cmd11;
			if ($tree_app eq "fasttree") {
				$cmd11="./bin/$fast $TEMdir/TEM6_AlnTrimmed/$file11 > $TEMdir/TEM7_GeneTrees/$file11.tree";	##Different from Windows version
				&FastreeResponse ($cmd11);
			} elsif ($tree_app eq "iqtree") {
				$cmd11="./bin/$iq -s $TEMdir/TEM6_AlnTrimmed/$file11 -T $thread";##Different from Windows version
				system ($cmd11); 				
				
				##Different from Windows version
				my $mv1="mv -f $TEMdir/TEM6_AlnTrimmed/$file11.ckp.gz $TEMdir/TEM7_GeneTrees/$file11.ckp.gz";####
				system ($mv1);
				my $mv2="mv -f $TEMdir/TEM6_AlnTrimmed/$file11.model.gz $TEMdir/TEM7_GeneTrees/$file11.model.gz";####
				system ($mv2);
				my $mv3="mv -f $TEMdir/TEM6_AlnTrimmed/$file11.bionj $TEMdir/TEM7_GeneTrees/$file11.bionj";####
				system ($mv3);
				my $mv4="mv -f $TEMdir/TEM6_AlnTrimmed/$file11.contree $TEMdir/TEM7_GeneTrees/$file11.contree";####
				system ($mv4);
				my $mv5="mv -f $TEMdir/TEM6_AlnTrimmed/$file11.iqtree $TEMdir/TEM7_GeneTrees/$file11.iqtree";####
				system ($mv5);
				my $mv6="mv -f $TEMdir/TEM6_AlnTrimmed/$file11.mldist $TEMdir/TEM7_GeneTrees/$file11.mldist";####
				system ($mv6);
				my $mv7="mv -f $TEMdir/TEM6_AlnTrimmed/$file11.splits.nex $TEMdir/TEM7_GeneTrees/$file11.splits.nex";####
				system ($mv7);
				my $mv8="mv -f $TEMdir/TEM6_AlnTrimmed/$file11.treefile $TEMdir/TEM7_GeneTrees/$file11.iqtreefile";####
				system ($mv8);
				my $mv9="mv -f $TEMdir/TEM6_AlnTrimmed/$file11.log $TEMdir/TEM7_GeneTrees/$file11.log";####
				system ($mv9);
				my $mv10="mv -f $TEMdir/TEM6_AlnTrimmed/$file11.uniqueseq.phy $TEMdir/TEM7_GeneTrees/$file11.uniqueseq.phy";####
				system ($mv10) if -e "$TEMdir/TEM6_AlnTrimmed/$file11.uniqueseq.phy";
				##Different from Windows version

			}
	
			
		}
		print LOG "Totally, $num11 phylogenies has been written to '$TEMdir/TEM7_GeneTrees'.\n\n";
		print  "Totally, $num11 phylogenies has been written to '$TEMdir/TEM7_GeneTrees'.\n\n";
			
	#============= Step 8: Generate consensus Tree ==============#		
		print "\n\n\n#============= Step 8: Generate consensus Tree ==============#\n\n";
		print LOG "\n\n\n#============= Step 8: Generate consensus Tree ==============#\n\n";
		open (TRE, "<$TEMdir/intreeInfo.txt") or die "Step 8: Can't open '$TEMdir/intreeInfo.txt': $!If it is absent, you should run 'task refine' again.\n";
		my %idddd;
		foreach my $in11 (<TRE>) {
			$in11=~ s/\n// unless $in11=~s/\r\n//;
			next unless $in11=~ /\d/;
			my @in7b=split /	/, $in11;
			my $keya=$in7b[1];
			
			$idddd{$keya}=$in7b[0];
		}
		close TRE;
		
		open (OU11, ">intree") or die "Step 8: Can't open 'intree': $!\n";
		foreach my $file11 (@in11) {
			$file11 .="\.tree" if $tree_app eq "fasttree";
			$file11 .="\.contree" if $tree_app eq "iqtree";
			open (FIL11, "<$TEMdir/TEM7_GeneTrees/$file11") or die "Step 8: Can't open '$TEMdir/TEM7_GeneTrees/$file11': $!\n";
			foreach my $in11 (<FIL11>) {
			$in11=~ s/\n// unless $in11=~s/\r\n//;
				next unless $in11=~ /\d/;
				foreach my $reid (keys %idddd) {
					my $vv=$idddd{$reid};
					$in11=~ s/$reid/$vv/; 
				}
				print OU11 "$in11\n";
			}
			close FIL11;
		}
		close OU11;
		my $cmd11t;
		$cmd11t="./bin/consenseM ";##Different from Linux version
		my $response11= `$cmd11t`;  
		print "$response11\n";
		open (OU11, "<outtree") or die "Step 8: Can't open 'outtree': $!\n";
		open (OU12, ">$inputDir.supertree.$tree_app.tree") or die "Step 8: Can't open '$inputDir.supertree.$tree_app.tree': $!\n";
		my $mytree;
		foreach my $in11 (<OU11>) {
			$in11=~ s/\n// unless $in11=~s/\r\n//;
			next unless $in11 =~ /\w/;
			$mytree .= $in11 if $mytree;
			$mytree = $in11 unless $mytree;
		}
		close OU11;
		$mytree=~ s/\)\:/\)/g; ############# replace '):'
		$mytree=~ s/\:\d+\.0+\,/\,/g; ############# replace ':10.00,'
		$mytree=~ s/\:\d+\.0+\)/\)/g; ############# replace ':10.00('
		$mytree =reverse $mytree;
		$mytree =~ s/EPS\,0+\.\d+\)/EPS\,\)/;
		$mytree =reverse $mytree;
		foreach my $reid (keys %idddd) {
			$mytree=~ s/$idddd{$reid}/$reid/; 
		}	
			
		print OU12 "$mytree\n";
		
		close OU12;
		print "The supertree has been written successfully in '$inputDir.supertree.$tree_app.tree'.\n\n";
		print LOG "The supertree has been written successfully in '$inputDir.supertree.$tree_app.tree'.\n\n";
	} 
}

if ($tree_app eq "iqtree") { 
	my $sig4;
	if ($iq =~ m/-m\s+MF/) {
		$sig4=1 unless $iq =~ m/-m\s+MFP/;
	} elsif ($iq =~ m/-m/) {
	} else {
		$sig4=1;
	}
	$sig4=0 unless $task eq "all" || $task eq "tree_infer";
	
	if ($sig4) {

	
		if ($tree_type eq "sm") {
			print "\n\n*******************************\n\nUsefull Information from IQ-TREE:\n\nBest-fit model	-LnL	BIC\n\n"; 
			print LOG "\n\n*******************************\n\nUsefull Information from IQ-TREE:\n\nBest-fit model	-LnL	BIC\n\n"; 
			open (INFO, "<$inputDir.concatenation.fas.log") or die "Step 8: Can't open 'outtree': $!\n";
			my $isig=0;
			my %modelinfo;
			my $best;
			foreach my $inf (<INFO>) {
				$inf=~ s/\n// unless $inf=~s/\r\n//;
				next unless $inf =~ /\w/;
				$isig=1 if $inf =~ /No. Model/;
				$isig=0 if $inf =~ /Corrected Akaike Information Criterion/;
				if ($isig) {
					$inf =~ s/^\s+//;
					my @modin = split /\s+/, $inf;
					my $key1111= $modin[1];
					$modelinfo{$key1111} ="$modin[2]	$modin[6]" if $modin[6];
				}
				if ($inf =~ /Best-fit model/) {
					$best = $inf;
					$best=~ s/Best-fit model: //; 
					$best=~ s/ chosen according to BIC//; 
					last;
				}
		
			}
			close INFO;
			print "$best	$modelinfo{$best}\n\n"; 
			print LOG "$best	$modelinfo{$best}\n\n"; 
	
	
		} elsif ($tree_type eq "st") {
			print "\n\n*******************************\n\nUsefull Information from IQ-TREE:\n\nNo.	Gene	Best-fit model	-LnL	BIC\n\n"; 
			print LOG "\n\n*******************************\n\nUsefull Information from IQ-TREE:\n\nNo.	Gene	Best-fit model	-LnL	BIC\n\n"; 
			my $num22=0;
			opendir(DIR, "$TEMdir/TEM7_GeneTrees") || die "Step 10: Can't open input directory '$TEMdir/TEM7_GeneTrees': $!\n";
			my @in11 = readdir(DIR);
			closedir(DIR);
			die "Input directory $TEMdir/TEM7_GeneTrees does not contain any files.\n" unless scalar(@in11);
			foreach my $file11 ( @in11) {
				next unless $file11 =~ /\.log$/;
				$num22++;
				my $genet =$file11;
				$genet =~ s/\.fasta\.log//;
				open (INFO, "<$TEMdir/TEM7_GeneTrees/$file11") or die "Can't open '$TEMdir/TEM7_GeneTrees/$file11': $!\n";
				my $isig=0;
				my %modelinfo;
				my $best;
				foreach my $inf (<INFO>) {
					$inf=~ s/\n// unless $inf=~s/\r\n//;
					next unless $inf =~ /\w/;
					$isig=1 if $inf =~ /No. Model/;
					$isig=0 if $inf =~ /Corrected Akaike Information Criterion/;
					if ($isig) {
						$inf =~ s/^\s+//;
						my @modin = split /\s+/, $inf;
						my $key1111= $modin[1];
						$modelinfo{$key1111}="$modin[2]	$modin[6]" if $modin[6];
					}
					if ($inf =~ /Best-fit model/) {
						$best = $inf;
						$best=~ s/Best-fit model: //; 
						$best=~ s/ chosen according to BIC//; 
						last;
					}
		
				}
				close INFO;
				print "$num22	$genet	$best	$modelinfo{$best}\n\n"; 
				print LOG "$num22	$genet	$best	$modelinfo{$best}\n\n"; 
					
			}

		}
		print "*******************************\n\n"; 
		print LOG "*******************************\n\n"; 
	}
}

############### Ending Information ###############
print "\n\n\n############### Ending Information ###############\n\n";
print LOG "\n\n\n############### Ending Information ###############\n\n";
my ($sec1,$min1,$hour1,$mday1,$mon1,$year1,$wday1,$yday1,$isdst1) = localtime(time);
my $rmon1=$mon1 +1;
my $ryear1=$year1+1900;


print LOG "Job was finished at: $hour1:$min1:$sec1,$ryear1-$rmon1-$mday1.\n";
print "Job was finished at: $hour1:$min1:$sec1,$ryear1-$rmon1-$mday1.\n";


my $sig=0;
my ($sec0,$min0,$hour0,$yday0);
if ($sec1>=$sec) {
	$sec0=$sec1-$sec;
	$sig=0;
} else {
	$sec0=$sec1-$sec+60;
	$sig=1;
}
	$min1=$min1-$sig;
if ($min1>=$min) {
	$min0=$min1-$min;
	$sig=0;
} else {
	$min0=$min1-$min+60;
	$sig=1;
}
	$hour1=$hour1-$sig;
if ($hour1>=$hour) {
	$hour0=$hour1-$hour;
	$sig=0;
} else {
	$hour0=$hour1-$hour+24;
	$sig=1;
}
	$yday1=$yday1-$sig;
if ($yday1>=$yday) {
	$yday0=$yday1-$yday;
	$sig=0;
} else {
	$yday0=$yday1-$yday+365;
	$sig=1;
}	
print "Running time: $yday0 d $hour0 h $min0 min $sec0 sec.\n\n\nEasyCGTree Version 4.0 by Dao-Feng Zhang\n\n"; ##########################注意更新版本和日期
print LOG "Running time: $yday0 d $hour0 h $min0 min $sec0 sec.\n\n\nEasyCGTree Version 4.0 by Dao-Feng Zhang\n\n";
close LOG;

sub givePercentage {
	my $g1 =shift;
	my $g2 =shift;	
	my $rar=$g1*100/$g2;
		
	$rar =int(100*$rar)/100;
	$rar .="0" if $rar =~ /\.\d$/;
	my $ret="$rar%";
	return $ret;
}



sub FastreeResponse {
	my $cmdsub =shift;
	my $response11= `$cmdsub`;  
	my @res11=split /\n/, $response11;
	print "$response11\n";
	foreach my $res1 (@res11) {
		if  ($res1 =~ /truncated/i && $res1 =~ /be too long/i) {
			print LOG "Step 7 ERROR: EasyCGTree stopped when execute 'FastTree':\n$response11\n\n\n\nPlease check whether the genome/proteome files were formated correctly. \nPlease find more informations in section 3.1 of the Manual.\n";##Different from Linux version
			die  "Step 7 ERROR: EasyCGTree stopped when execute 'FastTree':\n$response11\n\n\n\nPlease check whether the genome/proteome files were formated correctly. \nPlease find more informations in section 3.1 of the Manual.\n\n";##Different from Linux version
		} elsif ($res1 =~ /out of memory/i) {
			print LOG "Step 7 ERROR: EasyCGTree stopped when execute 'FastTree':\n$response11\n\n\n\nThis ERROR occured probably because 'TastTree' needs larger memory space than that of your computer. Please find a more powerful PC or server to complete the run.\n";##Different from Linux version
			die  "Step 7 ERROR: EasyCGTree stopped when execute 'FastTree':\n$response11\n\n\n\nThis ERROR occured probably because 'TastTree' needs larger memory space than that of your computer. Please find a more powerful PC or server to complete the run.\n";##Different from Linux version
		}
	}
}
sub countSeq {
    my $dna = shift;
    
    open (FILE, "<$dna") || die "Can't open '$dna':$!.\n";

	my $count=0;	
	foreach my $in (<FILE>) {
		$in=~ s/\n// unless $in=~s/\r\n//;
		if ($in =~ /\>/) {
			$count++;
		} 
	}
	return $count;				
}


sub reverse_complement {
        my $dna = shift;
        # reverse the DNA sequence
        my $revcomp = reverse($dna);
        $revcomp =~tr/tgca/TGCA/;
        # complement the reversed DNA sequence
        $revcomp =~tr/ACGT/TGCA/;
		return $revcomp;
}

sub rename {
	my $in = shift;
	opendir(DIR, $in) || die "Can't open input directory '$in'\n";
	my @files = readdir(DIR);
	closedir(DIR);
	die "Input directory $inputDir does not contain any files" unless scalar(@files);

	my $number1 = 0;
	foreach my $file (@files) {
		next unless $file =~ /\w/;
		my $file1=$file;
		$file1 =~ s/\s/-/g;
		$file1 =~ s/_/-/g;
		$file1=~ s/\.\w+\z/\.fas/;
		$file1 .=".fas" unless $file1=~ /\./;
		$file1=~ s/\(//;
		$file1=~ s/\)//;

		die "Can't rename './$inputDir/$file': $!" unless rename("./$inputDir/$file","./$inputDir/$file1"); 
		$number1++;

	}
	print "\nTotally, $number1 files have been renamed successfully!\n";
}


sub checkDataType {
	my $dna = shift;
	my $out;	
	unless (opendir(DIR, $dna)) {
		$out= "Can't open input directory '$dna': $!\n";
		return $out;
		last;
	}
	my @querysub = readdir(DIR);
	closedir(DIR);
	my @name;
	foreach my $in00 (@querysub) {
		next unless $in00 =~ /\w/;	
		unless ($in00 =~ /\.fas$/) {
			$out="it contain file(s) of which the name ends without '.fas'. Please formate them by using script 'FormatNames.pl'.\n";
			last;
		} else {
			push @name, $in00;
		}
	} 
	if ($out) {
		return $out;
	} else {
		unless (@name) {
			$out="Directory $dna does not contain any files.\n" ;
			return $out;
		}else {
			my @geneloc;
			while (1) {
				my $nn =int(rand($#name+1));
				push @geneloc, $nn unless $nn ~~ @geneloc;
				last if $#geneloc == $#name || $#geneloc == 12;
			}
			my @gene;	
			foreach my $file (@geneloc) {
				open (FILE, "./$dna/$name[$file]");
				my $type;
				my $count=0;	
				my $sig=0;
				my $seq;
				foreach my $in (<FILE>) {
					$in=~ s/\n// unless $in=~s/\r\n//;
					if ($in =~ /\>/) {
						$count++;
						$sig=1;
					} else {
						unless ($seq) {
						$seq="$in";
						} else {
							$seq .="$in";
						}
					}
					last if $count ==10;
				}
				close FILE;
				unless ($sig) {
					$out = "./$dna/$name[$file] seems not a fasta-formated file. Please check others.\n";
					last;
				} else {
					my $len = length($seq);
					my ($cntA,$cntG,$cntC,$cntT);
					$seq =~ s/A/"X".$cntA++/ge;
					$seq =~ s/G/"X".$cntG++/ge;
					$seq =~ s/C/"X".$cntC++/ge;
					$seq =~ s/T/"X".$cntT++/ge;					
					my $rio;
					if ($cntC) {
						$rio= ($cntA+$cntG+$cntC+$cntT)/$len;
					} else {
						$rio= ($cntA+$cntG+$cntT)/$len;
					}
					if ($rio>0.9) {
						$type = "nucl";	
					} else {
						$type ="prot";
					}
					$out =$type unless $out;
				}
				unless ($type eq $out) {
					$out = "it containing files of two types of sequence (both DNA and protein).\n" ;
					last; 
				}
			}
		}		
	}
	return $out;
}
