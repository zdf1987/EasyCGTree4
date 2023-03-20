#!/usr/bin/perl
use warnings;
use strict; 
no warnings 'experimental::smartmatch';
use File::Copy qw(move mv);
use Getopt::Long;
use List::MoreUtils qw/uniq/;



my @usage=qq(
====== EasyCGTree ======
     Version 4.0 by Dao-Feng Zhang (Hohai University, China)
	 Update 2022-08-16

Usage: perl BuildHMM.pl [Options]
Options:
-proteome <String>	(Essential)
	Input data (proteomes) directory.
	
-hmm <String, 'bac120', 'rp1', et al.> (Optional)
	Profile HMM (in folder 'HMM') used for gene calling. [default: bac120]
	It will be ignored when a '~_TEM/TEM1_HMMsearch_out' folder is present. 
-thread <Int> (Optional)		
	Number of threads to be used by hmmsearch. [default: 4]
	It will be ignored when a '~_TEM/TEM1_HMMsearch_out' folder is present.
-evalue <Real, 10..0> (Optional)
	Expect value for screening hmmsearch hits. [default: 1e-10]
	
-help (Optional)				
	Display this message.

);



########## Default parameters ################

my $hmm ="bac120";
my $thread=4;
my $evalue=1e-10;

##########################################
my ($inputDir, $fn1,$check1,$check2,$check3,@query);
my %opt=qw();
GetOptions(\%opt,"proteome:s","hmm:s","thread:i","evalue:f","help!");

if (scalar(keys %opt )==0 || exists($opt{"help"})) {
	print join("\n",@usage)."\n\n";
	exit;
}

unless (exists($opt{"proteome"})) {
	print join("\n",@usage)."\n\nERROR: '-proteome' must be set to start.\n\n";
	exit;	
} else {
	$inputDir=$opt{"proteome"};
}

if (exists($opt{"hmm"})) {
	$hmm = $opt{"hmm"};
	my $ex="./HMM/$opt{'hmm'}.hmm";
	unless (-e $ex) {
		print join("\n",@usage)."\nERROR: Argument '-hmm $hmm': the profile HMM '$hmm' cannot be found in the directory './HMM'.\n\n";
		exit;
	}
} 


if (exists($opt{"thread"})) {	
	my $BNM=$opt{"thread"};
	if ($opt{"thread"} =~ /^\d+$/) {
		$thread=$opt{"thread"};
	} else {
		print join("\n",@usage)."\nERROR: Argument '-thread $BNM': '$BNM' is not an integer.\n\n";
		exit;
	}
}

$evalue=$opt{"evalue"} if exists($opt{"evalue"});



{
	unless ($inputDir=~/\w/ || $inputDir=~/\d/) {
		print join("\n",@usage)."\n\nERROR: Argument '-Gene_fam $inputDir': it should be a directory name, e.g., myGenome.\n\n";
		exit;
	}
	&rename($inputDir); 
	$check1= &checkDataType($inputDir);
	if ($check1 eq "nucl") {
		print join("\n",@usage)."\n\nERROR: directory '$inputDir' should contains file(s) of protein sequence, while some file(s) of nucleotide sequence were found.\n\n";
		exit;
	} elsif ($check1 eq "prot") {
	} elsif ($check1=~ /\.fas/ || $check1=~ /open input/) {
		print join("\n",@usage)."\n\nERROR: $check1\n\n";
		exit;
	} else {
		print join("\n",@usage)."\n\nERROR: directory '$inputDir' should only contains fasta-formated file(s) of protein sequence, while $check1\n\n";
		exit;		
	}
	opendir(DIR, $inputDir);
	@query = readdir(DIR);
	closedir(DIR);
	splice @query, 1, 1 unless $query[1]=~ /\.fas/;
	splice @query, 0, 1 unless $query[0]=~ /\.fas/;
	$fn1=$#query+1;
	unless ($fn1 >=1) {
		
		print join("\n",@usage)."\n\nERROR: Number of files (taxa on a tree) in Input directory '$inputDir' ($fn1) must be >=1 to start a run.\n\n";
		exit;
	}
}
my $TEMdir="$inputDir\_TEM";

unless ( -e "$TEMdir/TEM1_HMMsearch_out") {
	my $cmd="perl EasyCGTree.pl -proteome $inputDir -hmm $hmm -thread $thread -task hmmsearch";
	system ($cmd);
}

open(OUT, ">$hmm\_prevalence_$inputDir.txt") or die "Can't open '$hmm\_prevalence_$inputDir.txt': $!\n";
print OUT "Proteome	Pre_genes";
open(MMIN, "<$TEMdir/HMMinfo.txt") or die "Can't open '$TEMdir/HMMinfo.txt': $!\n";
my @genes;
my %score;
foreach my $hmmin (<MMIN>) {
	$hmmin=~ s/\n// unless $hmmin=~s/\r\n//;
	next unless $hmmin =~ /\w/;
	my @hmmi =split /\t/, $hmmin;
	push @genes, $hmmi[0];
	print OUT "	$hmmi[0]";
	$score{$hmmi[0]} = $hmmi[1];
}
close MMIN;
print OUT "\n";



opendir(DIR, "$TEMdir/TEM1_HMMsearch_out") || die "Can't open input directory '$TEMdir/TEM1_HMMsearch_out': $!\n";
my @in11 = readdir(DIR);
closedir(DIR);
die "Input directory '$TEMdir/TEM1_HMMsearch_out' does not contain any files.\n" unless scalar(@in11);		
splice @in11, 1, 1 unless $in11[1]=~ /\.fasta/;
splice @in11, 0, 1 unless $in11[0]=~ /\.fasta/;
		
		
foreach my $in (@in11) {
	next unless $in =~ /\.fas$/; 
	open (FIL11, "<$TEMdir/TEM1_HMMsearch_out/$in") or die "Can't open '$TEMdir/TEM1_HMMsearch_out/$in': $!\n";
	##
	
	my $geno = $in; 
	$geno =~ s/\.fas//; 
	my (%hmmSin,%hmmE,%hmmSinO);
	my $geneC=0;
	foreach my $in11 (<FIL11>) {
		$in11=~ s/\n// unless $in11=~s/\r\n//;
		next if $in11=~ /^\#/;
		next unless $in11=~ /\d/;
		my @in3 =split /\s+/, $in11;                           
			
		next if $in3[4] > $evalue;#####################
		
		my $sco = $score{$in3[2]}/4;      #####Drop it? If not, how to use it?
		next if $in3[5] < $sco;
		my $keyr=$in3[0];
		unless ($hmmSin{$keyr}) {
			$hmmSin{$keyr}=$in11;
			$hmmE{$keyr}=$in3[4];#print "$keyr\n";
		} else {
			$hmmSin{$keyr}=$in11 if $hmmE{$keyr}>$in3[4];
			$hmmE{$keyr}=$in3[4] if $hmmE{$keyr}>$in3[4];
		}
	}
		
	foreach my $ch (keys %hmmE) {
		my @in3 =split /\s+/, $hmmSin{$ch};                           
			
		my $keyr=$in3[2];
		unless ($hmmSinO{$keyr}) {
			
			$hmmSinO{$keyr}=$in3[0];##print "$hmmSin{$keyr}   $keyr\n";
			$geneC++;
		} else {
			$hmmSinO{$keyr} .="/$in3[0]";
		}
	}
	close FIL11;	
	print OUT "$geno	$geneC";
	foreach my $key (@genes) {
		if ($hmmSinO{$key}) {
			print OUT "	$hmmSinO{$key}";
		} else {
			print OUT "	-";
		}
	}
	print OUT "\n";
}

print "\nThe prevalence of $hmm gene set in proteomes of $inputDir has been wroten in '$hmm\_prevalence_$inputDir.txt'.\n\n";


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
	print "\nIn total, $number1 files have been renamed successfully!\n";
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
