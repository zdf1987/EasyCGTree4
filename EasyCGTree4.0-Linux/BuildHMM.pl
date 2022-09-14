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

Usage: perl EasyCGTree.pl inputDir
Aligned gene clusters: perl EasyCGTree.pl inputDir noAln
);

die join("\n",@usage) unless $ARGV[0];
my ($inputDir,$queryDir,$fn1,$check1,$check2,$check3,@query);


{
	$inputDir = $ARGV[0];
	unless ($inputDir=~/\w/ || $inputDir=~/\d/) {
		print join("\n",@usage)."\n\nERROR: Argument '$inputDir': it should be a directory name, e.g., myGenome.\n\n";
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




my $TEMdir= $inputDir."_BuildHMM_TEM";
if (-e $TEMdir) {
	my $del1="rm -rf $TEMdir";                ###########Different from Windows version
	system ($del1);
	die "Step 1: Directory '$TEMdir' already exist, but permission denied when delete it.\n" if (-e "$TEMdir"); 
}


mkdir "$TEMdir" || die " Permission denied to create directories. Please contact the Administrator of your System\n"; 
	

my $aln;
if ($ARGV[1] && $ARGV[1] eq "noAln") {
$aln="./$inputDir"
} else {
############### Do Alignment #################
mkdir "$TEMdir/TEM1_Alignment" || die " Permission denied to create directories. Please contact the Administrator of your System\n";
my $num7 = 0;

foreach (@query) {
	my $program = "./bin/clustalo";##Different from Windows version
	my @para = ("-i", "./$inputDir/$_", "-o", "./$TEMdir/TEM1_Alignment/$_", "--outfmt=fasta", "--output-order=tree-order","--threads=8","-v","--force");###Different from Windows version
	
	unless (system "$program", @para) {
		$num7++;
		print "\n********$_ has been completed! It is the $num7/$fn1 file.\n\n";
	} else {
		die "	Can't execute '.\\bin\\muscle.exe': $!.\n";   ##Different from Linux version
	}
}
$aln="./$TEMdir/TEM1_Alignment";
print "\n$num7 alignments has been created in '$TEMdir/TEM1_Alignment'!\n\n";
}

############### Build HMMs #################

mkdir "$TEMdir/TEM2_HMMs";

my $num8 = 0;

foreach my $filex (@query) {
	my $out=$filex;
	$out =~ s/\.fas//;
	my $cmd = "./bin/hmmbuild.exe --informat afa $TEMdir/TEM2_HMMs/$out.hmm $aln/$filex";##Different from Linux version
	
	system ($cmd); ##Different from Linux version
	$num8++;
	print "\n********$TEMdir/TEM1_Alignment/$filex has been processed! It is the $num8/$fn1 file.\n\n";
}
print "Totally, $num8 profile HMMs has been built successfully and written to '$TEMdir/TEM2_HMMs'!\n\n";


open (OUT, ">./HMM/$inputDir.hmm");
foreach my $filex (@query) {
	my $out=$filex;
	$out =~ s/\.fas//;
	open (FILE, "$TEMdir/TEM2_HMMs/$out.hmm");

	foreach my $in (<FILE>) {
		print OUT "$in";
	}
	close FILE;

}
close OUT;
print "The profile HMMs for the gene clusters in '$inputDir' has been built successfully and written to 'HMM/$inputDir.hmm'!\n\n";




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
