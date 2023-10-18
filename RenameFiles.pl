  #!/usr/bin/perl  
use warnings;
use strict;
my @usage=qq(
====== EasyCGTree ======
     Version 2.0

Usage: perl formatGenomes.pl Folder style
styles:
1. GCA_000423565.1_ASM42356v1_genomic...  TO GCA_000423565.1.fas
2. GCA_000423565.1_ASM42356v1_genomic...  TO GCA_000423565.fas
3. XXX.asd...                             TO XXX.hgt
);

die join("\n",@usage)."\n" unless $ARGV[1];      
my $inputDir = $ARGV[0];
opendir(DIR, $inputDir) || die "Can't open input directory '$inputDir'\n";
my @files = readdir(DIR);
closedir(DIR);
die "Input directory $inputDir does not contain any files" unless scalar(@files);

my $number1 = 0;

if ($ARGV[1] ==1) {
	foreach my $file (@files) {
	next unless $file =~ /\w/;
	my $file1=$file;
	my @in=split /_/, $file1;
	my $outname = "$in[0]\_$in[1].fas";
	die "Can't rename './$inputDir/$file': $!" unless rename("./$inputDir/$file","./$inputDir/$outname"); 
	$number1++;

	}
} elsif ($ARGV[1] ==2) {

	foreach my $file (@files) {
	next unless $file =~ /\w/;
	my $file1=$file;
	my @in=split /\./, $file1;


	die "Can't rename './$inputDir/$file': $!" unless rename("./$inputDir/$file","./$inputDir/$in[0].fas"); 
	$number1++;

	}
} elsif ($ARGV[1] ==3) {

	foreach my $file (@files) {
	next unless $file =~ /\w/;
	my $file1=$file;
	$file1 =~ s/\.nuc/\.function/;


	die "Can't rename './$inputDir/$file': $!" unless rename("./$inputDir/$file","./$inputDir/$file1"); 
	$number1++;

	}
}
print "\nTotally, $number1 files have been renamed successfully!\n";