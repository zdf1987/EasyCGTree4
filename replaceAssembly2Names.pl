#!/usr/bin/perl
use warnings;
use strict; 
die "Usage: perl XXX.pl namelist tree > outTree\nExample: Cerasicoccus arenae	KCTC 12870	GCA_014651635.1\n\n" unless $ARGV[1];
open (PARA, "<$ARGV[0]") or die "Can't open '$ARGV[0]': $!\n";
my %name;
foreach my $in3 (<PARA>) {
	$in3=~ s/\n// unless $in3=~s/\r\n//;
	next unless $in3;
	my @in = split /\t/, $in3;
	my $key = $in[2];
	$key=~ s/_/-/; #
	my $value = $in[0];
	$value .=" $in[1]\/$in[2]";##/
	$name{$key}=$value;
}
close PARA;

open (PAR, "<$ARGV[1]") or die "Can't open '$ARGV[1]': $!\n";

foreach my $in3 (<PAR>) {
	$in3=~ s/\n// unless $in3=~s/\r\n//;
	next unless $in3;
	foreach my $ee (keys %name) {
		my $rr= $ee;
		$in3 =~ s/$ee/$name{$ee}/;####print "$name{$ee}\n" if $ee 
		$rr =~ s/GCA/GCF/;
		$in3 =~ s/$rr/$name{$ee}/;###print "$name{$ee}\n" if $rr
		
	}
	print "$in3\n";
}
	
			
			
			
			
			
			


