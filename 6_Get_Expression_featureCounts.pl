use strict;
use warnings;

if (@ARGV != 2) {die "Please provide directory of featurecounts output and an output file name\n";}

my $dir = $ARGV[0];
my $outfile = $ARGV[1];

open(my $ofh, ">", $outfile) or die $!;

my %Gene2ID2FragCount = ();
my %Gene2Length = ();
my @IDs = ();

foreach my $file (glob("$dir/*.fragmentcounts")) {
	my $ID = "ERR";
#	if ($file =~ /([ATCG]{5,})A/) {
#	if ($file =~ /([^\/]*)_1/) {
#	if ($file =~ /(Plate\d_Cell\d+)/) {
#	if ($file =~ /(_1_\d+)/) {
#	if ($file =~ /(E\w+)\.sort/) {
#	if ($file =~ /_([^_]+_Cell\d\d)/) {
	if ($file =~ /([^\/]*?)\./) {
		$ID = $1;
	} else {
		die "$file does not match\n";
	}	
	push(@IDs,$ID);
	my $fail = 0;
	open(my $ifh, $file) or die $!;
	while (<$ifh>) {
		chomp;
		if ($_ =~ /^#/ || $_ =~ /^Geneid/) {next;} #skip header & comments
		my @record=split(/\t/);
		my $gene = $record[0]; $gene =~ s/\s+//g;
		if (!defined($record[6])) {
			$fail = 1; last;
		}
		$Gene2ID2FragCount{$gene}->{$ID} = $record[6]; 
		$Gene2Length{$gene} = $record[5];
	} close ($ifh);
	if ($fail) {
		print STDERR "Fail: $file\n";
		#system("rm $file");
		pop(@IDs);
	}
}

print $ofh "Length\t".join("\t",@IDs)."\n";
foreach my $gene (keys(%Gene2ID2FragCount)) {
	print $ofh "$gene";
	print $ofh "\t".$Gene2Length{$gene};
	foreach my $ID (@IDs) {
		my $count = "NA";
		if (exists($Gene2ID2FragCount{$gene}->{$ID})) {
			$count = $Gene2ID2FragCount{$gene}->{$ID};
		} else { 
			$count = "0";
		}
		print $ofh "\t".$count;
	}
	print $ofh "\n";
}

my %ID2Unassigned = ();
foreach my $file (glob("$dir/*.fragmentcounts.summary")) {
	my $ID = "ERR";
	#if ($file =~ /_([^_]+_Cell\d\d)/) {
	#if ($file =~ /([ATCG]{5,})A/) {
	#if ($file =~ /([^\/]*)_1/) {
	#if ($file =~ /(Plate\d_Cell\d+)/) {
	#if ($file =~ /(E\w+)\.sort/) {
	#if ($file =~ /(_1_\d+)/) {
	if ($file =~ /([^\/]*?)\./) {
		$ID = $1;
	} else {
		die "$file does not match\n";
	}	
	my $fail = 0;
	open(my $ifh, $file) or die $!;
	<$ifh>; # header
	<$ifh>; #Assigned
	while (<$ifh>) {
		chomp;
		my @record=split(/\t/);
		if (!defined($record[1])) {$fail = 1;last;}
		$ID2Unassigned{$ID} += $record[1];
	} close ($ifh);
	if ($fail) {
		system("rm $file");
	}
}

print $ofh "__Unassigned_Various\t1";
foreach my $ID (@IDs) {
	print $ofh "\t".$ID2Unassigned{$ID};
} 
print $ofh "\n";
