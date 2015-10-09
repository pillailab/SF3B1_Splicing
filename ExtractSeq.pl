# Program to extract sequence around splice sites of 
use warnings;

# Usage: perl ExtractSeq.pl
# INPUT:
# SJ file format= chr:start:end:gid:strand
$control_sj = 'control.junc.uniq'; 
$mutant_sj = 'mutant.junc.uniq'; 
$comm_sj = 'control-mutant.junc.uniq';

$gencode_sj = 'gencode.junc.uniq';
$genome = '/home/fas/pillai/akk29/scratch/data/Human/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa';

# Define the window size around ss
$upstream = 500; # sequnece length upstream of splice site
$downstream = 20; # sequence length downstream of splice site


&getGenome;

%ss_motif = ();
&ReadSJ($control_sj, 'wt');
&ReadSJ($mutant_sj, 'mut');
&ReadSJ($comm_sj, 'wt-mut');
&ReadSJ($gencode_sj, 'gencode');

sub ReadSJ {
	#------ read the splice junction file -----
	my ($sj_file, $sample) = @_;
	open SJFILE, "$sj_file" or die "cant open sj_file\n";
	<SJFILE>;
	open OUT, ">$sj_file.intron.fasta";
	print "Writing fasta file for intron sequence..\n";
	while ($sj_id = <SJFILE>) {
		chomp $sj_id;
		($chr, $sj_x, $sj_y, $strand, $shared_ss) = split /\:/, $sj_id;
		next if ($strand eq '.'); # ignore the SJ has no strand info 
		# get 5' and 3ss (0_based coord)
		if ($strand eq '+') {
			$ss5 = $sj_x; # sj_x refers to 1 nt upstream of 5' ss, so its already 0_based 5'ss   
			$ss3 = $sj_y - 1; # sj_y - 1 => to make 0_based coordininates
		} elsif($strand eq '-') {
			$ss3 = $sj_x;  
			$ss5 = $sj_y - 1; 
		}
		($motif5ss, $motif3ss, $intron) = &ExtractSeq ($ss5, $ss3, $strand, $chr);
		print OUT ">$sj_id\n$intron\n";
	
		# store and count each 5'/3' consensus dinucleotides (motifs)
		# count 5'ss and 3'ss dinucleotide motif if thier 5'ss, 3ss' or both are uniquwe to wt or mut 
		$ss_motif{$sample}{'motif5ss'}{$motif5ss} += 1 if ($shared_ss eq '01' or $shared_ss eq '00'); # 5' ss is unique  
		$ss_motif{$sample}{'motif3ss'}{$motif3ss} += 1 if ($shared_ss eq '10' or $shared_ss eq '00'); # 3' ss is unique

		# count 5'ss and 3'ss dinucleotide motif of SJ unique to wt and mut
		$motif53ss = $motif5ss.'/'.$motif3ss;
		$ss_motif{$sample}{'motif53ss'}{$motif53ss} += 1; 
	}
	close SJFILE;
	
	#$di_nt_type = scalar keys %{$ss_motif{$sammple}};
	#print "Total number of dinucleotides (or SJ) = $di_nt_count";
	
	print "For unique SJ:\n";
 	foreach $motif (sort keys %{$ss_motif{$sample}{'motif53ss'}}) {
		print "$sample:$motif:",$ss_motif{$sample}{'motif53ss'}{$motif}, "\n";
	}
	print "For unique 5SS:\n";
 	foreach $motif (sort keys %{$ss_motif{$sample}{'motif5ss'}}) {
		print "$sample:$motif:",$ss_motif{$sample}{'motif5ss'}{$motif}, "\n";
	}
	print "For unique 3SS:\n";
 	foreach $motif (sort keys %{$ss_motif{$sample}{'motif3ss'}}) {
		print "$sample:$motif:",$ss_motif{$sample}{'motif3ss'}{$motif}, "\n";
	}
	#print "For GENCODE unique SJ:\n";
 	#foreach $motif (sort keys %{$ss_motif{$sample}{'motif3ss'}}) {
	#	print "$sample:$motif:",$ss_motif{$sample}{'motif3ss'}{$motif}, "\n";
	#}
}
sub getGenome {
	%GENOME = ();
	open GENOME, "$genome";
	while ($line = <GENOME>) {
		chomp $line;
		if ($line =~ /^>/) {
			$line =~ s/>//g;
			$chr = $line;
			print "Preparing $chr sequence..\n";
		} else {
			$GENOME{$chr} .= $line;	
		}
		
	}
	close GENOME;
}

sub RevComp {
	# Reverse complement of sequence
	my $seq = $_[0]; 
	$revcomp_seq = reverse $seq;
	$revcomp_seq  =~ tr/ATGCatgc/TACGtacg/;
	return ($revcomp_seq);
}
sub ExtractSeq {
	# Extract 5'/3'SS dinucleotides
	# Extract intron seq upstream of 3'SS 

	my ($ss5, $ss3, $strand, $chr) = @_;

	# Get 5'/3'SS dinucleotides
	if ($strand eq '+') {
		$nt5ss = uc substr ($GENOME{$chr}, $ss5, 2);
		$nt3ss = uc substr ($GENOME{$chr}, $ss3 - 1, 2);
		$intron = uc substr ($GENOME{$chr}, $ss3 - 1 - $upstream, $upstream);
	} elsif ($strand eq '-') {
		$nt5ss = uc &RevComp (substr ($GENOME{$chr}, $ss5 - 1, 2));
		$nt3ss = uc &RevComp (substr ($GENOME{$chr}, $ss3, 2));
		$intron = uc &RevComp (substr ($GENOME{$chr}, $ss3, $upstream));
		#$exon = substr ($GENOME{$chr}, $ss - $downstream, $downstream);
	}
	return ($nt5ss, $nt3ss, $intron);
}
