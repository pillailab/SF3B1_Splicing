# PERL script to extract splice junction from GENCODE GTF

# usage: perl gencode_sj.pl genocode.gtf 
use warnings;
$gencode_gtf = $ARGV[0]; 
chomp $gencode_gtf;

open GTF, $gencode_gtf or die "cant open $gencode_gtf!\n";
print "Getting exon coordinates\n";
open GENE, ">gencode.gene";

%EXONS = ();
while ($line = <GTF>) {
	chomp $line;
	next if ($line =~ /^#/);
	@fields = split /\t/,$line;
	$feature = $fields[2];
	#print $fields[8],"\n\n\n";<>;

	if ($feature eq 'gene') {
		# print gene boundary
		$chr = $fields[0];
		$strand = $fields[6];
		$i = index ($fields[8], 'gene_id');
		$j = index ($fields[8], ';', $i + 1);
		$gid = substr ($fields[8], $i + 9, $j - $i - 10);
		print GENE "$fields[0]\t$gid\t$fields[3]\t$fields[4]\t$strand\n"; 
	} elsif ($feature eq 'transcript') {
		$i = index ($fields[8], 'transcript_id');
		$j = index ($fields[8], ';', $i + 1);
		$tid = substr ($fields[8], $i + 15, $j - $i - 16);
		#print "$chr\t$gid\t$tid\t$strand\n";<>;
	} elsif ($feature eq 'exon') {
		$exon_beg = $fields[3];
		$exon_end = $fields[4];
		if ($strand eq '+') {
			push @{$EXONS{$gid}{$tid}{'exon'}},[$exon_beg,$exon_end];
			$EXONS{$gid}{$tid}{'strand'} = $strand;
			$EXONS{$gid}{$tid}{'chr'} = $chr;
		} elsif ($strand eq '-') {
			# push in reverse order, as coords go higher to lower for '-' strand 
			unshift @{$EXONS{$gid}{$tid}{'exon'}},[$exon_beg,$exon_end];
			$EXONS{$gid}{$tid}{'strand'} = $strand;
			$EXONS{$gid}{$tid}{'chr'} = $chr;
		}	 
	}	
}
close GTF;
close GENE;
print "gencode.gene has been generated in current directory\n";

&GetSJ;
sub GetSJ {
	print "Getting splice junctions\n";
	open OUT,  ">gencode.junc.uniq";
	print OUT "Gencode unique SJ\n";
	%seen = ();
	foreach my $gid (keys %EXONS) {
		foreach my $tid (keys %{$EXONS{$gid}}) {
			$exon_count = scalar @{$EXONS{$gid}{$tid}{'exon'}};
			next if ($exon_count == 1); # single-exon transcript
			foreach $i (0..$exon_count - 2) {
				$left = $EXONS{$gid}{$tid}{'exon'}->[$i][1]; # last nt of exon
				$right = ($EXONS{$gid}{$tid}{'exon'}->[$i+1][0]) - 1; # last nt of intron
				$strand = $EXONS{$gid}{$tid}{'strand'};
				$chr = $EXONS{$gid}{$tid}{'chr'};
				#print "$chr\t$gid\t$tid\t$left\t$right\$strand\n";

				# print unique splice junction ----
				$junction = $left.':'.$right;
				# To make SJ format as same as wt and mut, print '11'(means 5' & 3' ss can belong to both wt and mut))
				print OUT "$chr:$left:$right:$strand:11\n" unless ($seen{$gid}{$junction});
				$seen{$gid}{$junction} = 1;
			}
		}
	}
	close OUT;
	print "gencode.junc.uniq has been generated in current location\n";
}

