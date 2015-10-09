# Program to define cryptic 3'ss usin SJ infomation 
# from sammple and gencode annotation
use warnings;
use Storable;
use Statistics::R;

# Usage: perl CrypticSS.pl
# INPUT:
# SJ file format= chr:start:end:gid:strand
$control_sj = 'control.junc.uniq'; 
$mutant_sj = 'mutant.junc.uniq'; 
$comm_sj = 'control-mutant.junc.uniq';
$gencode_sj = 'gencode.junc.uniq';

%SJCount = %{retrieve ('SJCount.bin')};

# !!!!! Dont forget to define the counts of replicates
$control_repl = 5; 
$mutant_repl = 8; 
#$genome = '/home/fas/pillai/akk29/scratch/data/Human/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa';
#&getGenome;

%SS5 = %SS3 = %TrackSJ = ();
&ReadSJ($control_sj, 'wt');
&ReadSJ($mutant_sj, 'mut');
&ReadSJ($comm_sj, 'wt-mut');
&ReadSJ($gencode_sj, 'gencode');

&Cryptic3SS('wt');
&Cryptic3SS('mut');
&Cryptic3SS('wt-mut');

sub ReadSJ {
	#------ read the splice junction file and store in hash -----
	my ($sj_file, $sample) = @_;
	open SJFILE, "$sj_file" or die "cant open sj_file\n";
	<SJFILE>;
	while ($sj_id = <SJFILE>) {
		chomp $sj_id;
		($chr, $sj_x, $sj_y, $strand, $shared_ss, $gname) = split /\:/, $sj_id;
		next if ($strand eq '.'); # ignore the SJ that has no strand info 
		#-----------------------------------------------------
		#----------  get 5'and 3'ss (0_based coord) ----------
		#-----------------------------------------------------
		# Adjust 1 nt position of SJ coords to get 5'ss and 3'ss coords y 
		if ($strand eq '+') {
			$ss5 = $sj_x; # sj_x refers to 1 nt upstream of 5' ss, so its already 0_based 5'ss   
			$ss3 = $sj_y - 1; # sj_y - 1 => to make 0_based coordininates
		} elsif($strand eq '-') {
			$ss3 = $sj_x;  
			$ss5 = $sj_y - 1; # same as above 
		}
		# Store all 3'ss for a 5'ss
		push @{$SS5{$sample}{$ss5}}, $ss3;
		# Also, store all 5'ss for a 3'ss
		push @{$SS3{$sample}{$ss3}}, $ss5;

		# to keep track to sj_id of corresponding 5'/3'ss 
		$TrackSJ{$sample}{"$ss5:$ss3"} = $sj_id;	
	}
	close SJFILE;
}

sub Cryptic3SS {
	# find nearest 3'ss in gencode downstream (preferably), otherwise upstream of 
	# sample 3'ss (cryptic)
	my $sample = $_[0];
	my ($ss5_sample, $ss3_sample, $ss3_gencode, $ss3_cano, $sj_id_sample, $sj_id_gencode, $constant_wt, $constant_mut);
	my $R = Statistics::R->new();
	open OUT, ">$sample.cryptic3ss";
	foreach $ss5_sample (sort {$a <=> $b} keys %{$SS5{$sample}}) {
		# to define cryptic 3'ss, the 5'ss of sample should be present
		# in gencode, skip otherwise 
		next unless (exists $SS5{'gencode'}{$ss5_sample});
		# iterate through all 3'ss for this 5'ss in sample
		foreach $ss3_sample (@{$SS5{$sample}{$ss5_sample}}) {
			# get sj id of corresponding 5' and 3' ss and the strand
			$sj_id_sample = $TrackSJ{$sample}{"$ss5_sample:$ss3_sample"}; 
			$strand = (split /\:/, $sj_id_sample)[-3];
			# find non-shared 3'ss for shared 5'ss of sample in gencode
			%dist3ss = ();
			$found_non_shared_3ss = 0;
			foreach $ss3_gencode (@{$SS5{'gencode'}{$ss5_sample}}) {
				# get distance b/w ss3_sample and ss3_gencode to determine non-shared 3'ss 
				$dist = $ss3_sample - $ss3_gencode;
				next if ($dist == 0); # skip 3'ss shared in gencode and sample
				# store abs. dist and the corresponding gencode 3'ss
				$dist3ss{'abs'}{abs $dist} = $ss3_gencode;
				# also store pos(+) and neg(-) dist to keep track on upstream/downstream 3'ss
				# of geneocode w.r.t. sample 3'ss
				$dist3ss{'pos'}{abs $dist} = $ss3_gencode if ($dist > 0); 
				$dist3ss{'neg'}{abs $dist} = $ss3_gencode if ($dist < 0);
				$found_non_shared_3ss = 1;
			}
			if ($found_non_shared_3ss) {
				#print "Strand:$strand";<>;
				# find nearest 3'ss for sample in gencode (in case multiple non-shared 3'ss)
				@tmp = sort {$a <=> $b} keys %{$dist3ss{'abs'}};	
				$shortest_dist = $tmp[0];
				# check if there are two identical shortest_dist; referring to two gencode 3'ss 
				# where one is upstream and another is downstream of sample 3'ss	
				if (exists $dist3ss{'pos'}{$shortest_dist} && exists $dist3ss{'neg'}{$shortest_dist}){
					# then prefer the gencode 3'ss downstream of sample 3'ss
					if ($strand eq '+') {
						$ss3_cano = $dist3ss{'neg'}{$shortest_dist};
					} elsif ($strand eq '-') {
						$ss3_cano = $dist3ss{'pos'}{$shortest_dist};
					} 
					$cryptic3ss_position = 'upstream';
				} elsif (!exists $dist3ss{'pos'}{$shortest_dist} && exists $dist3ss{'neg'}{$shortest_dist}) {
					$ss3_cano = $dist3ss{'abs'}{$shortest_dist};
					$cryptic3ss_position = 'upstream' if ($strand eq '+');	
					$cryptic3ss_position = 'downstream' if ($strand eq '-');	
				} elsif (exists $dist3ss{'pos'}{$shortest_dist} && !exists $dist3ss{'neg'}{$shortest_dist}) {
					$ss3_cano = $dist3ss{'abs'}{$shortest_dist};
					$cryptic3ss_position = 'downstream' if ($strand eq '+');	
					$cryptic3ss_position = 'upstream' if ($strand eq '-');	
				}	
				$sj_id_gencode = $TrackSJ{'gencode'}{"$ss5_sample:$ss3_cano"};
				#print "$sample  $sj_id_sample => $ss5_sample:$ss3_sample\n";
				#print "Gencode $sj_id_gencode => $ss3_cano\nDist = $shortest_dist ($cryptic3ss_position)\n";

				#--------------------------------- GET PSI --------------------------------
				if ($sample eq 'wt') {
					($psi_wt, $constant_wt) = &PSI ($sj_id_sample, $sj_id_gencode, $sample, $control_repl);
					# create list (referenced) of '0' for corresponding sample of comparison
					# Here, map function created array of 0, and the use of [] referenced that array
					$psi_mut = [map {0} 0..$mutant_repl - 1]; 
					$constant_mut = 1;
				} elsif ($sample eq 'mut') {
					($psi_mut, $constant_mut) = &PSI ($sj_id_sample, $sj_id_gencode, $sample, $mutant_repl);
					# create list (referenced) of '0' for corresponding sample of comparison
					$psi_wt = [map {0} 0..$control_repl - 1];
					$constant_wt = 1;
				} elsif ($sample eq 'wt-mut') {
					($psi_wt, $constant_wt) = &PSI ($sj_id_sample, $sj_id_gencode, 'wt', $control_repl);
					($psi_mut, $constant_mut)  = &PSI ($sj_id_sample, $sj_id_gencode, 'mut', $mutant_repl);
				}
				#print "WT($constant_wt):@$psi_wt:MUT($constant_mut)@$psi_mut\n";
				# get p-value and difference b/w av. psi(wt-mut)
				# using PERL-R interface
				if ($constant_wt && $constant_mut) {
					# do not perform t-test if  all elements (psi) are same 
					# in wt and also all psi's are same in mut. Regardless,
					# if the t-test is still perfomed, 'R' will throw stderr
					$pval = 1;
					$av_psi_wt = 'NA';
					$av_psi_mut = 'NA';
				} else {
					$R->set( 'x', $psi_wt);
					$R->set( 'y', $psi_mut);
					$R->run(qq `ttest <- t.test(x,y)`);
					$ttest = $R->get('ttest');
					$pval = $$ttest[16];
					$av_psi_wt = $$ttest[-2];
					$av_psi_mut = $$ttest[-1];
					$delta_psi = $av_psi_wt - $av_psi_mut;
				}	
				next if ($pval >= 0.05);
				#print "pval = $pval; DeltaPSI = $delta_psi\n";
				#print "\n";
				print OUT "$sample:$sj_id_sample\tGencode:$sj_id_gencode\t$shortest_dist\t$cryptic3ss_position\t$pval\t$delta_psi\n";
				#print "\n";
			}
				
			# if 3'ss of sample is present in gencode
			# it can't be a cryptic 3' ss, so skip this site
			#next if (exists $SS3{'gencode'}{$ss3});
			# otherwise, look for the nearest up or downstream 3'ss
		}
	} 
	close OUT;
	print "$sample.cryptic3ss has been created in the current directory\n";
}

sub PSI {
	my ($sj_id_sample, $sj_id_gencode, $sample, $repl_count) = @_;
	# rename sj_id for sample and gencode in order to 
	# access the read counts
	$sj_id_sample = join ":", (split (/\:/, $sj_id_sample))[0..2];
	$sj_id_gencode = join ":", (split (/\:/, $sj_id_gencode))[0..2];
	#print "$sj_id_sample, $sj_id_gencode, $sample\n";<>;

	# prepare read count array in the same order as sample replicates  
	my @sample_reads = my @ref_reads = ();
	foreach my $i (0..$repl_count - 1) {
		# read count array for sample junctions
		if (defined $SJCount{$sample}{$sj_id_sample}{'read_count'}->[$i]) {
			push @sample_reads, $SJCount{$sample}{$sj_id_sample}{'read_count'}->[$i];
		} else {
			push @sample_reads, 0;
		}
		# read count array for corresponding gencode junction
		if (defined $SJCount{$sample}{$sj_id_gencode}{'read_count'}->[$i]) {
			push @ref_reads, $SJCount{$sample}{$sj_id_gencode}{'read_count'}->[$i];
		} else {	
			push @ref_reads, 0;
		}	
	} 
	# calculate PSI
	my @psi_list = ();
	foreach $i (0..$repl_count - 1) {
		$count_a = $sample_reads[$i];
		$count_b = $ref_reads[$i];
		if ($count_a != 0 && $count_b != 0) {
			$psi = $count_a/($count_a + $count_b);
		} elsif ( ($count_a == 0 && $count_b == 0) || ($count_a == 0 && $count_b != 0) ) {
			$psi = 0;
		} elsif ( $count_a != 0 && $count_b == 0) {
			$psi = 1;
		}
		push @psi_list, $psi;
	} 
	#print "Sample Reads: @sample_reads\n";
	#print "Ref Reads: @ref_reads\n";
	#print "PSI: @psi_list\n";
	#print "Mean PSI: $mean_psi\n";	

	# Check if all elements of the array are indentical, as indentical values throws 
	# error during t-test calculation by R function
	my $constant = 1;	
	foreach my $value (@psi_list) {
		if ($value != $psi_list[0]) {
			$constant = 0;	
			last;
		}
	}
	return (\@psi_list, $constant);
}

=begin comment
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
=end comment
