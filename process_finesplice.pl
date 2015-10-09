# This program reads 2 SJ files (chr:start:end)- WT and MUT, and outputs
# the corresponding files containing unique junctions with gid and strand

#Usage: perl process_finesplice.pl

use warnings;
use Storable;


# INPUT FILE:
# list of finesplice output files  (*.accepted.junc)
$control_list = 'finesplice_list.healthy';
$mutant_list = 'finesplice_list.mds';

# list of junction.bed files generated from TopHat alignment
$bedlist='bedlist'; 

# GTF file to link SJ to corresponding gid
#$gtf='/home/fas/pillai/akk29/scratch/data/Human/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf';

# OUTPUT FILES
# $control.uniq and $mutant.uniq

###########################################################################

#---- set read count cutoff ---
$wt_read_cutoff = 10;
$mut_read_cutoff = 10;
$com_read_cutoff = 10;


&ProcessTopHatSJ; # to get the strand of SJ
%SJ = %SS5 = %SS3 = ();	
&ProcessFineSplice ('wt', $control_list);
&ProcessFineSplice ('mut', $mutant_list);
store (\%SJ, "SJCount.bin"); print "SJCount.bin has been generated in current directory\n"; 
&GetUniqSJ;

sub ProcessFineSplice {
	# sum the SJ read counts over all replicates
	# format of sj_id = chr:start:end
	# store in hash: SJ{sample}{sj_id} = summed_read_count

	my ($sample, $finesplice_list) = @_;
	
	print "Reading FineSplice files for $sample...\n";
	open LIST, $finesplice_list or die "cant open $finesplice_list\n";
	$lib_index = 0;
	while ($filename = <LIST>) {
		chomp $filename;
		open FILE, $filename or die "cant open $filename\n";
		<FILE>;
		my @tmp = ();
		while ($sj_info = <FILE>) {
			chomp $sj_info;
			@tmp = split (/\t/, $sj_info);
			$sj_id = "$tmp[0]:$tmp[1]:$tmp[2]"; 
			$SJ{$sample}{$sj_id}{'summed_read_count'} += $tmp[4]; # read counts summed over replicates
			#push @{$SJ{$sample}{$sj_id}{'read_count'}}, $tmp[4]; # store read count in array 
			# store read count using array index corresponding to specific library 
			$SJ{$sample}{$sj_id}{'read_count'}->[$lib_index] = $tmp[4]; 
			#----------  get read strand ---------
			if (defined $TopHatSJ{$sj_id}) {
				$strand = $TopHatSJ{$sj_id}{'strand'};
			} else {
				$strand = '.';
			}
			# add strand information to SJ
			$SJ{$sample}{$sj_id}{'strand'} = $strand;

			#------ Store 5'& 3' SS information -----
			if ($strand eq '+') {
				$SS5{$sample}{$tmp[1]} = 1; 
				$SS3{$sample}{$tmp[2]} = 1; 
				$SJ{$sample}{$sj_id}{'5ss'} = $tmp[1]; 
				$SJ{$sample}{$sj_id}{'3ss'} = $tmp[2];
			} elsif ($strand eq '-') {
				$SS5{$sample}{$tmp[2]} = 1; 
				$SS3{$sample}{$tmp[1]} = 1; 
				$SJ{$sample}{$sj_id}{'5ss'} = $tmp[2];
				$SJ{$sample}{$sj_id}{'3ss'} = $tmp[1]; 
			} 
		}
		close FILE; 
		$lib_index++;
	}
	close LIST;
	$count = scalar keys %{$SJ{$sample}}; print "SJ in $sample = $count\n";
}

sub GetUniqSJ {
	print "Getting SJ uniq to WT, MUT and those common in both...\n";


	# find SJ uniq to WT, MUT and those common in both, sum the read counts
	# from replicates
	%SJ_GROUP = (); # wt, mut and wt-mut groups
	my $sj_id = "";

	#-------------------  iterarate through wt SJ -----------------
	foreach $sj_id ( keys %{$SJ{'wt'}} ) {
		#print $SJ{'wt'}{$sj_id}{'summed_read_count'};<>;
		if ($SJ{'mut'}{$sj_id})	{
			# if sj is common in wt and mut, the sum the supporting read counts from both
			$SJ_GROUP{'wt-mut'}{$sj_id}{'summed_read_count'} = $SJ{'wt'}{$sj_id}{'summed_read_count'} + $SJ{'mut'}{$sj_id}{'summed_read_count'};
		} else {
			# unique to wt
			$SJ_GROUP{'wt'}{$sj_id}{'summed_read_count'} = $SJ{'wt'}{$sj_id}{'summed_read_count'};
		}
	}
	#-------------------  iterarate through mut SJ -----------------
	foreach $sj_id (keys %{$SJ{'mut'}} ) {
		if ($SJ{'wt'}{$sj_id})	{
			# if sj is common in wt and mut, the sum the supporting read counts from both
			$SJ_GROUP{'wt-mut'}{$sj_id}{'summed_read_count'} = $SJ{'wt'}{$sj_id}{'summed_read_count'} + $SJ{'mut'}{$sj_id}{'summed_read_count'};
		} else {
			# unique to mut
			$SJ_GROUP{'mut'}{$sj_id}{'summed_read_count'} = $SJ{'mut'}{$sj_id}{'summed_read_count'};
		}
	}
	$count = scalar keys %{$SJ_GROUP{'wt'}}; print "SJ unique to WT (absent in MUT) = $count\n";
	$count = scalar keys %{$SJ_GROUP{'mut'}}; print "SJ unique to MUT (absent WT) = $count\n";
	$count = scalar keys %{$SJ_GROUP{'wt-mut'}}; print "SJ common in WT and MUT = $count\n";
	#-------------------------------------------------------------------
	# Filter SJ by read count cutoff
	# Check wt and mut share 5' and 3' junction site
	open WT_OUT, ">control.junc.uniq";
	open MUT_OUT, ">mutant.junc.uniq";
	open COM_OUT, ">control-mutant.junc.uniq";

	print WT_OUT "SJ with at least $wt_read_cutoff\n";
	print MUT_OUT "SJ with at least $mut_read_cutoff\n";
	print COM_OUT "SJ with at least $com_read_cutoff\n";

	# wt uniq 
	foreach $sj_id (keys %{$SJ_GROUP{'wt'}}) {
		next unless ($SJ_GROUP{'wt'}{$sj_id}{'summed_read_count'} >= $wt_read_cutoff);
		print WT_OUT "$sj_id:",$SJ{'wt'}{$sj_id}{'strand'},":",&SharedSS('wt', 'mut', $sj_id),"\n"; 
	}	
	# mut uniq
	foreach $sj_id (keys %{$SJ_GROUP{'mut'}}) {
		next unless ($SJ_GROUP{'mut'}{$sj_id}{'summed_read_count'} >= $mut_read_cutoff);
		print MUT_OUT "$sj_id:",$SJ{'mut'}{$sj_id}{'strand'},":",&SharedSS('mut','wt', $sj_id),"\n";
	}	
	# comm
	foreach $sj_id (keys %{$SJ_GROUP{'wt-mut'}}) {
		next unless ($SJ_GROUP{'wt-mut'}{$sj_id}{'summed_read_count'} >= $com_read_cutoff);
		print COM_OUT "$sj_id:",$SJ{'wt'}{$sj_id}{'strand'},":11\n"; # 11 means both 5' and 3' ss are common in wt and mut  
	}	
	close WT_OUT;
	close MUT_OUT;
	close COM_OUT;

	print "3 files (*.junc.uniq) generated in current location\n";
	#-----------------------------------------------------------------
}

sub SharedSS {
	# find if WT and MUT 5'& 3'are shared
	my ($sample1, $sample2, $sj_id) = @_;
	
	# get 5' and 3'ss for sample 1
	$five_ss = $SJ{$sample1}{$sj_id}{'5ss'}; 
	$three_ss = $SJ{$sample1}{$sj_id}{'3ss'}; 

	# compare with sample 2 
	my $shared5ss = my $shared3ss = 0;
	$shared5ss = 1 if (exists $SS5{$sample2}{$five_ss});
	$shared3ss = 1 if (exists $SS3{$sample2}{$three_ss});

	#print "$sample1:$sample2:$sj_id:5SS=$five_ss:3ss=$three_ss:$shared5ss:$shared3ss";<>;

	# get 5' and 3'ss for sample 2 and compare with sample 1
	#my $shared5ss = my $shared3ss = 0;
	#foreach my $sj_id2 (keys %{$SJ{$sample2}} ) {
	#	$five_ss2 = $SJ{$sample2}{$sj_id}{'5ss'}; 
	#	$three_ss2 = $SJ{$sample2}{$sj_id}{'3ss'}; 
	#	
	#	$shared5ss = 1 if ($five_ss1 == $five_ss2);
	#	$shared3ss = 1 if ($three_ss1 == $three_ss2);
	#}
	return ($shared5ss,$shared3ss);
}

sub ProcessTopHatSJ {
	print "Getting read strand from junctions.bed (TopHat) file...\n";
	# process junctions.bed to hash and store corresponding strand
	# (becuase finesplice does provide strand information)
	%TopHatSJ = ();
	open BEDLIST, "$bedlist";
	while ($filename = <BEDLIST>) {
		chomp $filename;
		open BED, $filename;
		<BED>;
		while ($line = <BED>) {
			chomp $line;
			@tmp = split /\t/, $line;
			($len1, $len2) = split /\,/, $tmp[10];
			$sj_x = $tmp[1] + $len1;
			$sj_y = $tmp[2] - $len2;
			$sj_id = "$tmp[0]:$sj_x:$sj_y";
			$TopHatSJ{$sj_id}{'strand'} = $tmp[5]; # sj_id with strand
			#print "$sj_id\n" if ($tmp[5] eq "");
			#print "$sj_id\t$tmp[5]\n";
		}	
		close BED;
	}
	close BEDLIST;
}	
