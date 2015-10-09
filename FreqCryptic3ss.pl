# script to calculate frequency of cryptic3ss on different positions in 
# upstream and downstream of canonocal 3'ss
use warnings;

# Input Files: 
#&Freq ('wt.cryptic3ss');
&Freq ('mut.cryptic3ss');
#&Freq ('wt-mut.cryptic3ss');

sub Freq {
	my $file = $_[0];
	%SS3 = ();
	open FILE, $file or die "cant open input $file\n";
	<FILE>;
	while ($line = <FILE>) {
		chomp $line;
		# pos; upstream or downstream postion of sample 3'ss w.r.t. ref. 3'ss position 
		# delta_psi; mut - wt
		($sample3ss, $ref3ss, $dist, $pos, $pval, $delta_psi) = split /\t/, $line;
		if ($delta_psi > 0) { 		
			push @{$SS3{$pos}{'level_up'}}, $dist; 
		} elsif ($delta_psi < 0 ) {
			push @{$SS3{$pos}{'level_down'}}, $dist; 
		}	
	}
	close FILE;
	# bin the distribution of position of 3'ss
	$window = 10;
	foreach $pos (sort keys %SS3) {
		foreach $level (sort keys %{$SS3{$pos}}) {
			print "$pos:$level\n";
			@sorted_dist = sort {$a <=> $b} @{$SS3{$pos}{$level}};
			$longest_dist = $sorted_dist[-1];
			#--------------  get frequencies ------------
			%FREQ = @RANGE = ();
			$i = 0;
			# group frequencies in varying size windows
			while ($i < $longest_dist) {
				$x = $i + 1;
				$y = $i + $window;
				$range = "$x..$y"; 
				push @RANGE, $range;
				foreach $dist (@sorted_dist) {
					$FREQ{$range} += 1 if ($dist > $i && $dist <= $i + $window); 
				}
				$i += $window;
			}
			# print freq's from different windows 
			foreach $range (@RANGE) {
				next unless (exists $FREQ{$range});
				print "$range\t", $FREQ{$range},"\n";
			}

			#print "$pos\t$longest_dist";<>;
		}	
	}

}
