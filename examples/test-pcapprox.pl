#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util qw/prime_count prime_count_approx prime_count_lower prime_count_upper LogarithmicIntegral RiemannR/;
$| = 1;  # fast pipes
my $use64 = Math::Prime::Util::_maxbits > 32;


my %pivals = (
                  10 => 4,
                 100 => 25,
                1000 => 168,
               10000 => 1229,
              100000 => 9592,
             1000000 => 78498,
            10000000 => 664579,
           100000000 => 5761455,
          1000000000 => 50847534,
         10000000000 => 455052511,
        100000000000 => 4118054813,
       1000000000000 => 37607912018,
      10000000000000 => 346065536839,
     100000000000000 => 3204941750802,
    1000000000000000 => 29844570422669,
   10000000000000000 => 279238341033925,
  100000000000000000 => 2623557157654233,
 1000000000000000000 => 24739954287740860,
10000000000000000000 => 234057667276344607,
);

printf("  N    %12s  %12s  %12s\n", "pc_approx", "Li", "R");
printf("-----  %12s  %12s  %12s\n", '-'x12,'-'x12,'-'x12);
foreach my $n (sort {$a<=>$b} keys %pivals) {
  my $pin  = $pivals{$n};
  my $pca  = prime_count_approx($n);

  my $pcli = ($n < 2) ? 0 : int(LogarithmicIntegral($n)-LogarithmicIntegral(2)+0.5);

  my $r = int(RiemannR($n)+0.5);

  printf "10^%2d  %12d  %12d  %12d\n", length($n)-1,
         abs($pca-$pin), abs($pcli-$pin), abs($r-$pin);
}

# Also see http://empslocal.ex.ac.uk/people/staff/mrwatkin/zeta/encoding1.htm
# for some ideas one how this could be made even more accurate.

print "\n";
print "Lower / Upper bounds.  Percentages.\n";
print "\n";

printf("  N    %12s  %12s  %12s  %12s\n", "lower", "upper", "SchoenfeldL", "SchoenfeldU");
printf("-----  %12s  %12s  %12s  %12s\n", '-'x12,'-'x12,'-'x12,'-'x12);
foreach my $n (sort {$a<=>$b} keys %pivals) {
  my $pin  = $pivals{$n};
  my $pcl  = prime_count_lower($n);
  my $pcu  = prime_count_upper($n);
  my ($scl,$scu) = schoenfeld($n);

  #printf "10^%2d  %12d  %12d\n", length($n)-1, $pin-$pcl, $pcu-$pin;
  printf "10^%2d  %12.7f  %12.7f  %12.7f  %12.7f\n",
         length($n)-1, 100*($pin-$pcl)/$pin, 100*($pcu-$pin)/$pin,
         100*($pin-$scl)/$pin, 100*($scu-$pin)/$pin;
}

sub schoenfeld {
  my $x = shift;
  my $lix = LogarithmicIntegral($x);
  my $bound = ((sqrt($x)*log($x)) / 8*3.1415926535);
  ($lix-$bound,$lix+$bound);
}

