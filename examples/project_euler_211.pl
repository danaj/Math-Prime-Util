#!/usr/bin/env perl
use warnings;
use strict;
use Math::Prime::Util qw/:all/;

# Brute force using MPU's divisor_sum.
# MPU v0.38   1.5 minutes
# Pari        3.5 minutes:
#   s=0; for(n=1,64000000-1,if(issquare(sigma(n,2)),s=s+n;))

my $n = shift || 64_000_000;

my $sum = 0;
foreach my $i (0 .. $n-1) {
  $sum += $i if is_power( divisor_sum($i, 2) , 2);
}
print "$sum\n";
