#!/usr/bin/env perl
use warnings;
use strict;
use Math::Prime::Util qw/is_prime/;

my $limit = shift || 1000000;

# Any prime p where n^3 + n^2*p = m^3 must be the difference of (i+1)^3 - i^3.
# So we'll just walk them looking for primes.

my $sum = 0;
foreach my $i (1 .. 2650070) {
  my $j = $i+1;
  my $p = $j*$j*$j - $i*$i*$i;
  last if $p > $limit;
  $sum++ if is_prime($p);
}
print "$sum\n";
