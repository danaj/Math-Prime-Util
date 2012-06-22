#!/usr/bin/env perl
use strict;
use warnings;

use Math::Prime::Util qw/next_prime is_prime/;

my $count = shift || 20;

# Simple little loop looking for Sophie Germain primes (numbers where
# p and 2p+1 are both prime).  Calculating the first 100k runs 17x faster than
# Math::NumSeq::SophieGermainPrimes (about 3x faster if the latter's algorithm
# is changed).
my $prime = 2;
for (1..$count) {
  $prime = next_prime($prime) while (!is_prime(2*$prime+1));
  print "$prime\n";
  $prime = next_prime($prime);
}
