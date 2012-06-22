#!/usr/bin/env perl
use strict;
use warnings;

use Math::Prime::Util qw/next_prime is_prime/;

my $count = shift || 20;

# Simple little loop looking for twin primes (numbers where p and p+2 are
# both prime).  Print them both.  About 3x faster than Math::NumSeq::TwinPrimes.
my $prime = 2;
my $next;
for (1..$count) {
  while (1) {
    $next = next_prime($prime);
    last if ($next-2) == $prime;
    $prime = $next;
  }
  print $prime, ",", $prime+2, "\n";
  $prime = $next;
}
