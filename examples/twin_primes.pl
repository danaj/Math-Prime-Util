#!/usr/bin/env perl
use strict;
use warnings;

use Math::Prime::Util qw/prime_iterator nth_prime_upper prime_precalc/;

my $count = shift || 20;

# Find twin primes (numbers where p and p+2 are prime)

# Time for the first 300k:
#   3m28s  Math::NumSeq::TwinPrimes
#   2.5s   this iterator
#   9.1s   this iterator without the precalc

# This speeds things up, but isn't necessary.
my $estimate = 5000 + int( nth_prime_upper($count) * 1.4 * log($count) );
prime_precalc($estimate);

sub get_twin_prime_iterator {
  my $p = shift || 2;
  my $it = prime_iterator($p);
  my $prev = $it->();    # prev = 2
  $p = $it->();          # p = 3
  return sub {
    do {
      ($prev, $p) = ($p, $it->())
    } while ($p-$prev) != 2;
    $prev;
  };
}
my $twinit = get_twin_prime_iterator();
for (1..$count) {
  print $twinit->(), "\n";
}
