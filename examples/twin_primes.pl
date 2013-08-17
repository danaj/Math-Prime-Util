#!/usr/bin/env perl
use strict;
use warnings;

use Math::Prime::Util qw/-nobigint
                         prime_iterator next_prime
                         nth_prime_upper prime_precalc/;

my $count = shift || 20;

# Find twin primes (numbers where p and p+2 are prime)

# Time for the first 300k:
#   3m36s  Math::NumSeq::TwinPrimes (Perl 5.19.2 with v61)
#   1.4s   this iterator
#   7.9s   this iterator without the precalc

# This speeds things up, but isn't necessary.
my $estimate = 5000 + int( nth_prime_upper($count) * 1.4 * log($count) );
prime_precalc($estimate);

# Create a twin prime iterator using the prime_iterator construct
sub get_twin_prime_iterator1 {
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

# Create a twin prime iterator using the next_prime function
# A bit faster than the prime_iterator version.
sub get_twin_prime_iterator2 {
  my $start = shift || 2;
  my $prev = next_prime($start-1);
  my $p = next_prime($prev);
  return sub {
    do {
      ($prev, $p) = ($p, next_prime($p))
    } while ($p-$prev) != 2;
    $prev;
  };
}
my $twinit = get_twin_prime_iterator2();
for (1..$count) {
  print $twinit->(), "\n";
}
