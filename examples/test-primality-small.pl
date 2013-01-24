#!/usr/bin/env perl
use strict;
use warnings;
$| = 1;  # fast pipes

# Make sure the is_prob_prime functionality is working for small inputs.
# Good for making sure the first few M-R bases are set up correctly.

use Math::Prime::Util qw/is_prob_prime/;
use Math::Prime::Util::PrimeArray;

my @primes;  tie @primes, 'Math::Prime::Util::PrimeArray';

# Test just primes
if (0) {
  foreach my $i (1 .. 10000000) {
    my $n = shift @primes;
    die unless is_prob_prime($n);
    #print "." unless $i % 16384;
    print "$i $n\n" unless $i % 262144;
  }
}

# Test every number up to the 100Mth prime (about 2000M)
if (1) {
  die "2 should be prime" unless is_prob_prime(2);
  shift @primes;
  my $n = shift @primes;
  foreach my $i (2 .. 100_000_000) {
    die "$n should be prime" unless is_prob_prime($n);
    print "$i $n\n" unless $i % 262144;
    my $next = shift @primes;
    my $diff = ($next - $n) >> 1;
    if ($diff > 1) {
      foreach my $d (1 .. $diff-1) {
        my $cn = $n + 2*$d;
        die "$cn should be composite" if is_prob_prime($cn);
      }
    }
    $n = $next;
  }
}
