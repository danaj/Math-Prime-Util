#!/usr/bin/env perl
use strict;
use warnings;
$| = 1;  # fast pipes

# Make sure the is_prob_prime functionality is working for small inputs.
# Good for making sure the first few M-R bases are set up correctly.
my $limit = 800_000_000;

use Math::Prime::Util qw/is_prob_prime/;
# Use another code base for comparison.
# Math::Prime::FastSieve is very fast -- far faster than Math::Primality
use Math::Prime::FastSieve;
my $sieve = Math::Prime::FastSieve::Sieve->new($limit + 10_000);

if (0) {  # just primes using Math::Prime::FastSieve
  my $n = 2;
  my $i = 1;
  while ($n < $limit) {
    die "$n" unless is_prob_prime($n);
    $n = $sieve->nearest_ge( $n+1 );
    print "$i $n\n" unless $i++ % 16384;
  }
}

# Test every number up to the 100Mth prime (about 2000M)
if (1) {
  my $n = 2;
  my $i = 1;
  while ($n <= $limit) {
    die "$n should be prime" unless is_prob_prime($n);
    print "$i $n\n" unless $i++ % 262144;
    my $next = $sieve->nearest_ge( $n+1 );
    my $diff = ($next - $n) >> 1;
    if ($diff > 1) {
      foreach my $d (1 .. $diff-1) {
        my $cn = $n + 2*$d;
        die "$cn should be composite" if is_prob_prime($cn);
      }
    }
    $n = $next;
  }
  print "Success to $limit!\n";
}
