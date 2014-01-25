#!/usr/bin/env perl
use strict;
use warnings;

# Find the first N abundant, deficient, or perfect numbers.

use Math::Prime::Util qw/divisor_sum next_prime is_prime/;

my $count = shift || 20;
my $type = lc(shift || 'abundant');

my $p = 0;
if ($type eq 'abundant') {
  while ($count-- > 0) {
    do { $p++ } while divisor_sum($p)-$p <= $p;
    print "$p\n";
  }
} elsif ($type eq 'deficient') {
  while ($count-- > 0) {
    do { $p++ } while divisor_sum($p)-$p >= $p;
    print "$p\n";
  }
} elsif ($type eq 'perfect') {
  # We'll use the chain of work by Euclid, Ibn al-Haytham, Euler, and others.
  # We just look for 2^(p-1)*(2^p-1) where 2^p-1 is prime.
  # Basically we're just finding Mersenne primes.
  # It's possible there are odd perfect numbers larger than 10^1500.
  do { require Math::BigInt;  Math::BigInt->import(try=>"GMP,Pari"); };
  while ($count-- > 0) {
    while (1) {
      $p = next_prime($p);
      last if is_prime(Math::BigInt->new(2)->bpow($p)->bdec);
    }
    print Math::BigInt->from_bin( '0b' . '1'x$p . '0'x($p-1) ), "\n";
  }
} else {
  die "Unknown type: $type\n";
}
