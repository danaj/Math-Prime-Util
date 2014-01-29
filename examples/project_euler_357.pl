#!/usr/bin/env perl
use warnings;
use strict;
use Math::Prime::Util qw/:all/;

# Takes ~ 10 seconds
my $sum = 0;
forprimes {
  my $n = $_-1;                  # 1+$n/1 is prime  (hence n=1 or even)
  if (is_prime(2+($n>>1))) {     # 2+$n/2 is prime  (noting n is even or 1)
    $sum += $n unless scalar grep { !is_prime($_+$n/$_) } divisors($n);
  }
} 100_000_000;
print "$sum\n";
