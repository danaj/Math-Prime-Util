#!/usr/bin/env perl
use warnings;
use strict;
use Math::Prime::Util qw/:all/;
use List::MoreUtils qw/all/;

my $maxn = shift || 100_000_000;

# Takes ~ 5 seconds
my $sum = 0;
forprimes {
  my $n = $_-1;                  # 1+$n/1 is prime  (hence n=1 or even)
  if (is_prime(2+($n>>1))) {     # 2+$n/2 is prime  (noting n is even or 1)
    if (moebius($n) != 0) {      # n should be square free
      $sum += $n if all { is_prime($_+$n/$_) } divisors($n);
    }
  }
} $maxn;
print "$sum\n";

# We could additionally check these:
#   if ( (($n+2) % 4) == 0 || $n == 1) {
# Using all is more efficient, but this works:
#   $sum += $n unless scalar grep { !is_prime($_+$n/$_) } divisors($n);
# We could alternately generate primes to $maxn/2,
#    my $n = 2*$_-4; if (is_prime($n+1)) { ... }
