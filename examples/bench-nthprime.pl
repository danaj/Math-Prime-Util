#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util qw/nth_prime prime_precalc/;
use Devel::TimeThis;
#prime_precalc(100000000);

foreach my $e (3 .. 9) {
  my $n = 10 ** $e;
  my $t = Devel::TimeThis->new("nth_prime(10^$e)");
  nth_prime($n);
}
