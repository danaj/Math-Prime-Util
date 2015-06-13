#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/sum_primes vecsum primes/;

my %sums = (
  "189695660 to 189695892" => 0,
  "0 to 300000" => 3709507114,
  "12345 to 54321" => 132980191,
  "10000000 to 10001000" => 610034659,
);

plan tests => 1 + scalar(keys %sums);

{
  my @sum;
  my @exp;
  for (0..1000) {
    push @exp, vecsum( @{primes($_)} );
    push @sum, sum_primes($_);
  }
  is_deeply( \@sum, \@exp, "sum_primes for 0 to 1000" );
}
while (my($range, $expect) = each (%sums)) {
  my($low,$high) = $range =~ /(\d+) to (\d+)/;
  is( sum_primes($low,$high), $expect, "sum primes from $low to $high" );
}
