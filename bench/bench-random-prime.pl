#!/usr/bin/env perl
use strict;
use warnings;

use Math::Prime::Util qw/-nobigint random_prime random_ndigit_prime/;
use Benchmark qw/:all/;
use List::Util qw/min max/;
my $count = shift || -3;
my $maxdigits = (~0 <= 4294967295) ? 10 : 20;

srand(29);
test_at_digits($_) for (2 .. $maxdigits);

sub test_at_digits {
  my $digits = shift;
  die "Digits must be > 0" unless $digits > 0;

  cmpthese($count,{
    "$digits digits" => sub { random_ndigit_prime($digits) for (1..1000) },
  });
}
