#!/usr/bin/env perl
use strict;
use warnings;

use Math::Prime::Util qw/random_prime random_ndigit_prime/;
use Benchmark qw/:all/;
use List::Util qw/min max/;
my $count = shift || -3;

srand(29);
test_at_digits($_) for (2 .. 10);

sub test_at_digits {
  my $digits = shift;
  die "Digits must be > 0" unless $digits > 0;

  cmpthese($count,{
    "$digits digits" => sub { random_ndigit_prime($digits) for (1..10000) },
  });
}
