#!/usr/bin/env perl
use strict;
use warnings;

use Math::Prime::Util qw/random_nbit_prime/;
use Math::BigInt try=>'GMP';
use Benchmark qw/:all/;
use List::Util qw/min max/;
my $count = shift || -3;

srand(29);
test_at_bits($_) for (15, 30, 60, 128, 256, 512, 1024, 2048, 4096);

sub test_at_bits {
  my $bits = shift;
  die "Digits must be > 0" unless $bits > 0;

  cmpthese($count,{
    "$bits bits" => sub { random_nbit_prime($bits); },
  });
}
