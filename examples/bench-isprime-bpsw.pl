#!/usr/bin/env perl
use strict;
use warnings;
$| = 1;  # fast pipes

use Math::Prime::Util;
use Math::Primality;

srand(500);
use bigint try=>'GMP';
use Math::BigInt::Random::OO;
#my $gen = Math::BigInt::Random::OO -> new(length => 50);
my $gen = Math::BigInt::Random::OO -> new(length => 25);

my @rns;
push @rns, $gen->generate() for (1 .. 100);

use Benchmark qw/:all/;
cmpthese(-.5, {
  "MP  MR" => sub { Math::Primality::is_strong_pseudoprime("$_","2") for @rns; },
  "MPU MR" => sub { Math::Prime::Util::PP::miller_rabin($_,2) for @rns; },
  "MP  LP" => sub { Math::Primality::is_strong_lucas_pseudoprime("$_") for @rns;},
  "MPU LP" => sub { Math::Prime::Util::PP::is_strong_lucas_pseudoprime($_) for @rns;},
  "MP  IP" => sub { Math::Primality::is_prime("$_") for @rns;},
  "MPU IP" => sub { Math::Prime::Util::PP::is_prime($_) for @rns;},
});
