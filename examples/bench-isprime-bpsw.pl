#!/usr/bin/env perl
use strict;
use warnings;
$| = 1;  # fast pipes

use Math::Prime::Util;
use Math::Primality;

srand(500);
#use bigint lib=>'';
use Math::BigInt::Random::OO;
my $gen = Math::BigInt::Random::OO -> new(length => 180);
#my $gen = Math::BigInt::Random::OO -> new(length => 8);

my @rns;
push @rns, $gen->generate() for (1 .. 50);
#my @rns;
#push @rns, 1000000000 + int(rand(1000000000)) for (1..100);

print "Verifying";
for my $n (@rns) {
  die "bad MR for $n" unless Math::Prime::Util::PP::miller_rabin($n,2) == Math::Primality::is_strong_pseudoprime("$n","2");
  die "bad LP for $n" unless Math::Prime::Util::PP::is_strong_lucas_pseudoprime($n) == Math::Primality::is_strong_lucas_pseudoprime("$n");
  die "bad IP for $n" unless (Math::Prime::Util::PP::is_prime($n)?1:0) == (Math::Primality::is_prime("$n")?1:0);
  print ".";
}
print "OK\n";

use Benchmark qw/:all/;
my $sum = 0;
cmpthese(-.5, {
  "MP  MR" => sub { $sum += Math::Primality::is_strong_pseudoprime("$_","2") for @rns; },
  "MPU MR" => sub { $sum += Math::Prime::Util::PP::miller_rabin($_,2) for @rns; },
  #"MPUxMR" => sub { Math::Prime::Util::miller_rabin($_,2) for @rns; },
  "MP  LP" => sub { $sum += Math::Primality::is_strong_lucas_pseudoprime("$_") for @rns;},
  "MPU LP" => sub { $sum += Math::Prime::Util::PP::is_strong_lucas_pseudoprime($_) for @rns;},
  "MP  IP" => sub { $sum += Math::Primality::is_prime("$_") for @rns;},
  "MPU IP" => sub { $sum += Math::Prime::Util::PP::is_prime($_) for @rns;},
  #"MPUxIP" => sub { Math::Prime::Util::is_prime($_) for @rns;},
});
