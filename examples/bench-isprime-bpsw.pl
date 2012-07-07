#!/usr/bin/env perl
use strict;
use warnings;
$| = 1;  # fast pipes

use Math::Prime::Util;
use Math::Primality;

# GMP is ~3x faster than Calc or Pari for these operations
use bigint try=>'GMP';
srand(500);
use Math::BigInt::Random::OO;
my $gen = Math::BigInt::Random::OO -> new(length => 80);
#my $gen = Math::BigInt::Random::OO -> new(length => 8);

my @rns;
while (@rns < 50) {
  # Ensure $n is an object of our bigint class, not MBROO's choice.
  my $n = Math::BigInt->new( $gen->generate()->bstr );
  $n++ if ($n % 2) == 0;    # Math::BigInt::Random::OO keeps making evens (bug?)
  next unless ($n % 2) != 0;
  push @rns, $n;
}
map { $_ = int($_->bstr) if $_ <= ~0 } @rns;
#print "$_\n" for @rns;
no bigint;  # Benchmark doesn't work with bigint on.

print "Verifying";
for my $n (@rns) {
  die "bad MR for $n" unless Math::Prime::Util::is_strong_pseudoprime($n,2) == Math::Primality::is_strong_pseudoprime("$n","2");
  die "bad LP for $n" unless Math::Prime::Util::is_strong_lucas_pseudoprime($n) == Math::Primality::is_strong_lucas_pseudoprime("$n");
  die "bad IP for $n" unless (Math::Prime::Util::is_prime($n)?1:0) == (Math::Primality::is_prime("$n")?1:0);
  print ".";
}
print "OK\n";

use Benchmark qw/:all/;
my $sum = 0;
cmpthese(-2, {
  "MP  MR" => sub { $sum += Math::Primality::is_strong_pseudoprime("$_","2") for @rns; },
  "MPU MR" => sub { $sum += Math::Prime::Util::is_strong_pseudoprime($_,2) for @rns; },
  #"MPUxMR" => sub { Math::Prime::Util::miller_rabin($_,2) for @rns; },
  "MP  LP" => sub { $sum += Math::Primality::is_strong_lucas_pseudoprime("$_") for @rns;},
  "MPU LP" => sub { $sum += Math::Prime::Util::is_strong_lucas_pseudoprime($_) for @rns;},
  "MP  IP" => sub { $sum += Math::Primality::is_prime("$_") for @rns;},
  "MPU IP" => sub { $sum += Math::Prime::Util::is_prime($_) for @rns;},
  #"MPUxIP" => sub { Math::Prime::Util::is_prime($_) for @rns;},
});
