#!/usr/bin/env perl
use strict;
use warnings;

use Math::Prime::Util;
use Math::Prime::Util::GMP;
use Bytes::Random::XS;
use Bytes::Random;
use Bytes::Random::Secure;
use Math::Random::MTwist;
#use Data::Entropy::Algorithms;
use Crypt::PRNG;
#use Crypt::Random;
use Benchmark qw/:all/;

# On a Macbook early 2015, random_bytes() can pump out over 1GB/s in C.

my $trial = -1;

if (0) {
 print "#  8 random bytes\n";
 cmpthese($trial,{
  "MPU"      => sub { Math::Prime::Util::random_bytes(8); },
  "MPU::GMP" => sub { Math::Prime::Util::GMP::random_bytes(8); },
  "BRXS"     => sub { Bytes::Random::XS::random_bytes(8); },
  "BR"       => sub { Bytes::Random::random_bytes(8); },
  "BRS"      => sub { Bytes::Random::Secure::random_bytes(8); },
  "MRMT"     => sub { Math::Random::MTwist::_randstr(8); },
#  "DEA"      => sub { Data::Entropy::Algorithms::rand_bits(8*8); },
  "Crypt::PRNG"   => sub { Crypt::PRNG::random_bytes(8); },
  "rand"     => sub { pack('C*', map { int(rand 256) } 1..8); },
#  "Crypt::Random"   => sub { Crypt::Random::makerandom_octet(Length=>8,Strength=>0); },
 });
}
if (1) {
 print "#  256 random bytes\n";
 cmpthese($trial,{
  "MPU"      => sub { Math::Prime::Util::random_bytes(256); },
  "MPU::GMP" => sub { Math::Prime::Util::GMP::random_bytes(256); },
  "BRXS"     => sub { Bytes::Random::XS::random_bytes(256); },
  "BR"       => sub { Bytes::Random::random_bytes(256); },
  "BRS"      => sub { Bytes::Random::Secure::random_bytes(256); },
  "MRMT"     => sub { Math::Random::MTwist::_randstr(256); },
#  "DEA"      => sub { Data::Entropy::Algorithms::rand_bits(8*256); },
  "CryptX"   => sub { Crypt::PRNG::random_bytes(256); },
  "rand"     => sub { pack('C*', map { int(rand 256) } 1..256); },
#  "Crypt::Random"   => sub { Crypt::Random::makerandom_octet(Length=>256,Strength=>0); },
 });
}
if (0) {
 print "#  16384 random bytes\n";
 cmpthese($trial,{
  "MPU"      => sub { Math::Prime::Util::random_bytes(16384); },
  "MPU::GMP" => sub { Math::Prime::Util::GMP::random_bytes(16384); },
  "BRXS"     => sub { Bytes::Random::XS::random_bytes(16384); },
  "BR"       => sub { Bytes::Random::random_bytes(16384); },
  "BRS"      => sub { Bytes::Random::Secure::random_bytes(16384); },
  "MRMT"     => sub { Math::Random::MTwist::_randstr(16384); },
#  "DEA"      => sub { Data::Entropy::Algorithms::rand_bits(8*16384); },
  "CryptX"   => sub { Crypt::PRNG::random_bytes(16384); },
  "rand"     => sub { pack('C*', map { int(rand 256) } 1..16384); },
#  "Crypt::Random"   => sub { Crypt::Random::makerandom_octet(Length=>16384,Strength=>0); },
 });
}
