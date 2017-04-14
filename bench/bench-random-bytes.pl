#!/usr/bin/env perl
use strict;
use warnings;

use Math::Prime::Util;
use Math::Prime::Util::GMP;
use Math::Prime::Util::ChaCha;
use Math::Prime::Util::ISAAC;
use Bytes::Random::XS;          # XS version calling Drand01()
use Bytes::Random;              # Just a loop around CORE::rand!
use Bytes::Random::Secure;
use Math::Random::MTwist;
use Crypt::PRNG;
use Crypt::Random;
use Data::Entropy::Algorithms;
use Benchmark qw/:all/;

Math::Prime::Util::ISAAC::srand;
Math::Prime::Util::ChaCha::srand;
use Math::Random::ISAAC::PP; my $mripp=Math::Random::ISAAC::PP->new();
use Math::Random::ISAAC::XS; my $mrixs=Math::Random::ISAAC::XS->new();

# On a Macbook early 2015, the fastest XS methods can pump out over 1 GB/s.
# In theory we could see 4+ GB/s from a module.

my $trial = shift || -1;

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
if (0) {
 print "#  256 random bytes\n";
 cmpthese($trial,{
  "MPU"      => sub { Math::Prime::Util::random_bytes(256); },
  "MPU X"    => sub { Math::Prime::Util::GMP::random_bytes(256); },
  "BRXS"     => sub { Bytes::Random::XS::random_bytes(256); },
  "PP MPU ISAAC" => sub { Math::Prime::Util::ISAAC::random_bytes(256); },
  "PP MPU ChaCha" => sub { Math::Prime::Util::ChaCha::random_bytes(256); },
  "PP MR ISAAC" => sub { pack("L*",map{$mripp->irand}1..64); },
  "XS MR ISAAC" => sub { pack("L*",map{$mrixs->irand}1..64); },
  "XS MR ISAAC2" => sub { my$s='';$s.=pack("L",$mrixs->irand)for 1..64;$s; },
  "CryptX"   => sub { Crypt::PRNG::random_bytes(256); },
  "BRS"      => sub { Bytes::Random::Secure::random_bytes(256); },
  "MRMTwist" => sub { Math::Random::MTwist::_randstr(256); },
  "rand"     => sub { pack('L*', map { int(rand 4294967296) } 1..64); },
  "DEA"      => sub { Data::Entropy::Algorithms::rand_bits(8*256); },
  "Crypt::Random"   => sub { Crypt::Random::makerandom_octet(Length=>256,Strength=>0); },
  "BR"       => sub { Bytes::Random::random_bytes(256); },
 });
}
if (0) {
 print "#  16384 random bytes\n";
 cmpthese($trial,{
  "MPU"      => sub { Math::Prime::Util::random_bytes(16384); },
  "MPU X"    => sub { Math::Prime::Util::xrandom_bytes(16384); },
  #"MPU::GMP" => sub { Math::Prime::Util::GMP::random_bytes(16384); },
  "BRXS"     => sub { Bytes::Random::XS::random_bytes(16384); },
  #"BR"       => sub { Bytes::Random::random_bytes(16384); },
  #"BRS"      => sub { Bytes::Random::Secure::random_bytes(16384); },
  #"MRMT"     => sub { Math::Random::MTwist::_randstr(16384); },
#  "DEA"      => sub { Data::Entropy::Algorithms::rand_bits(8*16384); },
  #"CryptX"   => sub { Crypt::PRNG::random_bytes(16384); },
  #"rand"     => sub { pack('C*', map { int(rand 256) } 1..16384); },
#  "Crypt::Random"   => sub { Crypt::Random::makerandom_octet(Length=>16384,Strength=>0); },
 });
}
if (1) {
 print "#  64k random bytes\n";
 cmpthese($trial,{
  "entropy"  => sub { Math::Prime::Util::entropy_bytes(64*1024); },
  "MPU"      => sub { Math::Prime::Util::random_bytes(64*1024); },
  "MPU X"    => sub { Math::Prime::Util::GMP::random_bytes(64*1024); },
  "BRXS"     => sub { Bytes::Random::XS::random_bytes(64*1024); },
  "PP MPU ISAAC" => sub { Math::Prime::Util::ISAAC::random_bytes(64*1024); },
  "PP MPU ChaCha" => sub { Math::Prime::Util::ChaCha::random_bytes(64*1024); },
  "PP MR ISAAC" => sub { pack("L*",map{$mripp->irand}1..16384); },
  "XS MR ISAAC" => sub { pack("L*",map{$mrixs->irand}1..16384); },
  "XS MR ISAAC2" => sub { my$s='';$s.=pack("L",$mrixs->irand)for 1..16384;$s; },
  "CryptX"   => sub { Crypt::PRNG::random_bytes(64*1024); },
  "BRS"      => sub { Bytes::Random::Secure::random_bytes(64*1024); },
  "MRMTwist" => sub { Math::Random::MTwist::_randstr(64*1024); },
  "rand"     => sub { pack('L*', map { int(rand 4294967296) } 1..16384); },
  #"DEA"      => sub { Data::Entropy::Algorithms::rand_bits(8*65536); },
  #"Crypt::Random"   => sub { Crypt::Random::makerandom_octet(Length=>65536,Strength=>0); },
  "BR"       => sub { Bytes::Random::random_bytes(64*1024); },
 });
}
