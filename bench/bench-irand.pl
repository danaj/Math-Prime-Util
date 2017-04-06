use strict;
use warnings;

use Benchmark qw/cmpthese/;
use Math::Random::ISAAC;
use Math::Random::MT;
use Math::Random::MT::Auto;
use Math::Random::Xorshift;
use Math::Random::MTwist;
use Math::Prime::Util::GMP;
use Math::Prime::Util;
use Bytes::Random::Secure;
use Bytes::Random::Secure::Tiny;
use Crypt::PRNG;

my $time = time;
srand($time);

my $isaac = Math::Random::ISAAC->new($time);
my $mt = Math::Random::MT->new($time);

Math::Random::Xorshift::srand($time);
my $xor = Math::Random::Xorshift->new($time);

my $brs = Bytes::Random::Secure->new(NonBlocking=>1);
my $brst = Bytes::Random::Secure::Tiny->new(NonBlocking=>1);

cmpthese(-1, {
  'CORE::rand' => sub { int(CORE::rand(4294967295)) },
  'MRMT::irand' => sub { $mt->irand },
  'MRMTA::irand' => sub { Math::Random::MT::Auto::irand },
  'M::R::ISAAC#irand' => sub { $isaac->irand },
  'M::R::Xorshift#irand' => sub { $xor->irand },
  'BRS#irand' => sub { $brs->irand },
  'BRST#irand' => sub { $brst->irand },
  'M::R::Xorshift::irand' => sub { Math::Random::Xorshift::irand },
  'M::R::MTwist::irand32' => sub { Math::Random::MTwist::_irand32 },
  'MPUGMP::irand' => sub { Math::Prime::Util::GMP::irand },
  'MPU::irand' => sub { Math::Prime::Util::irand },
  'Crypt::PRNG::irand' => sub { Crypt::PRNG::irand },
}) if 0;

cmpthese(-1, {
  'CORE::rand' => sub { CORE::rand },
  #'M::R::ISAAC#rand' => sub { $isaac->rand },
  #'M::R::MT#rand' => sub { $mt->rand },
  #'M::R::Xorshift#rand' => sub { $xor->rand },
  'M::R::Xorshift::rand' => sub { Math::Random::Xorshift::rand },
  #'MRMTA::rand' => sub { Math::Random::MT::Auto::rand },
  'M::R::MTwist::rand' => sub { Math::Random::MTwist::_rand },
  #'MPUGMP::drand' => sub { Math::Prime::Util::GMP::drand },
  'MPU::drand' => sub { Math::Prime::Util::drand },
  #'Crypt::PRNG::rand' => sub { Crypt::PRNG::rand },
}) if 1;

# TestU01 SmallCrush on floating point output
# Passes
#   ntheory
#   Crypt::PRNG
#   Math::Random::MTwist
#   Math::Random::MT
#   Math::Random::MT::Auto
#   Math::Random::ISAAC   (32-bit, [0,1])
#   Math::Random::ISAAC::XS   (32-bit, [0,1])
#   Data::Entropy::Algorithms (AES, 48-bit)
#   Math::Random::Secure (ISAAC, 32-bit)
# Fails
#   5 CORE::rand
#   2 Math::Random::Xorshift   (32-bit, [0,1])

# perl -MMath::Random::Xorshift=rand -E 'say rand for 1..52000000' >/tmp/fr.txt
# bat2
# bat2.c:  bbattery_SmallCrushFile("/tmp/fr.txt");
