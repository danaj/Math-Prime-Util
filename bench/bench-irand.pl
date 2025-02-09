use strict;
use warnings;
no warnings 'void';

use Benchmark qw/cmpthese/;
use Math::Random::ISAAC;
use Math::Random::ISAAC::XS;
use Math::Random::MT;
use Math::Random::MT::Auto;
use Math::Random::Xorshift;
#use Math::Random::MTwist;
use Math::Random::PCG32;
use Math::Prime::Util::GMP;
use Math::Prime::Util;
use Bytes::Random::Secure;
use Bytes::Random::Secure::Tiny;
use Crypt::PRNG;

my $trials = shift || -1;

# There is no real point in seeding here.
my $time = time;
srand($time);
Math::Random::Xorshift::srand($time);
Math::Prime::Util::srand($time);
my $isaac = Math::Random::ISAAC::XS->new($time);
#my $mt = Math::Random::MT->new($time);
my $xor = Math::Random::Xorshift->new($time);
my $brs = Bytes::Random::Secure->new(NonBlocking=>1);
my $brst = Bytes::Random::Secure::Tiny->new(NonBlocking=>1);
my $pcg = Math::Random::PCG32->new(42,54);

cmpthese($trials, {
  'CORE::rand' => sub { int(CORE::rand(4294967295)) for 1..1000; },
  #'MRMT::irand' => sub { $mt->irand for 1..1000; },
  'M::R::ISAAC#irand' => sub { $isaac->irand for 1..1000; },
  'M::R::Xorshift::irand' => sub { Math::Random::Xorshift::irand for 1..1000; },
#  'M::R::Xorshift#irand' => sub { $xor->irand for 1..1000; },
  'BRS#irand' => sub { $brs->irand for 1..1000; },
  'BRST#irand' => sub { $brst->irand for 1..1000; },
  'Crypt::PRNG::irand' => sub { Crypt::PRNG::irand for 1..1000; },

  'MRMTA::irand' => sub { (0xFFFFFFFF & Math::Random::MT::Auto::irand) for 1..1000; },
  #'M::R::MTwist::irand32' => sub { Math::Random::MTwist::_irand32 for 1..1000; },
  'M::R::PCG32::rand' => sub { $pcg->rand for 1..1000; },
  'MPUGMP::irand' => sub { Math::Prime::Util::GMP::irand for 1..1000; },
  'MPU::irand' => sub { Math::Prime::Util::irand for 1..1000; },
}) if 1;

cmpthese($trials, {

  # These have only 32-bit irand, nothing else.
  #'CORE::rand' => sub { int(CORE::rand(4294967295)) for 1..1000; },
  #'MRMT::irand 32x2' => sub { (($mt->irand <<32)|$mt->irand) for 1..1000; },
  'M::R::ISAAC 32x2' => sub { (($isaac->irand <<32)|$isaac->irand) for 1..1000; },
  'M::R::PCG32 32x2' => sub { (($pcg->rand <<32)|$pcg->rand) for 1..1000; },
  #'M::R::Xorshift#irand' => sub { $xor->irand for 1..1000; },
  'M::R::Xorshift 32x2' => sub { ((Math::Random::Xorshift::irand <<32)|Math::Random::Xorshift::irand) for 1..1000; },

  # These don't have 64-bit irand functions, but have random bytes.
  # Select the fastest of the two options.
  #'BRS prb' => sub { unpack("Q",$brs->bytes(8)) for 1..1000; },
  'BRS 32x2' => sub { (($brs->irand << 32)|$brs->irand) for 1..1000; },
  #'BRST prb' => sub { unpack("Q",$brst->bytes(8)) for 1..1000; },
  'BRST 32x2' => sub { (($brst->irand << 32)|$brst->irand) for 1..1000; },
  'Crypt::PRNG prb' => sub { unpack("Q",Crypt::PRNG::random_bytes(8)) for 1..1000; },
  #'Crypt::PRNG 32x2' => sub { ((Crypt::PRNG::irand << 32)|Crypt::PRNG::irand) for 1..1000; },

  # These have 64-bit irand functions
  'MRMTA::irand' => sub { Math::Random::MT::Auto::irand for 1..1000; },
  #'M::R::MTwist::irand64' => sub { Math::Random::MTwist::_irand64 for 1..1000; },
  'MPUGMP::irand64' => sub { Math::Prime::Util::GMP::irand64 for 1..1000; },
  'MPU::irand64' => sub { Math::Prime::Util::irand64 for 1..1000; },
  #'MPU prb' => sub { unpack("Q",Math::Prime::Util::random_bytes(8)) for 1..1000; },
  #'MPU 32x2' => sub { ((Math::Prime::Util::irand << 32)|Math::Prime::Util::irand) for 1..1000; },
}) if 1;
