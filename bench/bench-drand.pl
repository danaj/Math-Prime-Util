use strict;
use warnings;

use Benchmark qw/cmpthese/;
use Math::Random::ISAAC;           #           32-bit 2^32-1
use Math::Random::MT;              #           32-bit 2^32
use Math::Random::MT::Auto;        #           52-bit (x>>12)*2^-52+2^-53
use Math::Random::Xorshift;        #           32-bit 2^32-1
use Math::Random::MTwist;          #  :rand    52-bit x*2^-53
use Math::Random::Secure;          #           32-bit 2^32
use ntheory;                       #  :rand    NV bit x*2^-64
use Math::Prime::Util::GMP;        #           53+bit x*2^-64
use Crypt::PRNG;                   #           53?    (a*2^32+b)/2^53
                                   #  core     48-bit strong periods
# Could also use Data::Entropy::Algorithms but:
#   1) its dependencies have been broken for a while
#   2) it's really slow
# It is a nice idea, using AES counters.  Doubles are filled with only 48 bits.

my $trials = shift || -1;

# There isn't any good reason to expressly seed.
my $time = time;
srand($time);
Math::Random::Xorshift::srand($time);
ntheory::srand($time);
Math::Random::Secure::srand($time.$time.$time.$time);
my $isaac = Math::Random::ISAAC->new($time);
my $mt = Math::Random::MT->new($time);
my $xor = Math::Random::Xorshift->new($time);

use Math::Random::ISAAC::XS;  my $mrixs = Math::Random::ISAAC::XS->new($time);
use Math::Random::ISAAC::PP;  my $mripp = Math::Random::ISAAC::PP->new($time);

#                      Performance / Quality:
#   CORE::rand    29000k/s    ++++ / ---  drand48 has many bad points
#   Xorshift      16000k/s    +++  / ---  32-bit, old alg, closed interval
#   MTwist        14000k/s    +++  /  ++
#   MPU::GMP      14000k/s    +++  / +++  ISAAC CSPRNG
#   ntheory       12000k/s    +++  / +++  ChaCha20 CSPRNG
#   MT::Auto       4800k/s    +    /  ++  MTwist is faster
#   ISAAC          2400k/s    -    /  --  32-bit, bad seeding, closed interval
#   MT             2200k/s    -    /  ++  32-bit, MTwist is faster
#   Crypt::PRNG     705k/s    --   / +++
#   Secure          426k/s    ---  / ---  32-bit
#   ntheory PP      110k/s    ---- / +++  ChaCha20, very very slow
#
#  Also see  http://www.pcg-random.org/statistical-tests.html
#            https://blogs.unity3d.com/2015/01/07/a-primer-on-repeatable-random-numbers/

cmpthese($trials, {
  # These are known to fail TestU01 SmallCrush
  'CORE::rand' => sub { CORE::rand for 1..1000 },
  'M::R::Xorshift->rand' => sub { $xor->rand for 1..1000 },
  'M::R::Xorshift::rand' => sub { Math::Random::Xorshift::rand for 1..1000 },

  # doubles with only 32-bits of random data
  'M::R::ISAAC::XS' => sub { $mrixs->rand for 1..1000 },
  'M::R::ISAAC::PP' => sub { $mripp->rand for 1..1000 },
  'M::R::ISAAC->rand' => sub { $isaac->rand for 1..1000 },
  'M::R::Secure::rand' => sub { Math::Random::Secure::rand for 1..1000 },
  'M::R::MT->rand' => sub { $mt->rand for 1..1000 },

  # 52-bit, 53-bit doubles
  'M::R::MT::A::rand' => sub { Math::Random::MT::Auto::rand for 1..1000 },
  'M::R::MTwist::rand' => sub { Math::Random::MTwist::_rand for 1..1000 },
  'Crypt::PRNG::rand' => sub { Crypt::PRNG::rand for 1..1000 },
  # 53-bit or 64-bit NVs
  'MPU::GMP' => sub { Math::Prime::Util::GMP::drand for 1..1000 },

  # Fill all NV significand bits (24,53,64,113)
  'ntheory::drand' => sub { ntheory::drand for 1..1000 },
});

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
