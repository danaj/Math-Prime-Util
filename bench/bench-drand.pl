use strict;
use warnings;

use Benchmark qw/cmpthese/;
use Math::Random::ISAAC;
use Math::Random::MT;
use Math::Random::MT::Auto;
use Math::Random::Xorshift;
use Math::Random::MTwist;
use Math::Random::Secure;
use ntheory;
use Crypt::PRNG;
# Could also use Data::Entropy::Algorithms but:
#   1) its dependencies have been broken for a while
#   2) it's really slow
# It is a nice idea, using AES counters.  Doubles are filled with 48 bits.

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

#                             Performance / Quality:
#   CORE::rand    29000k/s    ++++ / ---  drand48 has lots of issues
#   Xorshift      16000k/s    +++  / ---  old alg, only 32 bits
#   MTwist        14000k/s    +++  /  ++  
#   ntheory       13000k/s    +++  / +++  ISAAC CSPRNG
#   MT::Auto       4800k/s    +    /  ++  MTwist is faster
#   ISAAC          2400k/s    -    / ---  only 32 bits filled, wrong interval
#   MT             2200k/s    -    /  ++  MTwist is faster
#   Crypt::PRNG     705k/s    --   / +++  
#   Secure          426k/s    ---  / ---  only 32 bits filled
#
#  Also see  http://www.pcg-random.org/statistical-tests.html

cmpthese($trials, {
  # These are known to fail TestU01 SmallCrush
  'CORE::rand' => sub { CORE::rand for 1..1000 },
  'M::R::Xorshift->rand' => sub { $xor->rand for 1..1000 },
  'M::R::Xorshift::rand' => sub { Math::Random::Xorshift::rand for 1..1000 },

  # These only give 32 bits of output in their doubles!
  'M::R::ISAAC->rand' => sub { $isaac->rand for 1..1000 },
  'M::R::Secure::rand' => sub { Math::Random::Secure::rand for 1..1000 },

  'M::R::MT->rand' => sub { $mt->rand for 1..1000 },
  'M::R::MT::A::rand' => sub { Math::Random::MT::Auto::rand for 1..1000 },
  'M::R::MTwist::rand' => sub { Math::Random::MTwist::_rand for 1..1000 },
  'ntheory::drand' => sub { ntheory::drand for 1..1000 },
  'Crypt::PRNG::rand' => sub { Crypt::PRNG::rand for 1..1000 },
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
