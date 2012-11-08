#!/usr/bin/env perl
use strict;
use warnings;
$| = 1;  # fast pipes

use Math::Prime::Util;
use Math::Primality;
use Config;

my $nlinear = 10000;
my $nrandom = shift || 20000;
my $randmax = ~0;

# I was using Math::BigInt::Random::OO, but on my machine:
#     my $gen = Math::BigInt::Random::OO -> new(length => 23);
# generates only even numbers.
my $rgen = sub {
  my $range = shift;
  return 0 if $range <= 0;
  my $rbits = 0; { my $t = $range; while ($t) { $rbits++; $t >>= 1; } }
  while (1) {
    my $rbitsleft = $rbits;
    my $U = $range - $range;  # 0 or bigint 0
    while ($rbitsleft > 0) {
      my $usebits = ($rbitsleft > $Config{randbits}) ? $Config{randbits} : $rbitsleft;
      $U = ($U << $usebits) + int(rand(1 << $usebits));
      $rbitsleft -= $usebits;
    }
    return $U if $U <= $range;
  }
};
my $rand_ndigit_gen = sub {
  my $digits = shift;
  die "Digits must be > 0" unless $digits > 0;
  my $howmany = shift || 1;
  my ($base, $max);

  if ( 10**$digits < ~0) {
    $base = ($digits == 1) ? 0 : int(10 ** ($digits-1));
    $max = int(10 ** $digits);
    $max = ~0 if $max > ~0;
  } else {
    $base = Math::BigInt->new(10)->bpow($digits-1);
    $max = Math::BigInt->new(10)->bpow($digits) - 1;
  }
  my @nums = map { $base + $rgen->($max-$base) } (1 .. $howmany);
  return (wantarray) ? @nums : $nums[0];
};



if (1) {
print "OK for first 1";
my $dig = 1;
my $i = 9;
foreach my $n (2 .. $nlinear) {
  die "MR(2) failure for $n" unless Math::Prime::Util::is_strong_pseudoprime($n,2) == Math::Primality::is_strong_pseudoprime($n,2);
  die "SLPSP failure for $n" unless Math::Prime::Util::is_strong_lucas_pseudoprime($n) == Math::Primality::is_strong_lucas_pseudoprime($n);
  die "Prime failure for $n" unless Math::Prime::Util::is_prime($n) == Math::Primality::is_prime($n);
  if (--$i == 0) {
    print "0";
    $dig++;
    $i = (10 ** $dig) - (10 ** ($dig-1));
  }
}
print " numbers\n";
print "Testing random numbers from $nlinear to ", $randmax, "\n";

foreach my $r (1 .. $nrandom) {
  my $n = $nlinear + 1 + int(rand($randmax - $nlinear));
  my $rand_base = 2 + $rgen->($n-4);
  die "MR(2) failure for $n" unless Math::Prime::Util::is_strong_pseudoprime($n,2) == Math::Primality::is_strong_pseudoprime($n,2);
  die "MR($rand_base) failure for $n" unless Math::Prime::Util::is_strong_pseudoprime($n,$rand_base) == Math::Primality::is_strong_pseudoprime($n,$rand_base);
  die "SLPSP failure for $n" unless Math::Prime::Util::is_strong_lucas_pseudoprime($n) == Math::Primality::is_strong_lucas_pseudoprime($n);
  die "Prime failure for $n" unless (Math::Prime::Util::is_prime($n)?1:0) == Math::Primality::is_prime($n);
  print "." if ($r % 256) == 0;
}
print "\n";

use bigint try => 'GMP';
my $big_base = 2**64 + 1;
my $range = 2**1024 - 1;
my $end_base = $big_base + $range;
print "Testing random numbers from $big_base to $end_base\n";

foreach my $r (1 .. int($nrandom/50)) {
  my $n = $big_base + $rgen->($range);
  my $rand_base = 2 + $rgen->($n-4);
  die "MR(2) failure for $n" unless Math::Prime::Util::is_strong_pseudoprime($n,2) == Math::Primality::is_strong_pseudoprime("$n","2");
  die "MR($rand_base) failure for $n" unless Math::Prime::Util::is_strong_pseudoprime($n,$rand_base) == Math::Primality::is_strong_pseudoprime($n,$rand_base);
  die "SLPSP failure for $n" unless Math::Prime::Util::is_strong_lucas_pseudoprime($n) == Math::Primality::is_strong_lucas_pseudoprime("$n");
  die "Prime failure for $n" unless (Math::Prime::Util::is_prime($n)?1:0)  == Math::Primality::is_prime("$n");
  #print "SUCCESS with $n\n";
  print "." if ($r % 16) == 0;
}
print "\n";
}

use bigint try => 'GMP';
my $num_rns = 50;
my $len_rns = 100;
my $count = -.1;

my @rns;  # make the primality tests at least lift a finger.
while (@rns < $num_rns) {
  my $n = $rand_ndigit_gen->($len_rns);
  next unless $n%2 && $n%3 && $n%5 && $n%7 && $n%11 && $n%13;
  push @rns, $n;
}

use Benchmark qw/:all/;
print "Starting benchmarks, $num_rns $len_rns-digit random numbers...\n";

if (1) {
  print "\nMiller-Rabin, one base:\n";
  cmpthese($count, {
    "MPU:PP"  => sub { Math::Prime::Util::PP::miller_rabin($_,2) for @rns; },
    "MPU:GMP" => sub { Math::Prime::Util::GMP::is_strong_pseudoprime($_,2) for @rns; },
    "MPU"     => sub { Math::Prime::Util::is_strong_pseudoprime($_,2) for @rns; },
    "MP"      => sub { Math::Primality::is_strong_pseudoprime("$_","2") for @rns; },
  });
}

if (1) {
  print "\nStrong Lucas test:\n";
  cmpthese($count, {
    "MPU:PP"  => sub { Math::Prime::Util::PP::is_strong_lucas_pseudoprime($_) for @rns;},
    "MPU:GMP" => sub { Math::Prime::Util::GMP::is_strong_lucas_pseudoprime($_) for @rns;},
    "MPU"     => sub { Math::Prime::Util::is_strong_lucas_pseudoprime($_) for @rns;},
    "MP"      => sub { Math::Primality::is_strong_lucas_pseudoprime("$_") for @rns;},
  });
}

if (1) {
  print "\nBPSW test:\n";
  cmpthese($count, {
    "MPU:PP"  => sub { my $sum = 0;
               do { $sum += ( Math::Prime::Util::PP::miller_rabin($_, 2) &&
                      Math::Prime::Util::PP::is_strong_lucas_pseudoprime($_) )
                      ? 1 : 0 } for @rns; },
    "MPU:GMP" => sub { Math::Prime::Util::GMP::is_prob_prime($_) for @rns; },
    "MPU"     => sub { Math::Prime::Util::is_prob_prime($_) for @rns;},
    "MP"      => sub { Math::Primality::is_prime("$_") for @rns;},
  });
}
