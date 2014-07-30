#!/usr/bin/env perl
use strict;
use warnings;
$| = 1;  # fast pipes

use Math::Prime::Util;
use Math::Primality;

my $count = shift || -1;

# GMP is ~3x faster than Calc or Pari for these operations
use bigint try=>'GMP';
srand(500);
use Config;

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

my @rns;
while (@rns < 50) {
  my $n = $rgen->( Math::BigInt->new(2)->bpow(81) );
  $n++ if ($n % 2) == 0;
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
cmpthese($count, {
  "MP  MR" => sub { $sum += Math::Primality::is_strong_pseudoprime("$_","2") for @rns; },
  "MPU MR" => sub { $sum += Math::Prime::Util::GMP::is_strong_pseudoprime($_,2) for @rns; },
  #"MPUxMR" => sub { Math::Prime::Util::is_strong_pseudoprime($_,2) for @rns; },
  "MP  LP" => sub { $sum += Math::Primality::is_strong_lucas_pseudoprime("$_") for @rns;},
  "MPU LP" => sub { $sum += Math::Prime::Util::GMP::is_strong_lucas_pseudoprime($_) for @rns;},
  "MPU ELP" => sub { $sum += Math::Prime::Util::GMP::is_extra_strong_lucas_pseudoprime($_) for @rns;},
  "MPU AELP" => sub { $sum += Math::Prime::Util::GMP::is_almost_extra_strong_lucas_pseudoprime($_) for @rns;},
  "MP  IP" => sub { $sum += Math::Primality::is_prime("$_") for @rns;},
  "MPU IP" => sub { $sum += Math::Prime::Util::is_prime($_) for @rns;},
});
