#!/usr/bin/env perl
use strict;
use warnings;
$| = 1;  # fast pipes

use Math::Prime::Util qw/-nobigint/;
use Math::Pari;

Math::Prime::Util::prime_precalc(10_000_000);
my $nlinear = 100000;
my $nrandom = shift || 100000;
my $randmax = ~0;

# Looks like MPU is 4x faster than Pari for 1 .. 1M, 5x faster for 1M - 10_000M.
#
print "OK for first 1";
my $dig = 1;
my $i = 9;
foreach my $n (2 .. $nlinear) {
  die "failure for eulerphi($n)" unless Math::Prime::Util::euler_phi($n) == Math::Pari::eulerphi($n);
  die "failure for moebius($n)" unless Math::Prime::Util::moebius($n) == Math::Pari::moebius($n);
  if (--$i == 0) {
    print "0";
    $dig++;
    $i = (10 ** $dig) - (10 ** ($dig-1));
  }
}
print " numbers\n";
print "Testing random numbers from $nlinear to ", $randmax, "\n";
my $mod = int($nrandom / 80);

while ($nrandom-- > 0) {
  my $n = $nlinear + 1 + int(rand($randmax - $nlinear));
  die "failure for eulerphi($n)" unless Math::Prime::Util::euler_phi($n) == Math::Pari::eulerphi($n);
  die "failure for moebius($n)" unless Math::Prime::Util::moebius($n) == Math::Pari::moebius($n);
  print "." if ($nrandom % $mod) == 0;
}
print "\n";
