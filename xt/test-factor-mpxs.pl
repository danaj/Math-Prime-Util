#!/usr/bin/env perl
use strict;
use warnings;
$| = 1;  # fast pipes

use Math::Prime::Util qw/factor urandomm/;
use Math::Factor::XS qw/prime_factors/;
use Config;

my $nlinear = 1000000;
my $nrandom = shift || 1000000;
my $randmax = ~0;
# MFXS is so slow on 17+ digit numbers, skip them.
$randmax = int(2**55) if $randmax > 2**55;

print "OK for first 1";
my $dig = 1;
my $i = 9;
foreach my $n (2 .. $nlinear) {
  my @mfxs = prime_factors($n);
  my @mpu  = factor($n);
  die "failure for $n" unless scalar @mfxs == scalar @mpu;
  for (0 .. $#mfxs) { die "failure for $n" unless $mfxs[$_] == $mpu[$_]; }
  if (--$i == 0) {
    print "0";
    $dig++;
    $i = (10 ** $dig) - (10 ** ($dig-1));
  }
}
print " numbers\n";
print "Testing random numbers from $nlinear to ", $randmax, "\n";

while ($nrandom-- > 0) {
  my $n = $nlinear + 1 + urandomm($randmax - $nlinear);
  my @mfxs = prime_factors($n);
  my @mpu  = factor($n);
  die "failure for $n" unless scalar @mfxs == scalar @mpu;
  for (0 .. $#mfxs) { die "failure for $n" unless $mfxs[$_] == $mpu[$_]; }
  print "." if ($nrandom % 256) == 0;
}
print "\n";
