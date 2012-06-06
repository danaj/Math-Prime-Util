#!/usr/bin/perl -w
use strict;
use warnings;
$| = 1;  # fast pipes

use Math::Prime::Util qw/factor/;
use Math::Factor::XS qw/prime_factors/;

my $nlinear = 1000000;
my $nrandom = ~0;
my $randmax = 1_000_000_000_000;

print "OK for first 1";
my $dig = 1;
my $i = 9;
foreach my $n (2 .. $nlinear) {
  my @mfxs = sort { $a<=>$b } prime_factors($n);
  my @mpu  = sort { $a<=>$b } factor($n);
  die "failure for $n" unless scalar @mfxs == scalar @mpu;
  for (0 .. $#mfxs) { die "failure for $n" unless $mfxs[$_] == $mpu[$_]; }
  if (--$i == 0) {
    print "0";
    $dig++;
    $i = (10 ** $dig) - (10 ** ($dig-1));
  }
}
print " numbers\n";
print "Testing random numbers from $nlinear to ", $nlinear+$randmax, "\n";

while ($nrandom-- > 0) {
  my $n = $nlinear + int(rand($randmax));
  my @mfxs = sort { $a<=>$b } prime_factors($n);
  my @mpu  = sort { $a<=>$b } factor($n);
  die "failure for $n" unless scalar @mfxs == scalar @mpu;
  for (0 .. $#mfxs) { die "failure for $n" unless $mfxs[$_] == $mpu[$_]; }
  print "." if ($nrandom % 1024) == 0;
}
