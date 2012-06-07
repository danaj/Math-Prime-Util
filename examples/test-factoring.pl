#!/usr/bin/perl -w
use strict;
use warnings;
$| = 1;  # fast pipes

use Math::Prime::Util qw/factor/;
use Math::Factor::XS qw/prime_factors/;

my $nlinear = 1000000;
my $nrandom = shift || 1000000;
my $randmax = ~0;

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
print "Testing random numbers from $nlinear to ", $randmax, "\n";

while ($nrandom-- > 0) {
  my $n = $nlinear + 1 + int(rand($randmax - $nlinear));
  my @mfxs = sort { $a<=>$b } prime_factors($n);
  my @mpu  = sort { $a<=>$b } factor($n);
  die "failure for $n" unless scalar @mfxs == scalar @mpu;
  for (0 .. $#mfxs) { die "failure for $n" unless $mfxs[$_] == $mpu[$_]; }
  print "." if ($nrandom % 1024) == 0;
}
print "\n";
