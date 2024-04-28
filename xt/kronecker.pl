#!/usr/bin/env perl
use strict;
use warnings;
$| = 1;  # fast pipes

use Math::Prime::Util qw/kronecker urandomb negint/;
use Math::Prime::Util::PP;
use Math::Prime::Util::GMP;

my($lim,$klim) = (200000,800);

for (1 .. 74) {
  for (0 .. 20000) {
    my $au = urandomb(64);
    my $bu = urandomb(64);
    my $as = negint(urandomb(63));
    my $bs = negint(urandomb(63));
    next if $bu == 0 || $bs == 0;
    tkron($au,$bu);  tkron($as,$bu);
    tkron($au,$bs);  tkron($as,$bs);
  }
  print ".";
}
print " PASS\n";


sub tkron {
  my($a,$b) = @_;

  my $k1 = kronecker($a,$b);
  my $k2 = Math::Prime::Util::GMP::kronecker($a,$b);
  my $k3 = Math::Prime::Util::PP::kronecker($a,$b);

  die "FAIL kronecker($a,$b)  MPU $k1  GMP $k2  PP $k3\n"
    unless $k1 == $k2 && $k1 == $k3;
}
