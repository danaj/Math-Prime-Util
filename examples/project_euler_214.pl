#!/usr/bin/env perl
use warnings;
use strict;
use Math::Prime::Util qw/forprimes euler_phi/;

my $limit = shift || 40000000;
my $cl = shift || 25;

my @c;
sub totchainlen {
  my $n = shift;
  return $n if $n <= 2;
  $c[$n] //= 1 + totchainlen(euler_phi($n));
  return $c[$n];
}

my $sum = 0;
forprimes {
  $sum += $_ if totchainlen($_) == $cl;
} $limit;
print "$sum\n";
