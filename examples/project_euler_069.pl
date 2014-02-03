#!/usr/bin/env perl
use warnings;
use strict;
use Math::Prime::Util qw/euler_phi pn_primorial/;

# Better way
my $n = 0;
$n++ while pn_primorial($n+1) < 1000000;
print pn_primorial($n), "\n";

# Brute force
my ($maxn, $maxratio, $ratio) = (0, 0);
foreach my $n (1 .. 1000000) {
  $ratio = $n / euler_phi($n);
  ($maxn, $maxratio) = ($n, $ratio) if $ratio > $maxratio;
}
print "$maxn  $maxratio\n";
