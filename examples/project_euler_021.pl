#!/usr/bin/env perl
use warnings;
use strict;
use Math::Prime::Util qw/:all/;

my $sum = 0;
foreach my $a (1..10000) {
  my $b = divisor_sum($a)-$a;
  $sum += $a + $b if $b > $a && $a == divisor_sum($b)-$b;
}
print "$sum\n";
