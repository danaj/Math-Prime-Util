#!/usr/bin/env perl
use warnings;
use strict;
use Math::Prime::Util qw/euler_phi/;

my $sum = 0;
$sum += $_ for euler_phi(2,1000000);
print "$sum\n";
