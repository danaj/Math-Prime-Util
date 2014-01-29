#!/usr/bin/env perl
use warnings;
use strict;
use Math::Prime::Util qw/:all/;

my $n = pn_primorial(4);  # Start with the first 4-factor number
# factor_exp in scalar context returns the number of distinct prime factors
$n++ while (factor_exp($n) != 4 || factor_exp($n+1) != 4 || factor_exp($n+2) != 4 || factor_exp($n+3) != 4);
print "$n\n";
