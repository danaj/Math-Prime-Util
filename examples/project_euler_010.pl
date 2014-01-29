#!/usr/bin/env perl
use warnings;
use strict;
use Math::Prime::Util qw/:all/;

my $sum = 0;
forprimes { $sum += $_ } 2_000_000;
print "$sum\n";
