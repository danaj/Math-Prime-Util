#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/random_factored_integer irand factor vecprod/;

#my $use64 = (~0 > 4294967295);
#my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
#my $maxbits = $use64 ? 64 : 32;

plan tests => 4;

########

my($n, $factors) = random_factored_integer(1000000);
ok($n > 0, "random_factored_integer did not return 0");
ok($n <= 1000000, "random_factored_integer in requested range");
my @sfactors = sort {$a<=>$b} @$factors;
is_deeply( \@sfactors, [factor($n)], "factors match factor routine");
is( vecprod(@$factors), $n, "product of factors = n");
