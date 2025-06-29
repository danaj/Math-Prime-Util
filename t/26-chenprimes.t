#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/ is_chen_prime next_chen_prime /;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};

my @chen200 = (2,3,5,7,11,13,17,19,23,29,31,37,41,47,53,59,67,71,83,89,101,107,109,113,127,131,137,139,149,157,167,179,181,191,197,199);

plan tests => 4;

is_deeply( [grep { is_chen_prime($_) } 0..200],
           \@chen200,
           "is_chen_prime(0..200)" );

is_deeply( [map { next_chen_prime($chen200[$_]) } 0..$#chen200-1],
           [@chen200[1..$#chen200]],
           "next_chen_prime for small values" );

is(is_chen_prime("10000000000687"), 1, "is_chen_prime(10^13+687)");

is("".next_chen_prime("10000000000000000"), "10000000000000079", "next_chen_prime(10^16)");
