#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_mersenne_prime/;
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};

my @A000043 = (2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127);
push @A000043, (512, 607, 1279, 2203, 2281) if $extra;
#push @A000043, (3217, 4253, 4423, 9689, 9941) if $extra;
#push @A000043, (11213, 19937, 21701, 23209, 44497, 86243) if $extra;

plan tests => 1;

is_deeply( [grep { is_mersenne_prime($_) } 0 .. $A000043[-1]],
           \@A000043,
           "Find Mersenne primes from 0 to $A000043[-1]" );
