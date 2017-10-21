#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_fundamental/;
#my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
#my $usegmp = Math::Prime::Util::prime_get_config->{'gmp'};
#my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};

plan tests => 0
            + 4
            ;

is_deeply( [grep { is_fundamental($_) } -50 .. 0],
           [-47,-43,-40,-39,-35,-31,-24,-23,-20,-19,-15,-11,-8,-7,-4,-3],
           "is_fundamental(-50 .. 0)" );
is_deeply( [grep { is_fundamental($_) } 0 .. 50],
           [1,5,8,12,13,17,21,24,28,29,33,37,40,41,44],
           "is_fundamental(0 .. 50)" );

is( is_fundamental("147573952589676412929"), 1, "is_fundamental(2^67+1)" );
is( is_fundamental("-147573952589676412920"), 1, "is_fundamental(-2^67+8)" );
