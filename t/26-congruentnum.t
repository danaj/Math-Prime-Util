#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_congruent_number/;
#my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
#my $usegmp = Math::Prime::Util::prime_get_config->{'gmp'};
#my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};

my @cn200 = (5,6,7,13,14,15,20,21,22,23,24,28,29,30,31,34,37,38,39,41,45,46,47,52,53,54,55,56,60,61,62,63,65,69,70,71,77,78,79,80,84,85,86,87,88,92,93,94,95,96,101,102,103,109,110,111,112,116,117,118,119,120,124,125,126,127,133,134,135,136,137,138,141,142,143,145,148,149,150,151,152,154,156,157,158,159,161,164,165,166,167,173,174,175,180,181,182,183,184,188,189,190,191,194,197,198,199);

my @cn1e6 = (1,5,6,7,9);

plan tests => 0
            + 2
            ;

# This covers all current special cases in the code
is_deeply([grep { is_congruent_number($_) } 1..200], \@cn200, "congruent numbers to 200");

is_deeply([grep { is_congruent_number(1000000+$_) } 1..10], \@cn1e6, "congruent numbers 10^6 + (1..10)");
