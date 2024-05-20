#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_perfect_number/;

#my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
#my $usexs  = Math::Prime::Util::prime_get_config->{'xs'};
#my $usegmp = Math::Prime::Util::prime_get_config->{'gmp'};

my @true  = (qw/8128 33550336 8589869056 137438691328 2305843008139952128 2658455991569831744654692615953842176/);
my @false = (qw/8505 12285 19845 28665 31185 198585576189 8 32 2096128 35184367894528 144115187807420416 9444732965670570950656/);


plan tests => 1 + 2;
            ;

is_deeply([grep { is_perfect_number($_) } -10 .. 500], [6,28,496], "is_perfect_number(-10 .. 500)");

is_deeply([map { is_perfect_number($_) } @true], [map { 1 } @true], "perfect: [@true]");
is_deeply([map { is_perfect_number($_) } @false], [map { 0 } @false], "not perfect: [@false]");
