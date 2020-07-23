#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_practical/;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};

my @a005153 = (1, 2, 4, 6, 8, 12, 16, 18, 20, 24, 28, 30, 32, 36, 40, 42, 48, 54, 56, 60, 64, 66, 72, 78, 80, 84, 88, 90, 96, 100, 104, 108, 112, 120, 126, 128, 132, 140, 144, 150, 156, 160, 162, 168, 176, 180, 192, 196, 198, 200, 204, 208, 210, 216, 220, 224, 228, 234, 240, 252);
my @a174533 = (70, 350, 490, 770, 910, 945, 1190, 1330, 1575, 1610, 1750, 2030, 2170, 2205, 2450, 2584, 2590, 2835, 2870, 3010, 3128, 3290, 3430, 3465, 3710, 3850, 3944, 4095, 4130, 4216, 4270, 4550, 4690, 4725, 5355, 5390, 5775, 5950, 5985, 6370, 6615, 6650, 6825);

plan tests => 3;

is_deeply( [grep { is_practical($_) } 0..$a005153[-1]], \@a005153, "is_practical(0 .. $a005153[-1])");

is(is_practical(429606), 1, "is_practical(429606) = 1");

is_deeply( [map { is_practical($_) } @a174533], [map { 0 } @a174533], "is_practical(n) = 0 for almost practical numbers");
