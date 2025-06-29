#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/pisano_period factorial/;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
#my $usexs  = Math::Prime::Util::prime_get_config->{'xs'};
#my $usegmp = Math::Prime::Util::prime_get_config->{'gmp'};

my @tests = (
  [28657, 92],  # F_23
  [64079, 46],  # L_23
  [3590807, 3264380],
  [3628800, 86400],
  ["2980232238769531250","17881393432617187500"],
  ["14901161193847656250","89406967163085937500"],
  ["74505805969238281250","447034835815429687500"],
);

my @pisano = (qw/0 1 3 8 6 20 24 16 12 24 60 10 24 28 48 40 24 36 24 18 60 16 30 48 24 100 84 72 48 14 120 30 48 40 36 80 24 76 18 56 60 40 48 88 30 120 48 32 24 112 300 72 84 108 72 20 48 72 42 58 120 60 30 48 96 140 120 136 36 48 240 70 24 148 228 200 18 80 168 78 120 216 120 168 48 180 264 56 60 44 120 112 48 120 96 180 48 196 336 120 300 50 72 208 84 80 108 72 72 108 60 152 48 76 72 240 42 168 174 144 120 110 60 40 30 500 48 256 192 88 420 130 120 144 408 360 36 276 48 46 240 32 210 140 24 140 444 112 228 148 600 50 36 72 240 60 168 316 78 216 240 48 216 328 120 40 168 336 48 364 180 72 264 348 168 400 120 232 132 178 120/);

plan tests => 1
            + scalar(@tests)
            + 1
            ;

is_deeply([map {pisano_period($_)} 0..180], \@pisano, "pisano_period(0..180)");

for my $data (@tests) {
  my($n,$exp) = @$data;
  is("".pisano_period($n), $exp, "pisano_period($n) = $exp");
}

is("".pisano_period(factorial(30)), "204996473853050880000000", "pisano_period(30!)");
