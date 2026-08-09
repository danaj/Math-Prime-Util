#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_higgs_prime/;

my @higgs200 = qw(
  2 3 5 7 11 13 19 23 29 31 37 43 47 53 59 61 67 71 79 101 107 127
  131 139 149 151 157 173 181 191 197 199
);

is_deeply([grep { is_higgs_prime($_) } -5 .. 200], \@higgs200,
          "is_higgs_prime(-5..200)");

my @large_higgs = (
  '1695341',              # p-1 = 2^2 * 5 * 29 * 37 * 79
  '1122724972381',        # p-1 = 2^2 * 3 * 5 * 11^2 * 19 * 29 * 59 * 67 * 71
  '50465167359954207391301',
  # p-1 = 2^2 * 5^2 * 7 * 11^2 * 19^2 * 23 * 43^2 * 47 * 53^2 * 61^2 * 79
);
is_deeply([grep { is_higgs_prime($_) } @large_higgs], \@large_higgs,
          "large Higgs primes");

done_testing;
