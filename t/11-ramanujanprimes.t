#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/ramanujan_primes/;

my @a104272 = (2, 11, 17, 29, 41, 47, 59, 67, 71, 97, 101, 107, 127, 149, 151, 167, 179, 181, 227, 229, 233, 239, 241, 263, 269, 281, 307, 311, 347, 349, 367, 373, 401, 409, 419, 431, 433, 439, 461, 487, 491, 503, 569, 571, 587, 593, 599, 601, 607, 641, 643, 647, 653, 659, 677, 719, 727, 739, 751, 769, 809, 821, 823, 827, 853, 857, 881, 937, 941, 947, 967, 983);

my %small_range = (
  "182 to 226" => [],
  "11 to 16" => [11],
  "11 to 17" => [11,17],
  "11 to 18" => [11,17],
  "11 to 19" => [11,17],
  "11 to 20" => [11,17],
  "10 to 11" => [11],
  "11 to 29" => [11,17,29],
  "3 to 11" => [11],
  "2 to 11" => [2,11],
  "1 to 11" => [2,11],
  "0 to 11" => [2,11],
  "599 to 599" => [599],
  "10000 to 10100" => [10061,10067,10079,10091,10093],
);

plan tests => 1 + scalar(keys %small_range);

is_deeply( ramanujan_primes($a104272[-1]), \@a104272, "ramanujan_primes($a104272[-1])" );

while (my($range, $expect) = each (%small_range)) {
  my($low,$high) = $range =~ /(\d+) to (\d+)/;
  is_deeply( ramanujan_primes($low, $high), $expect, "ramanujan_primes($low,$high) should return [@{$expect}]");
}
