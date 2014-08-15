#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/twin_primes nth_twin_prime/;

my $use64 = ~0 > 4294967295 && ~0 != 18446744073709550592;

my @small_twins = (3, 5, 11, 17, 29, 41, 59, 71, 101, 107, 137, 149, 179, 191, 197, 227, 239, 269, 281, 311, 347, 419, 431, 461, 521, 569, 599, 617, 641, 659, 809, 821, 827, 857, 881, 1019, 1031, 1049, 1061, 1091, 1151, 1229, 1277, 1289, 1301, 1319, 1427, 1451, 1481, 1487, 1607);

my %small_range = (
  "6 to 10" => [],
  "5 to 10" => [5],
  "5 to 11" => [5,11],
  "5 to 12" => [5,11],
  "5 to 13" => [5,11],
  "5 to 16" => [5,11],
  "4 to 11" => [5,11],
  "3 to 11" => [3,5,11],
  "2 to 11" => [3,5,11],
  "1 to 11" => [3,5,11],
  "0 to 11" => [3,5,11],
  "29 to 31" => [29],
  "213897 to 213997" => [213947],
  "4294957296 to 4294957796" => [4294957307,4294957397,4294957697],
  "134217228 to 134217728" => [134217401,134217437],
);

plan tests => 2 + scalar(keys %small_range);

is_deeply( twin_primes($small_twins[-1]), \@small_twins, "twin_primes($small_twins[-1])" );
{
  my @tp = map { nth_twin_prime($_) } 1 .. scalar(@small_twins);
  is_deeply( \@tp, \@small_twins, "nth_twin_prime for small values" );
}

while (my($range, $expect) = each (%small_range)) {
  my($low,$high) = $range =~ /(\d+) to (\d+)/;
  is_deeply( twin_primes($low, $high), $expect, "twin_primes($low,$high) should return [@{$expect}]");
}
