#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/semi_primes
                         semiprime_count semiprime_count_approx
                         nth_semiprime nth_semiprime_approx/;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $usexs = Math::Prime::Util::prime_get_config->{'xs'};

my @small_semis = (4,6,9,10,14,15,21,22,25,26,33,34,35,38,39,46,49,51,55,57,58,62,65,69,74,77,82,85,86,87,91,93,94,95);

my %small_range = (
  "11 to 13" => [],
  "10 to 10" => [10],
  "10 to 11" => [10],
  "10 to 12" => [10],
  "10 to 13" => [10],
  "10 to 14" => [10,14],
  "5 to 16" => [6,9,10,14,15],
  "4 to 11" => [4,6,9,10],
  "3 to 11" => [4,6,9,10],
  "2 to 11" => [4,6,9,10],
  "1 to 11" => [4,6,9,10],
  "0 to 11" => [4,6,9,10],
  "26 to 33" => [26,33],
  "25 to 34" => [25,26,33,34],
  "184279943 to 184280038" => [184279943,184279969,184280038],
  "184279944 to 184280037" => [184279969],
  "8589990147 to 8589990167" => [8589990149,8589990157,8589990166],
);

my %small_semis = (
    1234 =>   4497,
   12345 =>  51019,
  123456 => 573355,
);

my %big_semis = (
          "2147483648" =>         "14540737711",
       "4398046511104" =>      "36676111297003",
  "100000000000000000" => "1030179406403917981",
  "288230376151711744" => "3027432768282284351",
);

my %small_counts = (
    1234 =>   363,
   12345 =>  3217,
  123456 => 28589,
);
my %range_counts = (
  "1000000 to 1000100" => 25,
  "1000000000 to 1000000100" =>     14,   # iterate
  "1000000000 to 1000010000" =>   1567,   # sieve
  "1000000000 to 1001000000" => 155612,   # count
);

my %big_counts = (
             "100000000" =>           "17427258",
          "100000000000" =>        "13959990342",
       "100000000000000" =>     "11715902308080",
  "10000000000000000000" => "932300026230174178",
);

plan tests => 2
            + scalar(keys %small_range)
            + scalar(keys %small_semis)
            + scalar(keys %small_counts)
            + scalar(keys %range_counts)
            + scalar(keys %big_counts)
            + scalar(keys %big_semis)
            + 1 * $extra;

is_deeply( semi_primes($small_semis[-1]), \@small_semis, "semi_primes($small_semis[-1])" );

{
  my @tp = map { nth_semiprime($_) } 1 .. scalar(@small_semis);
  is_deeply( \@tp, \@small_semis, "nth_semiprime for small values" );
}

while (my($range, $expect) = each (%small_range)) {
  my($low,$high) = $range =~ /(\d+) to (\d+)/;
  is_deeply( semi_primes($low, $high), $expect, "semi_primes($low,$high) should return [@{$expect}]");
}

while (my($n, $spc) = each (%small_counts)) {
  is( semiprime_count($n), $spc, "semiprime_count($n) = $spc");
}
while (my($range, $count) = each (%range_counts)) {
  my($low,$high) = $range =~ /(\d+) to (\d+)/;
  is_deeply( semiprime_count($low, $high), $count, "semiprime_count($low,$high) = $count");
}

while (my($n, $nth) = each (%small_semis)) {
  SKIP: {
    skip "PP nth_semiprime is slow",1 unless $n < 10000 || $usexs || $extra;
    is( nth_semiprime($n), $nth, "nth_semiprime($n) = $nth");
  }
}


while (my($n, $spc) = each (%big_counts)) {
  # XS routine is within 0.00001.  PP within 0.002.
  cmp_closeto( semiprime_count_approx($n), $spc, 0.002 * abs($spc), "semiprime_count_approx($n) ~ $spc");
}

while (my($n, $nth) = each (%big_semis)) {
  # XS routine is within 0.00001.  PP within 0.001.
  cmp_closeto( nth_semiprime_approx($n), $nth, 0.001 * abs($nth), "nth_semiprime_approx($n) ~ $nth");
}

if ($extra && $usexs) {
  # More than one interpolation needed
  is( nth_semiprime(12479400000), 89102597117, "nth_semiprime(12479400000) = 89102597117" );
}

sub cmp_closeto {
  my $got = shift;
  my $expect = shift;
  my $tolerance = shift;
  my $message = shift;
  cmp_ok( abs($got - $expect), '<=', $tolerance, $message );
}
