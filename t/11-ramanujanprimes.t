#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/ramanujan_primes nth_ramanujan_prime is_ramanujan_prime
                         nth_ramanujan_prime_upper nth_ramanujan_prime_lower
                         nth_ramanujan_prime_approx
                         ramanujan_prime_count_upper ramanujan_prime_count_lower
                         ramanujan_prime_count_approx/;

my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};

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

plan tests => 1 + scalar(keys %small_range) + 2 + 1 + 2 + 3 + 2;

is_deeply( ramanujan_primes($a104272[-1]), \@a104272, "ramanujan_primes($a104272[-1])" );

while (my($range, $expect) = each (%small_range)) {
  my($low,$high) = $range =~ /(\d+) to (\d+)/;
  is_deeply( ramanujan_primes($low, $high), $expect, "ramanujan_primes($low,$high) should return [@{$expect}]");
}

my @smalla = grep { $_ < ($usexs ? 1000 : 500) } @a104272;

{
  my @rp;
  push @rp, nth_ramanujan_prime($_) for 1 .. scalar(@smalla);
  is_deeply( \@rp, \@smalla, "nth_ramanujan_prime(1 .. ".scalar(@smalla).")");
  if ($usexs) {
    is( nth_ramanujan_prime(123456), 3657037, "The 123,456th Ramanujan prime is 3657037" );
  } else {
    is( nth_ramanujan_prime(1234), 24043, "The 1,234th Ramanujan prime is 24043" );
  }
}

{
  my @rp;
  for (0 .. $smalla[-1]) {
    push @rp, $_ if is_ramanujan_prime($_);
  }
  is_deeply( \@rp, \@smalla, "is_ramanujan_prime( 0 .. ".scalar(@smalla).")");
}

is(nth_ramanujan_prime(997), 19379, "997th Ramanujan prime is 19379");
SKIP: {
  skip "Without XS, Ramanujan primes are slow",1 unless $usexs || $extra;
  is(nth_ramanujan_prime(23744), 617759, "Rn[23744] is 617759");
}

is_deeply( [map{cmp_rn($_+1,$a104272[$_])} 0..$#a104272], \@a104272, "small ramanujan prime limits" );
is( cmp_rn(59643,1673993), 1673993, "ramanujan prime limits for 59643" );
is( cmp_rn(5964377,225792607), 225792607, "ramanujan prime limits for 5964377" );

is( approx_in_range(59643,1673993), 1673993, "ramanujan prime approx for 59643" );
is( approx_in_range(5964377,225792607), 225792607, "ramanujan prime approx for 5964377" );

sub cmp_rn {
  my($n,$rn) = @_;
  return 'nth lower' unless nth_ramanujan_prime_lower($n) <= $rn;
  return 'nth upper' unless nth_ramanujan_prime_upper($n) >= $rn;
  return 'pc lower' unless ramanujan_prime_count_lower($rn) <= $n;
  return 'pc upper' unless ramanujan_prime_count_upper($rn) >= $n;
  $rn;
}

sub approx_in_range {
  my($n,$rn) = @_;
  my $arn = nth_ramanujan_prime_approx($n);
  my $an  = ramanujan_prime_count_approx($rn);
  return 'nth approx too low' if "$arn" < ($rn-$rn/50);
  return 'nth approx too high' if "$arn" > ($rn+$rn/50);
  return 'count approx too low' if "$an" < ($n-$n/50);
  return 'count approx too high' if "$an" > ($n+$n/50);
  $rn;
}
