#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/almost_primes
                         almost_prime_count almost_prime_count_approx
                         almost_prime_count_lower almost_prime_count_upper
                         nth_almost_prime nth_almost_prime_approx
                         nth_almost_prime_lower nth_almost_prime_upper
                        /;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $usexs = Math::Prime::Util::prime_get_config->{'xs'};

my @small_kaps = (
  [],
  [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113],
  [4,6,9,10,14,15,21,22,25,26,33,34,35,38,39,46,49,51,55,57,58,62,65,69,74,77,82,85,86,87],
  [8,12,18,20,27,28,30,42,44,45,50,52,63,66,68,70,75,76,78,92,98,99,102,105,110,114,116,117,124,125],
  [16,24,36,40,54,56,60,81,84,88,90,100,104,126,132,135,136,140,150,152,156,184,189,196,198,204,210,220,225,228],
  [32,48,72,80,108,112,120,162,168,176,180,200,208,243,252,264,270,272,280,300,304,312,368,378,392,396,405,408,420,440],
);
my @counts_at_1e6 = (1,78498,210035,250853,198062,124465,68963,35585,17572,8491,4016,1878,865,400,179,79,35,14,7,2,0);

plan tests =>   5   # almost_primes sieve
              + 2   # count
              + 6   # nth
              + 2   # limits
              + 2   # approx
                ;

for my $k (1..5) {
  my $kap = $small_kaps[$k];
  is_deeply(almost_primes($k,1,$kap->[-1]), $kap, "small $k-almost-primes");
}

# TODO some special cases for k=0

# TODO: more testing.  Perhaps with extra or in xt

###### Test almost_prime_count
is_deeply([map { almost_prime_count($_, 1e6) } 0..20], \@counts_at_1e6, "k-almost prime counts at 1000000 for k=1..20");

is(almost_prime_count(17,1e9), 38537, "There are 38537 17-almost-primes <= 1,000,000,000");

###### Test nth_almost_prime
is(nth_almost_prime(1,2), 3, "2nd 1-almost-prime is 3");
is(nth_almost_prime(2,34), 95, "34th 2-almost-prime is 94");
is(nth_almost_prime(3,456), 1802, "456th 3-almost-prime is 1802");
is(nth_almost_prime(4,5678), 31382, "5678th 4-almost-prime is 31382");
is(nth_almost_prime(5,67890), 558246, "67890th 5-almost-prime is 558246");
is(nth_almost_prime(24,5555), "21678243840", "5555th 24-almost-prime is 21678243840");

###### Test limits
is( cmp_kap(3,59643,234618), 234618, "3-almost prime limits for 59643" );
is( cmp_kap(32,12,"26843545600"), "26843545600", "32-almost prime limits for 12" );

###### Test approx
is( approx_in_range(3,59643,234618), 234618, "approx 3-almost prime");
is( approx_in_range(32,12,"26843545600"), "26843545600", "32-almost prime limits for 12" );


sub cmp_closeto {
  my $got = shift;
  my $expect = shift;
  my $tolerance = shift;
  my $message = shift;
  cmp_ok( abs($got - $expect), '<=', $tolerance, $message );
}

sub cmp_kap {
  my($k,$n,$rn) = @_;
  return 'nth lower' unless nth_almost_prime_lower($k,$n) <= $rn;
  return 'nth upper' unless nth_almost_prime_upper($k,$n) >= $rn;
  return 'pc lower' unless almost_prime_count_lower($k,$rn) <= $n;
  return 'pc upper' unless almost_prime_count_upper($k,$rn) >= $n;
  $rn;
}

sub approx_in_range {
  my($k,$n,$rn) = @_;
  my $arn = nth_almost_prime_approx($k,$n);
  my $an  = almost_prime_count_approx($k,$rn);
  return 'nth approx too low' if $arn < ($rn-$rn/20);
  return 'nth approx too high' if $arn > ($rn+$rn/20);
  return 'count approx too low' if $an < ($n-$n/20);
  return 'count approx too high' if $an > ($n+$n/20);
  $rn;
}
