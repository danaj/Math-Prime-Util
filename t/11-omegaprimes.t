#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/omega_primes omega_prime_count nth_omega_prime/;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $usexs = Math::Prime::Util::prime_get_config->{'xs'};

my @small_kops = (
  [],
  [2,3,4,5,7,8,9,11,13,16,17,19,23,25,27,29,31,32,37,41,43,47,49,53,59,61,64,67,71,73,79,81,83,89,97,101,103,107,109,113],
  [6,10,12,14,15,18,20,21,22,24,26,28,33,34,35,36,38,39,40,44,45,46,48,50,51,52,54,55,56,57,58,62,63,65,68,69,72,74,75,76],
  [30,42,60,66,70,78,84,90,102,105,110,114,120,126,130,132,138,140,150,154,156,165,168,170,174,180,182,186,190,195,198,204,220,222,228,230,231,234,238,240],
  [210,330,390,420,462,510,546,570,630,660,690,714,770,780,798,840,858,870,910,924,930,966,990,1020,1050,1092,1110,1122,1140,1155,1170,1190,1218,1230,1254,1260,1290,1302,1320,1326],
  [2310,2730,3570,3990,4290,4620,4830,5460,5610,6006,6090,6270,6510,6630,6930,7140,7410,7590,7770,7854,7980,8190,8580,8610,8778,8970,9030,9240,9282,9570,9660,9690,9870,10010,10230,10374,10626,10710,10920,11130],
);
my @counts_at_1e6 = (1,78734,288726,379720,208034,42492,2285,8,0,0,0,0,0,0,0,0,0,0,0,0,0);
my @counts_at_1e4 = (1,1280,4097,3695,894,33,0,0,0,0,0);

plan tests =>   5   # omega_primes sieve
              + 6   # count
              + 8   # nth_omega_prime
                ;

###### omega_primes

for my $k (1..5) {
  my $kop = $small_kops[$k];
  is_deeply(omega_primes($k,$kop->[-1]), $kop, "small $k-omega-primes");
}

###### omega_prime_count

is_deeply([map{omega_prime_count($_,206)}1..10],[60,113,32,0,0,0,0,0,0,0],"omega_prime_count n=206, k 1..10");

is_deeply([map { omega_prime_count($_, 1e4) } 0..10], \@counts_at_1e4, "k-omega prime counts at 10000 for k=1..10");

is(omega_prime_count(6,9e4), 19, "There are 19 6-omega-primes <= 90,000");
is(omega_prime_count(8,2e7), 10, "There are 10 8-omega-primes <= 20,000,000");

SKIP: {
  skip "Slow in PP", 2 unless $usexs || $extra;
  is_deeply([map { omega_prime_count($_, 1e6) } 0..20], \@counts_at_1e6, "k-omega prime counts at 1000000 for k=1..20");
  is_deeply([map{omega_prime_count($_,206111)}1..10],[18613,66025,79801,36364,5182,125,0,0,0,0],"omega_prime_count n=206111, k 1..10");
}

###### nth_omega_prime

is_deeply([map { nth_omega_prime(0,$_) } 0..4], [undef,1,undef,undef,undef],
          "nth_omega_prime(0,...)");
is_deeply([map { nth_omega_prime($_,0) } 1..10], [(undef) x 10],
          "nth_omega_prime(...,0)");
for my $k (1..5) {
  my $nops = scalar(@{$small_kops[$k]});
  $nops = 10 unless $usexs;
  my @expect = @{$small_kops[$k]}[0..$nops-1];
  is_deeply([map { nth_omega_prime($k,$_) } 1..$nops], \@expect,
          "nth_omega_prime($k, 1 .. $nops)");
}
SKIP: {
  skip "nth_omega_prime is very slow in PP", 1 unless $usexs;
  is(nth_omega_prime(8,122), 46692030, "The 122nd 8-omega prime is 46692030");
}
