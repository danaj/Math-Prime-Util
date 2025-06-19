#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/primes nth_prime nth_twin_prime
                         nth_prime_lower nth_prime_upper
                         nth_prime_approx nth_twin_prime_approx
                         nth_semiprime is_semiprime
                         inverse_li inverse_li_nv/;

my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};

my $nsmallprimes = 1000;
my $nth_small_prime = 7919;  # nth_prime(1000)

my %pivals32 = (
                  1 => 0,
                 10 => 4,
                100 => 25,
               1000 => 168,
              10000 => 1229,
             100000 => 9592,
            1000000 => 78498,
);

#  Powers of 10: http://oeis.org/A006988/b006988.txt
#  Powers of  2: http://oeis.org/A033844/b033844.txt
my %nthprimes32 = (
                  1 => 2,
                 10 => 29,
                100 => 541,
               1000 => 7919,
              10000 => 104729,
             100000 => 1299709,
            1000000 => 15485863,
           10000000 => 179424673,
          100000000 => 2038074743,
  # Some values that estimate right around the value
            6305537 => 110040379,
            6305538 => 110040383,
            6305539 => 110040391,
            6305540 => 110040407,
            6305541 => 110040467,
            6305542 => 110040499,
            6305543 => 110040503,
);
my %nthprimes64 = (
         1000000000 => 22801763489,
        10000000000 => 252097800623,
       100000000000 => 2760727302517,
      1000000000000 => 29996224275833,
     10000000000000 => 323780508946331,
    100000000000000 => 3475385758524527,
#   1000000000000000 => 37124508045065437,
#  10000000000000000 => 394906913903735329,
# 100000000000000000 => 4185296581467695669,
);
if ($usexs || $extra) {
  $nthprimes64{  "1000000000000000"} =   "37124508045065437";
  $nthprimes64{ "10000000000000000"} =  "394906913903735329";
  $nthprimes64{"100000000000000000"} = "4185296581467695669";
}
my %nthprimes_small = map { $_ => $nthprimes32{$_} }
                      grep { ($_ <= 10_000_000) || $extra }
                      keys %nthprimes32;

my @small_primes = (undef, @{primes($nth_small_prime)});

my %ntpcs = (
             5 =>                   29,
            50 =>                 1487,
           500 =>                32411,
          5000 =>               557519,
         50000 =>              8264957,
        500000 =>            115438667,
       5000000 =>           1523975909,
      50000000 =>          19358093939,
     500000000 =>         239211160649,
);

my %nthsemi = (
          1234 =>          4497,
         12345 =>         51019,
        123456 =>        573355,
);
$nthsemi{1234567}     =      6365389  if $usexs || $extra;
$nthsemi{12345678}    =     69914722  if $usexs || $extra;
$nthsemi{123456789}   =    760797011  if $usexs && $extra;
$nthsemi{1234567890}  =   8214915893  if $usexs && $extra && $use64;
$nthsemi{8589934592}  =  60662588879  if $usexs && $extra && $use64;
$nthsemi{17179869184} = 123806899739  if $usexs && $extra && $use64;

plan tests => 1   # nth_prime
            + 1   # nth_prime lower/upper/approx
            + 1   # nth_twin_prime and approx
            + 1   # nth_semiprime
            + 1   # inverse_li and inverse_li_nv
            + 0;

subtest 'nth_prime', sub {

  is(nth_prime(0), undef, "nth_prime(0) = undef");

  my @pri100 = @{primes(541)};  # First 100 primes
  my @nth100 = map { nth_prime($_) } 1 .. 100;
  is_deeply(\@nth100, \@pri100, "nth_prime(1..100)");

  is_deeply([map { nth_prime($_) } 1229,1230,9592,9593,78498,78499],
            [9973,10007,99991,100003,999983,1000003],
            "nth_prime(n) results around 10k, 100k, 1M");

  is_deeply([map { nth_prime($_) } 6305537, 6305540, 6305543],
            [110040379,110040407,110040503],
            "nth_prime(6305537, 6305540, 6305543)");

  is_deeply([map { nth_prime($_) } 1000,10000,100000,1000000,10000000],
            [7919,104729,1299709,15485863,179424673],
            "nth_prime(n) 1k,10k,100k,1M,10M");

  # Test an nth prime value that uses the binary-search-on-R(n) algorithm
  SKIP: {
    # TODO: Can we use the NV R to solve this with PP?
    skip "nth_prime(21234567890) for 64-bit XS EXTENDED_TESTING", 1 unless $extra && $use64 && $usexs;
    is( nth_prime(21234567890), 551990503367, "nth_prime(21234567890)" );
  }
};

subtest 'nth_prime upper/lower/approx', sub {
  while (my($n, $nth) = each (%nthprimes32)) {
    cmp_ok(nth_prime_upper($n),'>=',$nth, "nth_prime($n) <= upper estimate");
    cmp_ok(nth_prime_lower($n),'<=',$nth, "nth_prime($n) >= lower estimate");

    my $okpct = ($n >= 775) ? 1 : 2;
    cmp_within(nth_prime_approx($n), $nth, $okpct, "nth_prime_approx($n)");
  }
  if ($use64) {
    while (my($n, $nth) = each (%nthprimes64)) {
      cmp_ok(nth_prime_upper($n),'>=',$nth, "nth_prime($n) <= upper estimate");
      cmp_ok(nth_prime_lower($n),'<=',$nth, "nth_prime($n) >= lower estimate");
      cmp_within(nth_prime_approx($n), $nth, .001, "nth_prime_approx($n)");
    }
  }

  my $maxindex   = $use64 ? '425656284035217743' : '203280221';
  my $maxindexp1 = $use64 ? '425656284035217744' : '203280222';
  my $maxprime   = $use64 ? '18446744073709551557' : '4294967291';
  cmp_ok(nth_prime_lower($maxindex),'<=',$maxprime, "nth_prime_lower(maxindex) <= maxprime");
  cmp_ok(nth_prime_upper($maxindex),'>=',$maxprime, "nth_prime_upper(maxindex) >= maxprime");
  cmp_ok(nth_prime_lower($maxindexp1),'>=',nth_prime_lower($maxindex), "nth_prime_lower(maxindex+1) >= nth_prime_lower(maxindex)");
};

####################################
subtest 'nth_twin_prime and approx', sub {
  is( nth_twin_prime(0), undef, "nth_twin_prime(0) = undef" );
  is( nth_twin_prime(17), 239, "239 = 17th twin prime" );
  is( nth_twin_prime(1234), 101207, "101207 = 1234'th twin prime" );

  while (my($n, $nthtpc) = each (%ntpcs)) {
    cmp_within(nth_twin_prime_approx($n), $nthtpc, 2, "nth_twin_prime_approx($n)");
  }
};

####################################
subtest 'nth_semiprime', sub {

  is( nth_semiprime(0), undef, "nth_semiprime(0) = undef" );

  my $range = $extra ? 10000 : 500;
  my @semiprimes = grep { is_semiprime($_) } 0 .. $range;
  my $nsmall = scalar(@semiprimes);
  my @nth_semis = map { nth_semiprime($_) } 1 .. $nsmall;
  is_deeply(\@nth_semis, \@semiprimes, "nth_semiprime(1 .. $nsmall)");

  while (my($n, $nthsemi) = each (%nthsemi)) {
    is( nth_semiprime($n), $nthsemi, "nth_semiprime($n) = $nthsemi" );
  }
};

####################################
subtest 'inverse_li and inverse_li_nv', sub {
  is_deeply(
    [ map { inverse_li($_) } 0 .. 50 ],
    [qw/0 2 3 5 6 8 10 12 15 18 21 24 27 30 34 37 41 45 49 53 57 61 65 69 73 78 82 86 91 95 100 105 109 114 119 123 128 133 138 143 148 153 158 163 168 173 179 184 189 194 199/],
    "inverse_li: Li^-1(0..50)"
  );
  # Allow +/- 2 for floating point differences in LogarithmicIntegral
  like(inverse_li(1000000000), qr/^2280162741[34567]$/, "inverse_li(1e9)");
  like(inverse_li(1100000000000), qr/^3310443690704[01234]$/, "inverse_li(11e11)");

  cmp_within(inverse_li_nv(4), 5.60927669305089, .001, "inverse_li_nv(4)");
  cmp_within(inverse_li_nv(64.2731216921018), 277, .001, "inverse_li_nv(64.2731216921018)");
  cmp_within(inverse_li_nv(40000), 478956.000953764, .001, "inverse_li_nv(40000)");
  cmp_within(inverse_li_nv(1234567890123), 37301814610592.3, .1, "inverse_li_nv(1234567890123)");
};

sub cmp_within {
  my($got, $exp, $tolpct, $mess) = @_;

  my $diffok = 0.01 * $tolpct * abs("$exp") + 0.5;
  my($diff)  = map { 0.0+"$_" } ($got > $exp ? $got-$exp : $exp-$got);
  my $diffstr = ($got == $exp) ? "0" : ($got > $exp) ? "+$diff" : "-$diff";

  ok($diff <= $diffok, "$mess =~ $exp ($diffstr)")
    or diag "$mess got $got with diff $diff > $diffok";
}
