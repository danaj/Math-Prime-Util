#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/primes nth_prime nth_twin_prime
                         nth_prime_lower nth_prime_upper
                         nth_prime_approx nth_twin_prime_approx/;

my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $broken64 = (18446744073709550592 == ~0);

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
   1000000000000000 => 37124508045065437,
  10000000000000000 => 394906913903735329,
 100000000000000000 => 4185296581467695669,
);
my %nthprimes_small = map { $_ => $nthprimes32{$_} }
                      grep { ($_ <= 10_000_000) || $extra }
                      keys %nthprimes32;

my @small_primes = (0, @{primes($nth_small_prime)});

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


plan tests => 0 + 2*scalar(keys %pivals32)
                + 1
                + 3*scalar(keys %nthprimes32)
                + scalar(keys %nthprimes_small)
                + $use64 * 3 * scalar(keys %nthprimes64)
                + 3   # nth_prime_lower with max index
                + 3   # nth_twin_prime
                + scalar(keys %ntpcs)   # nth_twin_prime_approx
                + (($extra && $use64 && $usexs) ? 1 : 0);


while (my($n, $pin) = each (%pivals32)) {
  my $next = $pin+1;
  cmp_ok( nth_prime($pin), '<=', $n, "nth_prime($pin) <= $n");
  cmp_ok( nth_prime($next), '>=', $n, "nth_prime($next) >= $n");
}

{
  my @nth_primes = map { nth_prime($_) }  (0 .. $#small_primes);
  is_deeply( \@nth_primes, \@small_primes, "nth_prime for primes 0 .. $#small_primes" );
}

while (my($n, $nth) = each (%nthprimes32)) {
  cmp_ok( nth_prime_upper($n), '>=', $nth, "nth_prime($n) <= upper estimate" );
  cmp_ok( nth_prime_lower($n), '<=', $nth, "nth_prime($n) >= lower estimate" );

  my $approx = nth_prime_approx($n);
  my $percent_limit = ($n >= 775) ? 1 : 2;
  cmp_ok( abs($nth - $approx) / $nth, '<=', $percent_limit/100.0, "nth_prime_approx($n) = $approx within $percent_limit\% of $nth");
}
while (my($n, $nth) = each (%nthprimes_small)) {
  is( nth_prime($n), $nth, "nth_prime($n) = $nth" );
}

if ($use64) {
  while (my($n, $nth) = each (%nthprimes64)) {
    cmp_ok( nth_prime_upper($n), '>=', $nth, "nth_prime($n) <= upper estimate" );
    cmp_ok( nth_prime_lower($n), '<=', $nth, "nth_prime($n) >= lower estimate" );

    my $approx = nth_prime_approx($n);
    my $percent_limit = 0.001;
    cmp_ok( abs($nth - $approx) / $nth, '<=', $percent_limit/100.0, "nth_prime_approx($n) = $approx within $percent_limit\% of $nth");
  }
}

my $maxindex   = $use64 ? '425656284035217743' : '203280221';
my $maxindexp1 = $use64 ? '425656284035217744' : '203280222';
my $maxprime   = $use64 ? '18446744073709551557' : '4294967291';
cmp_ok( nth_prime_lower($maxindex), '<=', $maxprime, "nth_prime_lower(maxindex) <= maxprime");
cmp_ok( nth_prime_upper($maxindex), '>=', $maxprime, "nth_prime_upper(maxindex) >= maxprime");
cmp_ok( nth_prime_lower($maxindexp1), '>=', nth_prime_lower($maxindex), "nth_prime_lower(maxindex+1) >= nth_prime_lower(maxindex)");

my $overindex = ($broken64) ? 425656284035217843 : $maxindexp1;

if ($extra && $use64 && $usexs) {
  # Test an nth prime value that uses the binary-search-on-R(n) algorithm
  is( nth_prime(21234567890), 551990503367, "nth_prime(21234567890)" );
}

####################################3

is( nth_twin_prime(0), 0, "nth_twin_prime(0) = 0" );
is( nth_twin_prime(17), 239, "239 = 17th twin prime" );
is( nth_twin_prime(1234), 101207, "101207 = 1234'th twin prime" );

while (my($n, $nthtpc) = each (%ntpcs)) {
  my $approx = nth_twin_prime_approx($n);
  my $errorp = 100 * abs($nthtpc - $approx) / $nthtpc;
  my $estr = sprintf "%8.6f%%", $errorp;
  cmp_ok( $errorp, '<=', 2, "nth_twin_prime_approx($n) is $estr (got $approx, expected ~$nthtpc)");
}
