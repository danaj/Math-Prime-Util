#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/primes nth_prime nth_prime_lower nth_prime_upper nth_prime_approx/;

my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
my $extra = defined $ENV{RELEASE_TESTING} && $ENV{RELEASE_TESTING};
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
                      grep { ($_ <= 2000000) || $extra }
                      keys %nthprimes32;

my @small_primes = (0, @{primes($nth_small_prime)});

#plan tests => 7*2 + $nsmallprimes+1 + 9*3 + 7 + ($extra ? 9 : 7) + ($use64 ? 9*3 : 0);

plan tests => 0 + 2*scalar(keys %pivals32)
                + scalar @small_primes
                + 3*scalar(keys %nthprimes32)
                + scalar(keys %nthprimes_small)
                + $use64 * 3 * scalar(keys %nthprimes64)
                + 4
                + 3;


while (my($n, $pin) = each (%pivals32)) {
  my $next = $pin+1;
  cmp_ok( nth_prime($pin), '<=', $n, "nth_prime($pin) <= $n");
  cmp_ok( nth_prime($next), '>=', $n, "nth_prime($next) >= $n");
}

foreach my $n (0 .. $#small_primes) {
  is(nth_prime($n), $small_primes[$n], "The ${n}th prime is $small_primes[$n]");
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

my $maxindex = $use64 ? 425656284035217743 : 203280221;
my $maxprime = $use64 ? 18446744073709551557 : 4294967291;
cmp_ok( nth_prime_lower($maxindex), '<=', $maxprime, "nth_prime_lower(maxindex) <= maxprime");
cmp_ok( nth_prime_upper($maxindex), '>=', $maxprime, "nth_prime_upper(maxindex) >= maxprime");
cmp_ok( nth_prime_approx($maxindex), '==', $maxprime, "nth_prime_approx(maxindex) == maxprime");
cmp_ok( nth_prime_lower($maxindex+1), '>=', nth_prime_lower($maxindex), "nth_prime_lower(maxindex+1) >= nth_prime_lower(maxindex)");

my $add = ($broken64) ? 100 : 1;

eval { nth_prime_upper($maxindex+$add); };
like($@, qr/overflow/, "nth_prime_upper(maxindex+$add) overflows");

eval { nth_prime_approx($maxindex+$add); };
like($@, qr/overflow/, "nth_prime_approx(maxindex+$add) overflows");

eval { nth_prime($maxindex+$add); };
like($@, qr/overflow/, "nth_prime(maxindex+$add) overflows");
