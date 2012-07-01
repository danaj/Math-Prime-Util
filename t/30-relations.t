#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/primes
                         nth_prime nth_prime_lower nth_prime_upper nth_prime_approx
                         prime_count prime_count_lower prime_count_upper prime_count_approx
                         next_prime prev_prime
                        /;


my @trials = qw/1 2 3 4 5 6 7 17 57 89 102 1337 8573 84763 784357 1000001 2573622/;

my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
my $extra = defined $ENV{RELEASE_TESTING} && $ENV{RELEASE_TESTING};

plan tests => 5 * scalar @trials;

my $last = 0;
foreach my $n (@trials) {

  is( prime_count($n), scalar @{primes($n)}, "Prime count and scalar primes agree for $n" );

  is( prime_count($n) - prime_count($last),
      scalar @{primes( $last+1, $n )},
      "scalar primes($last+1,$n) = prime_count($n) - prime_count($last)" );

  is( prime_count(nth_prime($n)), $n, "Pi(pn)) = n for $n");

  is( nth_prime(prime_count($n)+1), next_prime($n), "p(Pi(n)+1) = next_prime(n) for $n" );

  is( nth_prime(prime_count($n)), prev_prime($n+1), "p(Pi(n)) = prev_prime(n) for $n" );

  $last = $n;
}
