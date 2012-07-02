#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util;

use Test::More  tests => 1;

my @functions =  qw(
                     prime_get_config
                     prime_precalc prime_memfree
                     is_prime is_prob_prime miller_rabin
                     primes
                     next_prime  prev_prime
                     prime_count prime_count_lower prime_count_upper prime_count_approx
                     nth_prime nth_prime_lower nth_prime_upper nth_prime_approx
                     random_prime random_ndigit_prime
                     factor all_factors moebius euler_phi
                     ExponentialIntegral LogarithmicIntegral RiemannR
                   );
can_ok( 'Math::Prime::Util', @functions);
