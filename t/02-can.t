#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util;

use Test::More  tests => 1;

my @functions =  qw(
      prime_get_config prime_set_config
      prime_precalc prime_memfree
      is_prime is_prob_prime is_provable_prime is_provable_prime_with_cert
      prime_certificate verify_prime
      is_pseudoprime is_strong_pseudoprime
      is_lucas_pseudoprime
      is_strong_lucas_pseudoprime
      is_extra_strong_lucas_pseudoprime
      is_almost_extra_strong_lucas_pseudoprime
      is_frobenius_underwood_pseudoprime
      is_aks_prime
      miller_rabin_random
      lucas_sequence
      primes
      forprimes forcomposites
      prime_iterator prime_iterator_object
      next_prime  prev_prime
      prime_count
      prime_count_lower prime_count_upper prime_count_approx
      nth_prime nth_prime_lower nth_prime_upper nth_prime_approx
      random_prime random_ndigit_prime random_nbit_prime random_strong_prime
      random_proven_prime random_proven_prime_with_cert
      random_maurer_prime random_maurer_prime_with_cert
      primorial pn_primorial consecutive_integer_lcm
      factor factor_exp divisors
      moebius mertens euler_phi jordan_totient exp_mangoldt liouville
      partitions
      chebyshev_theta chebyshev_psi
      divisor_sum
      carmichael_lambda znorder
      ExponentialIntegral LogarithmicIntegral RiemannZeta RiemannR
                   );
can_ok( 'Math::Prime::Util', @functions);
