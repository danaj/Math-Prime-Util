#!/usr/bin/env perl
use strict;
use warnings;

BEGIN {
  unless ($ENV{RELEASE_TESTING}) {
    require Test::More;
    Test::More::plan(skip_all => 'these tests are for release candidate testing');
  }
}

#---------------------------------------------------------------------


use Test::More;
eval "use Test::Pod::Coverage 1.08";
plan skip_all => "Test::Pod::Coverage 1.08 required for testing POD coverage"
  if $@;

my @modules = grep { $_ ne 'Math::Prime::Util::PPFE' }
              Test::Pod::Coverage::all_modules();

plan tests => scalar @modules;

#my $ppsubclass = { trustme => [mpu_public_regex()] };

foreach my $m (@modules) {
  my $param = {
    also_private => [ qr/^(erat|segment|trial|sieve|segment_twin)_primes$/ ],
  };
  $param->{trustme} = [mpu_public_regex(), mpu_factor_regex()]
    if $m eq 'Math::Prime::Util::PP';
  $param->{trustme} = [mpu_public_regex(), mpu_factor_regex()]
    if $m eq 'ntheory';
  pod_coverage_ok( $m, $param );
}

sub mpu_public_regex {
  my @funcs =
  qw/ prime_get_config prime_set_config
      prime_precalc prime_memfree
      is_prime is_prob_prime is_provable_prime is_provable_prime_with_cert
      prime_certificate verify_prime
      is_pseudoprime is_strong_pseudoprime
      is_lucas_pseudoprime
      is_strong_lucas_pseudoprime
      is_extra_strong_lucas_pseudoprime
      is_almost_extra_strong_lucas_pseudoprime
      is_frobenius_underwood_pseudoprime
      is_frobenius_pseudoprime
      is_perrin_pseudoprime
      is_aks_prime is_bpsw_prime
      is_power
      miller_rabin_random
      lucas_sequence
      primes twin_primes
      forprimes forcomposites foroddcomposites fordivisors
      forpart forcomb forperm
      prime_iterator prime_iterator_object
      next_prime  prev_prime
      prime_count
      prime_count_lower prime_count_upper prime_count_approx
      nth_prime nth_prime_lower nth_prime_upper nth_prime_approx
      twin_prime_count twin_prime_count_approx
      nth_twin_prime nth_twin_prime_approx
      random_prime random_ndigit_prime random_nbit_prime random_strong_prime
      random_proven_prime random_proven_prime_with_cert
      random_maurer_prime random_maurer_prime_with_cert
      random_shawe_taylor_prime random_shawe_taylor_prime_with_cert
      primorial pn_primorial consecutive_integer_lcm gcdext chinese
      gcd lcm factor factor_exp divisors valuation invmod
      vecsum vecprod vecmin vecmax
      moebius mertens euler_phi jordan_totient exp_mangoldt liouville
      partitions bernfrac bernreal
      chebyshev_theta chebyshev_psi
      divisor_sum carmichael_lambda
      kronecker binomial factorial znorder znprimroot znlog legendre_phi
      ExponentialIntegral LogarithmicIntegral RiemannZeta RiemannR LambertW Pi
  /;
  my $pattern = '^(' . join('|', @funcs) . ')$';
  return qr/$pattern/;
}

sub mpu_factor_regex {
  my @funcs = (qw/trial_factor fermat_factor holf_factor squfof_factor prho_factor pbrent_factor pminus1_factor pplus1_factor ecm_factor/);
  my $pattern = '^(' . join('|', @funcs) . ')$';
  return qr/$pattern/;
}
