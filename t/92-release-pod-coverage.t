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
    also_private => [
                      qr/^(erat|segment|trial|sieve|segment_twin)_primes$/,
                      'semi_prime_sieve',
                    ],
  };
  $param->{trustme} = [mpu_public_regex(), mpu_factor_regex()]
    if $m eq 'Math::Prime::Util::PP';
  $param->{trustme} = [mpu_public_regex(), mpu_factor_regex(), qw/rand srand/]
    if $m eq 'ntheory';
  pod_coverage_ok( $m, $param );
}

sub mpu_public_regex {
  my @funcs = qw(
      prime_get_config prime_set_config
      prime_precalc prime_memfree
      is_prime is_prob_prime is_provable_prime is_provable_prime_with_cert
      prime_certificate verify_prime
      is_pseudoprime is_euler_pseudoprime is_strong_pseudoprime
      is_euler_plumb_pseudoprime
      is_lucas_pseudoprime
      is_strong_lucas_pseudoprime
      is_extra_strong_lucas_pseudoprime
      is_almost_extra_strong_lucas_pseudoprime
      is_frobenius_pseudoprime
      is_frobenius_underwood_pseudoprime is_frobenius_khashin_pseudoprime
      is_perrin_pseudoprime is_catalan_pseudoprime
      is_aks_prime is_bpsw_prime is_ramanujan_prime is_mersenne_prime
      is_power is_prime_power is_pillai is_semiprime is_square is_polygonal
      is_square_free is_primitive_root is_carmichael is_quasi_carmichael
      is_fundamental is_totient
      sqrtint rootint logint
      miller_rabin_random
      lucas_sequence lucasu lucasv
      primes twin_primes semi_primes ramanujan_primes
      sieve_prime_cluster sieve_range
      forprimes forcomposites foroddcomposites forsemiprimes fordivisors
      forpart forcomp forcomb forperm forderange formultiperm forsetproduct
      forfactored forsquarefree
      lastfor
      numtoperm permtonum randperm shuffle
      prime_iterator prime_iterator_object
      next_prime  prev_prime
      prime_count
      prime_count_lower prime_count_upper prime_count_approx
      nth_prime nth_prime_lower nth_prime_upper nth_prime_approx inverse_li
      twin_prime_count twin_prime_count_approx
      nth_twin_prime nth_twin_prime_approx
      semiprime_count semiprime_count_approx
      nth_semiprime nth_semiprime_approx
      ramanujan_prime_count ramanujan_prime_count_approx
      ramanujan_prime_count_lower ramanujan_prime_count_upper
      nth_ramanujan_prime nth_ramanujan_prime_approx
      nth_ramanujan_prime_lower nth_ramanujan_prime_upper
      sum_primes print_primes
      random_prime random_ndigit_prime random_nbit_prime random_strong_prime
      random_proven_prime random_proven_prime_with_cert
      random_maurer_prime random_maurer_prime_with_cert
      random_shawe_taylor_prime random_shawe_taylor_prime_with_cert
      random_semiprime random_unrestricted_semiprime
      random_factored_integer
      primorial pn_primorial consecutive_integer_lcm gcdext chinese
      gcd lcm factor factor_exp divisors valuation hammingweight
      todigits fromdigits todigitstring sumdigits
      invmod sqrtmod addmod mulmod divmod powmod
      vecsum vecmin vecmax vecprod vecreduce vecextract
      vecany vecall vecnotall vecnone vecfirst vecfirstidx
      moebius mertens euler_phi jordan_totient exp_mangoldt liouville
      partitions bernfrac bernreal harmfrac harmreal
      chebyshev_theta chebyshev_psi
      divisor_sum carmichael_lambda kronecker hclassno
      ramanujan_tau ramanujan_sum
      binomial stirling znorder znprimroot znlog legendre_phi
      factorial factorialmod
      ExponentialIntegral LogarithmicIntegral RiemannZeta RiemannR LambertW Pi
      irand irand64 drand urandomb urandomm csrand random_bytes entropy_bytes
  );
  my $pattern = '^(' . join('|', @funcs) . ')$';
  return qr/$pattern/;
}

sub mpu_factor_regex {
  my @funcs = (qw/trial_factor fermat_factor holf_factor lehman_factor squfof_factor prho_factor pbrent_factor pminus1_factor pplus1_factor ecm_factor/);
  my $pattern = '^(' . join('|', @funcs) . ')$';
  return qr/$pattern/;
}
