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
                      qr/^(erat|segment|trial|sieve)_primes$/,
                    ],
  };
  $param->{trustme} = [mpu_public_regex(), mpu_factor_regex(), mpu_PPM_regex()]
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
      is_delicate_prime is_chen_prime
      is_odd is_even is_divisible is_congruent
      is_power is_prime_power is_perfect_power is_square
      is_square_free is_powerfree
      is_pillai is_polygonal is_congruent_number is_perfect_number
      is_semiprime is_almost_prime is_omega_prime
      is_primitive_root is_carmichael is_quasi_carmichael
      is_fundamental is_totient is_gaussian_prime is_sum_of_squares
      is_smooth is_rough is_powerful is_practical is_lucky
      sqrtint rootint logint lshiftint rshiftint rashiftint absint negint
      signint cmpint addint subint add1int sub1int mulint powint
      divint modint cdivint divrem fdivrem cdivrem tdivrem
      miller_rabin_random
      lucas_sequence
      lucasu lucasv lucasuv lucasumod lucasvmod lucasuvmod pisano_period
      primes twin_primes semi_primes almost_primes omega_primes ramanujan_primes
      sieve_prime_cluster sieve_range prime_powers lucky_numbers
      forprimes forcomposites foroddcomposites forsemiprimes foralmostprimes
      forpart forcomp forcomb forperm forderange formultiperm forsetproduct
      fordivisors forfactored forsquarefree
      lastfor
      numtoperm permtonum randperm shuffle
      prime_iterator prime_iterator_object
      next_prime prev_prime
      next_prime_power prev_prime_power
      next_perfect_power prev_perfect_power
      next_chen_prime
      prime_count prime_count_lower prime_count_upper prime_count_approx
      nth_prime nth_prime_lower nth_prime_upper nth_prime_approx inverse_li
      twin_prime_count twin_prime_count_approx
      nth_twin_prime nth_twin_prime_approx
      semiprime_count semiprime_count_approx
      nth_semiprime nth_semiprime_approx
      almost_prime_count almost_prime_count_approx
      almost_prime_count_lower almost_prime_count_upper
      nth_almost_prime nth_almost_prime_approx
      nth_almost_prime_lower nth_almost_prime_upper
      omega_prime_count nth_omega_prime
      ramanujan_prime_count ramanujan_prime_count_approx
      ramanujan_prime_count_lower ramanujan_prime_count_upper
      nth_ramanujan_prime nth_ramanujan_prime_approx
      nth_ramanujan_prime_lower nth_ramanujan_prime_upper
      powerful_count nth_powerful sumpowerful powerful_numbers
      prime_power_count prime_power_count_approx
      prime_power_count_lower prime_power_count_upper
      nth_prime_power nth_prime_power_approx
      nth_prime_power_lower nth_prime_power_upper
      perfect_power_count perfect_power_count_approx
      perfect_power_count_lower perfect_power_count_upper
      nth_perfect_power nth_perfect_power_approx
      nth_perfect_power_lower nth_perfect_power_upper
      nth_powerfree powerfree_count powerfree_sum squarefree_kernel
      powerfree_part powerfree_part_sum
      smooth_count rough_count powersum
      lucky_count lucky_count_approx lucky_count_lower lucky_count_upper
      nth_lucky nth_lucky_approx nth_lucky_lower nth_lucky_upper
      sum_primes print_primes
      random_prime random_ndigit_prime
      random_nbit_prime random_safe_prime random_strong_prime
      random_proven_prime random_proven_prime_with_cert
      random_maurer_prime random_maurer_prime_with_cert
      random_shawe_taylor_prime random_shawe_taylor_prime_with_cert
      random_semiprime random_unrestricted_semiprime
      random_factored_integer
      primorial pn_primorial consecutive_integer_lcm gcdext chinese chinese2
      gcd lcm factor factor_exp divisors valuation hammingweight
      frobenius_number
      todigits fromdigits todigitstring sumdigits
      tozeckendorf fromzeckendorf
      sqrtmod allsqrtmod rootmod allrootmod cornacchia
      negmod invmod addmod submod mulmod divmod powmod muladdmod mulsubmod
      vecsum vecmin vecmax vecprod vecreduce vecextract vecequal vecuniq
      vecany vecall vecnotall vecnone vecfirst vecfirstidx vecmex vecpmex
      setbinop sumset
      moebius mertens liouville sumliouville prime_omega prime_bigomega
      euler_phi jordan_totient exp_mangoldt sumtotient
      partitions bernfrac bernreal harmfrac harmreal
      chebyshev_theta chebyshev_psi
      divisor_sum carmichael_lambda hclassno inverse_totient
      kronecker is_qr qnr
      ramanujan_tau ramanujan_sum
      stirling fubini znorder znprimroot znlog legendre_phi
      factorial factorialmod subfactorial binomial binomialmod
      falling_factorial rising_factorial
      ExponentialIntegral LogarithmicIntegral RiemannZeta RiemannR LambertW Pi
      irand irand64 drand urandomb urandomm csrand random_bytes entropy_bytes
  );
  my $pattern = '^(' . join('|', @funcs) . ')$';
  return qr/$pattern/;
}

sub mpu_factor_regex {
  my @funcs = (qw/trial_factor fermat_factor holf_factor lehman_factor squfof_factor prho_factor pbrent_factor pminus1_factor pplus1_factor cheb_factor ecm_factor/);
  my $pattern = '^(' . join('|', @funcs) . ')$';
  return qr/$pattern/;
}

sub mpu_PPM_regex {
  my @funcs = qw(
      Maddint Msubint Mmulint Mdivint Mcdivint Mpowint Mabsint Mnegint
      Mmodint Mlogint Mrootint Msqrtint Mcmpint
      Mlshiftint Mrshiftint
      Maddmod Msubmod Mmulmod Mrootmod Mdivmod Mpowmod Minvmod
      Mmuladdmod Mmulsubmod
      Mbinomial Mchinese Mfactor Mfactor_exp Mdivisors
      Mfactorial Mfordivisors Mforprimes Mgcd Mlcm Mznorder
      Mfalling_factorial
      Mkronecker Mmoebius Mtotient Mstirling
      Mnext_prime Mprev_prime Mprime_count
      Mprimorial Mpn_primorial
      MLi Mpowersum
      Murandomb Murandomm Mvaluation
      Mprimes
      Mis_power Mis_prime Mis_prime_power Mis_odd Mis_even
      Mis_semiprime Mis_square_free
      Mvecall Mvecany Mvecmin Mvecmax Mvecnone Mvecprod Mvecsum
      Mtodigits
      Mlucasumod Mvecfirst
      Saddint Ssubint Smulint Sdivint Spowint
      reftyped
      validate_integer     validate_integer_nonneg
      validate_integer_abs validate_integer_positive
  );
  my $pattern = '^(' . join('|', @funcs) . ')$';
  return qr/$pattern/;
}
