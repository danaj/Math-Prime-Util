#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw(
  factor factor_exp is_prob_prime is_bpsw_prime is_semiprime
  is_square_free moebius liouville prime_omega prime_bigomega
  is_omega_prime is_almost_prime vecprod
);

$| = 1;

sub progress {
  my ($msg) = @_;
  diag($msg);
}

sub canonical_decimal {
  my ($n) = @_;
  $n = "$n";
  $n =~ s/\A\+//;
  $n =~ s/\A0+(?=\d)//;
  return $n;
}

sub factor_is {
  my ($n, $exp, $name) = @_;
  my @got = factor($n);
  is_deeply(\@got, $exp, $name);
  is("" . vecprod(@got), canonical_decimal($n), "$name product");
}

sub factor_exp_is {
  my ($n, $exp, $name) = @_;
  my @got = factor_exp($n);
  is_deeply(\@got, $exp, $name);
}

my $p70     = "1180591620717411303427";                  # 2^70 + 3
my $p70_fac = [ 13, 25087, "3619992029943217" ];

my $u64p1    = "18446744073709551617";                   # 2^64 + 1
my $u64p1_fe = [ [274177, 1], ["67280421310721", 1] ];

my $m127    = "170141183460469231731687303715884105727"; # 2^127 - 1
my $twom127 = "340282366920938463463374607431768211454"; # 2*(2^127-1)

my $p128    = "340282366920938463463374607431768211297"; # 128-bit prime
my $u128max = "340282366920938463463374607431768211455"; # 2^128 - 1
my $u128p1  = "340282366920938463463374607431768211456"; # 2^128

my $u128max_fac = [
  3, 5, 17, 257, 641, 65537, 274177,
  6700417, "67280421310721",
];

my @bad_bitint_examples = (
  "193188397108589712782271002720077813901",
  "242584098619172565144193751769573258437",
  "199259740502625294193397878845987884783",
);

progress("70-bit warmup and product checks");
factor_is($p70, $p70_fac, "factor 70-bit composite");
factor_exp_is($p70,
  [ [13,1], [25087,1], ["3619992029943217",1] ],
  "factor_exp 70-bit composite");

progress("64-bit boundary values");
factor_is($u64p1, [274177, "67280421310721"], "factor 2^64+1");
factor_exp_is($u64p1, $u64p1_fe, "factor_exp 2^64+1");
is(is_semiprime($u64p1), 1, "is_semiprime 2^64+1");
is(is_square_free($u64p1), 1, "is_square_free 2^64+1");
is(moebius($u64p1), 1, "moebius 2^64+1");

progress("128-bit parser and large factor output");
factor_is($p128, [$p128], "factor 128-bit prime");
factor_exp_is($p128, [ [$p128, 1] ], "factor_exp 128-bit prime");
factor_is("+000000000000000000000$p128", [$p128], "factor leading-zero 128-bit prime");
is(is_prob_prime($p128), 1, "is_prob_prime 128-bit prime");
is(is_bpsw_prime($p128), 1, "is_bpsw_prime 128-bit prime");
is(is_prob_prime("+000000000000000000000$p128"), 1, "is_prob_prime leading-zero 128-bit prime");
is(is_bpsw_prime("000000000000000000000$p128"), 1, "is_bpsw_prime leading-zero 128-bit prime");
is(is_semiprime($p128), 0, "is_semiprime 128-bit prime");
is(liouville($p128), -1, "liouville 128-bit prime");

progress("mixed small factor plus 127-bit large factor");
factor_is($twom127, [2, $m127], "factor 2*(2^127-1)");
factor_exp_is($twom127, [ [2,1], [$m127,1] ], "factor_exp 2*(2^127-1)");
is(is_semiprime($twom127), 1, "is_semiprime 2*(2^127-1)");
is(is_square_free($twom127), 1, "is_square_free 2*(2^127-1)");
is(is_square_free("-$twom127"), 1, "is_square_free negative 2*(2^127-1)");
is(moebius($twom127), 1, "moebius 2*(2^127-1)");
is(moebius("-$twom127"), 1, "moebius negative 2*(2^127-1)");
is(liouville($twom127), 1, "liouville 2*(2^127-1)");
is(prime_omega($twom127), 2, "prime_omega 2*(2^127-1)");
is(prime_bigomega($twom127), 2, "prime_bigomega 2*(2^127-1)");
is(is_omega_prime(2, $twom127), 1, "is_omega_prime 2*(2^127-1)");
is(is_almost_prime(2, $twom127), 1, "is_almost_prime 2*(2^127-1)");

progress("2^128-1 exact upper bound");
factor_is($u128max, $u128max_fac, "factor 2^128-1");
factor_exp_is($u128max, [ map { [$_, 1] } @$u128max_fac ], "factor_exp 2^128-1");
is(is_prob_prime($u128max), 0, "is_prob_prime 2^128-1");
is(is_bpsw_prime($u128max), 0, "is_bpsw_prime 2^128-1");
is(prime_omega($u128max), 9, "prime_omega 2^128-1");
is(prime_bigomega($u128max), 9, "prime_bigomega 2^128-1");
is(is_omega_prime(9, $u128max), 1, "is_omega_prime 2^128-1");
is(is_almost_prime(9, $u128max), 1, "is_almost_prime 2^128-1");

progress("just over 128-bit fallback on easy input");
is(scalar factor($u128p1), 128, "scalar factor 2^128 fallback count");
is(is_prob_prime($u128p1), 0, "is_prob_prime 2^128 fallback");
is(is_bpsw_prime($u128p1), 0, "is_bpsw_prime 2^128 fallback");

if (!defined $ENV{MPU_FACTOR128_SLOW} || $ENV{MPU_FACTOR128_SLOW}) {
  progress("_BitInt(128) regression examples");
  for my $n (@bad_bitint_examples) {
    is(is_semiprime($n), 1, "is_semiprime previous _BitInt regression $n");
  }
} else {
  progress("skipping slow _BitInt(128) regression examples; set MPU_FACTOR128_SLOW=1 to run");
  pass("slow _BitInt regression examples skipped");
}

done_testing();
