#!/usr/bin/env perl
use strict;
use warnings;

BEGIN { $ENV{MPU_NO_GMP} = 1 }

use blib;
use Time::HiRes qw(time);
use Math::Prime::Util qw(factor modint mulint powint prev_prime subint);

# tinyqs_ms measures one nontrivial split; factor_ms measures complete public
# factorization.  This is a crossover signal, not equal-work microbenchmarking.
die "tinyqs128 requires XS with uint128_t support\n"
  unless Math::Prime::Util->can('_XS_tinyqs128') &&
         Math::Prime::Util::_XS_has_uint128();

my @bits = @ARGV ? @ARGV : (65, 72, 80, 88, 96, 104, 112,
                            120, 124, 126, 127, 128);
my $cases_per_size = $ENV{TINYQS_CASES} || 4;

sub make_cases {
  my ($bits) = @_;
  my @cases;
  my $pbits = int($bits / 2);
  my $qbits = $bits - $pbits;

  for my $i (1 .. $cases_per_size) {
    my $p = prev_prime(subint(powint(2, $pbits),
                              1009 + 127*$bits + 7919*$i));
    my $q = prev_prime(subint(powint(2, $qbits),
                              100003 + 313*$bits + 104729*$i));
    ($p, $q) = ($q, $p) if $p > $q;
    push @cases, [mulint($p, $q), $p, $q];
  }
  return \@cases;
}

sub elapsed {
  my ($code, $cases) = @_;
  my $checksum = 0;
  my $start = time;
  $checksum += $code->($_->[0]) for @$cases;
  return (time-$start, $checksum);
}

printf "%4s  %9s  %9s  %8s\n", qw(bits tinyqs_ms factor_ms speedup);
for my $bits (@bits) {
  my $cases = make_cases($bits);

  for my $case (@$cases) {
    my ($n, $p, $q) = @$case;
    my $f = Math::Prime::Util::_XS_tinyqs128($n);
    die "tinyqs128 failed for $n\n"
      unless ($f == $p || $f == $q) && modint($n, $f) == 0;
    my @f = factor($n);
    die "factor failed for $n\n"
      unless @f == 2 && $f[0] == $p && $f[1] == $q;
  }

  my ($tqs) = elapsed(
    sub { length("" . Math::Prime::Util::_XS_tinyqs128($_[0])) }, $cases);
  my ($tcur) = elapsed(
    sub { my @f = factor($_[0]); scalar(@f) + length("" . $f[0]) }, $cases);

  printf "%4d  %9.3f  %9.3f  %7.2fx\n",
         $bits, 1000*$tqs/$cases_per_size,
         1000*$tcur/$cases_per_size, $tcur/$tqs;
}
