#!/usr/bin/env perl
use warnings;
use strict;
use Math::Prime::Util qw/:all/;
use List::Util qw/sum/;
use Math::BigInt try=>"GMP";

# This shows examples of many sequences from:
#   https://metacpan.org/release/Math-NumSeq
# Some of them are faster, some are much faster, a few are slower.
# This usually shows up once past ~ 10k values, or for large preds.
# These also do not have the limit of 2^32 of most Math::NumSeq implementations.
#
# Note that this completely lacks the framework of the module, and Math::NumSeq
# often implements various options that aren't always here.  It's just
# showing some examples of using MPU to solve these sort of problems.

# The lucas_sequence function covers about 45 different OEIS sequences,
# including Fibonacci, Lucas, Pell, Jacobsthal, Jacobsthal-Lucas, etc.

my $type = shift || 'AllPrimeFactors';
my $count = shift || 100;
my $arg = shift;  $arg = '' unless defined $arg;
my @n;

if      ($type eq 'Abundant') {
# TODO

} elsif ($type eq 'All') {
  print join " ", 1..$count;
} elsif ($type eq 'AllPrimeFactors') {
  my $i = 2;
  if ($arg eq 'descending') {
    push(@n, reverse factor($i++)) while scalar @n < $count;
  } else {
    push(@n, factor($i++)) while scalar @n < $count;
  }
  print join " ", @n[0..$count-1];
} elsif ($type eq 'AlmostPrimes') {
  $arg = 2 unless $arg =~ /^\d+$/;
  my $i = 1;
  while (@n < $count) {
    # use factor_exp for distinct
    $i++ while scalar factor($i) != $arg;
    push @n, $i++;
  }
  print join " ", @n;
} elsif ($type eq 'Cubes') {
  # Done via pred to show use
  my $i = 0;
  while (@n < $count) {
    $i++ while !is_power($i,3);
    push @n, $i++;
  }
  print join " ", @n;
} elsif ($type eq 'DedekindPsiCumulative') {
  my $c = 0;
  print join " ", map { $c += psi($_) } 1..$count;
} elsif ($type eq 'DivisorCount') {
  print join " ", map { scalar divisors($_) } 1..$count;
} elsif ($type eq 'DuffinianNumbers') {
  my $i = 0;
  while (@n < $count) {
    $i++ while !is_duffinian($i);
    push @n, $i++;
  }
  print join " ", @n;
} elsif ($type eq 'PolignacObstinate') {
  my $i = 1;
  while (@n < $count) {
    $i += 2 while !is_polignac_obstinate($i);
    push @n, $i;
    $i += 2;
  }
  print join " ", @n;
} elsif ($type eq 'Emirps') {
  my $i = 13;
  while (@n < $count) {
    $i = next_prime($i) while !is_prime(reverse $i) || $i eq reverse($i);
    push @n, $i;
    $i = next_prime($i);
  }
  print join " ", @n;
} elsif ($type eq 'Fibonacci') {
  # This is not a good way to do it, but does show a use for the function.
  my $lim = ~0;
  $lim = Math::BigInt->new(2) ** $count  if $count > 70;
  print join " ", map { (lucas_sequence($lim, 1, -1, $_))[0] } 0..$count-1;
} elsif ($type eq 'LemoineCount') {
  print join " ", map { lemoine_count($_) } 1..$count;
} elsif ($type eq 'LiouvilleFunction') {
  print join " ", map { liouville($_) } 1..$count;
} elsif ($type eq 'LucasNumbers') {
  # This is not a good way to do it, but does show a use for the function.
  my $lim = ~0;
  $lim = Math::BigInt->new(2) ** $count  if $count > 91;
  print join " ", map { (lucas_sequence($lim, 1, -1, $_))[1] } 1..$count;
} elsif ($type eq 'MobiusFunction') {
  print join " ", moebius(1,$count);
} elsif ($type eq 'MoranNumbers') {
  my $i = 1;
  while (@n < $count) {
    $i++ while !is_moran($i);
    push @n, $i++;
  }
  print join " ", @n;
} elsif ($type eq 'Pell') {
  # This is not a good way to do it, but does show a use for the function.
  my $lim = ~0;
  $lim = Math::BigInt->new(3) ** $count  if $count > 51;
  print join " ", map { (lucas_sequence($lim, 2, -1, $_))[0] } 0..$count-1;
} elsif ($type eq 'PowerFlip') {
  print join " ", map { powerflip($_) } 1..$count;
} elsif ($type eq 'Primes') {
  print join " ", @{primes($count)};
} elsif ($type eq 'PrimeFactorCount') {
  if ($arg eq 'distinct') {
    print join " ", map { scalar factor_exp($_) } 1..$count;
  } else {
    print join " ", map { scalar factor($_) } 1..$count;
  }
} elsif ($type eq 'PrimeIndexPrimes') {
  $arg = 2 unless $arg =~ /^\d+$/;
  print join " ", map { primeindexprime($_,$arg) } 1..$count;
} elsif ($type eq 'Primorials') {
  print join " ", map { pn_primorial($_) } 0..$count-1;
} elsif ($type eq 'SophieGermainPrimes') {
  my $estimate = sg_upper_bound($count);
  my $numfound = 0;
  forprimes {  push @n, $_ if is_prime(2*$_+1);  } $estimate;
  print join " ", @n[0..$count-1];
} elsif ($type eq 'Squares') {
  # Done via pred to show use
  my $i = 0;
  while (@n < $count) {
    $i++ while !is_power($i,2);
    push @n, $i++;
  }
  print join " ", @n;
} elsif ($type eq 'Totient') {
  print join " ", euler_phi(1,$count);
} elsif ($type eq 'TotientCumulative') {
  # pred:   sum(euler_phi(0,$_[0]));
  my $c = 0;
  print join " ", map { $c += euler_phi($_) } 0..$count-1;
} elsif ($type eq 'TotientPerfect') {
  my $i = 1;
  while (@n < $count) {
    $i += 2 while $i != totient_steps_sum($i,0);
    push @n, $i;
    $i += 2;
  }
  print join " ", @n;
} elsif ($type eq 'TotientSteps') {
  print join " ", map { totient_steps($_) } 1..$count;
} elsif ($type eq 'TotientStepsSum') {
  print join " ", map { totient_steps_sum($_) } 1..$count;
} elsif ($type eq 'TwinPrimes') {
  my $l = 2;
  my $upper = 400 + int(1.01 * nth_twin_prime_approx($count));
  $l=2; forprimes { push @n, $l if $l+2==$_; $l=$_; } $upper;
  print join " ", @n[0..$count-1];
} else {

# The following sequences, other than those marked TODO, do not exercise the
# features of MPU, hence there is little point reproducing them here.

# Abundant             TODO
# AlgebraicContinued
# AllDigits
# AsciiSelf
# BalancedBinary
# Base::IteraeIth
# Base::IteratePred
# BaumSweet
# Beastly
# Catalan
# CollatzSteps
# ConcatNumbers
# CullenNumbers
# DedekindPsiSteps     TODO
# DeleteablePrimes     TODO
# DigitCount
# DigitCountHigh
# DigitCountLow
# DigitLength
# DigitLengthCumulative
# DigitProduct
# DigitProductSteps
# DigitSum
# DigitSumModulo
# ErdosSelfridgeClass  TODO
# Even
# Expression
# Factorials
# Fibbinary
# FibbinaryBitCount
# FibonacciRepresentations
# FibonacciWord
# File
# FractionDigits
# GolayRudinShapiro
# GolayRudinShapiroCumulative
# GoldbachCount        TODO
# GolombSequence
# HafermanCarpet
# HappyNumbers
# HappySteps
# HarshadNumbers
# HofstadterFigure
# JugglerSteps
# KlarnerRado
# Kolakoski
# LuckyNumbers
# MaxDigitCount
# MephistoWaltz
# Modulo
# Multiples
# NumAronson
# OEIS
# OEIS::Catalogue
# OEIS::Catalogue::Plugin
# Odd
# Palindromes
# Perrin
# PisanoPeriod
# PisanoPeriodSteps
# Polygonal
# PowerPart          TODO
# Powerful           TODO
# PrimeIndexOrder    TODO
# Pronic
# ProthNumbers
# PythagoranHypots
# RadixConversion
# RadixWithoutDigit
# ReReplace
# ReRound
# RepdigitAny
# RepdigitRadix
# Repdigits
# ReverseAdd
# ReverseAddSteps
# Runs
# SelfLengthCumulative
# SpiroFibonacci
# SqrtContinued
# SqrtContinuedPeriod
# SqrtDigits
# SqrtEngel
# StarNumbers
# SternDiatomic
# Tetrahedral
# Triangular
# UlamSequence
# UndulatingNumbers
# WoodallNumbers

  die "sequence '$type' is not implemented here\n";
}
print "\n";
exit(0);


# DedekindPsi
sub psi { jordan_totient(2,$_[0])/jordan_totient(1,$_[0]) }

sub is_duffinian {
  my $n = shift;
  return 0 if $n < 4 || is_prime($n);
  my $dsum = divisor_sum($n);
  foreach my $d (divisors($n)) {
    return 0 unless $d == 1 || $dsum % $d;
  }
  1;
}

sub is_moran {
  my $n = shift;
  my $digsum = sum(split('',$n));
  return 0 if $n % $digsum;
  return 0 unless is_prime($n/$digsum);
  1;
}

sub is_polignac_obstinate {
  my $n = shift;
  return (0,1,0,0)[$n] if $n <= 3;
  return 0 unless $n & 1;
  my $k = 1;
  while (($n >> $k) > 0) {
    return 0 if is_prime($n - (1 << $k));
    $k++;
  }
  1;
}
  

# Lemoine Count (A046926)
sub lemoine_count {
  my($n, $count) = (shift, 0);
  return is_prime(($n>>1)-1) ? 1 : 0 unless $n & 1;
  forprimes { $count++ if is_prime($n-2*$_) } $n>>1;
  $count;
}

sub powerflip {
  my($n, $prod) = (shift, 1);
  # The spiffy log solution for bigints taken from Math::NumSeq
  my $log = 0;
  foreach my $pe (factor_exp($n)) {
    my ($p,$e) = @$pe;
    $log += $p * log($e);
    $e = Math::BigInt->new($e) if $log > 31;
    $prod *= $e ** $p;
  }
  $prod;
}

sub primeindexprime {
  my($n,$level) = @_;
  $n = nth_prime($n) for 1..$level;
  $n;
}

# TotientSteps
sub totient_steps {
  my($n, $count) = (shift,0);
  while ($n > 1) {
    $n = euler_phi($n);
    $count++;
  }
  $count;
}

# TotientStepsSum
sub totient_steps_sum {
  my $n = shift;
  my $sum = shift;  $sum = $n unless defined $sum;
  while ($n > 1) {
    $n = euler_phi($n);
    $sum += $n;
  }
  $sum;
}

# Sophie-Germaine primes upper bound.  Messy.
sub sg_upper_bound {
  my $count = shift;
  my $nth = nth_prime_upper($count);
  # For lack of a better formula, do this step-wise estimate.
  my $estimate = ($count <   5000) ? 150 + int( $nth * log($nth) * 1.2 )
               : ($count <  19000) ? int( $nth * log($nth) * 1.135 )
               : ($count <  45000) ? int( $nth * log($nth) * 1.10 )
               : ($count < 100000) ? int( $nth * log($nth) * 1.08 )
               : ($count < 165000) ? int( $nth * log($nth) * 1.06 )
               : ($count < 360000) ? int( $nth * log($nth) * 1.05 )
               : ($count < 750000) ? int( $nth * log($nth) * 1.04 )
               : ($count <1700000) ? int( $nth * log($nth) * 1.03 )
               :                     int( $nth * log($nth) * 1.02 );

  return $estimate;
}
