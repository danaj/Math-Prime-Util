#!/usr/bin/env perl
use strict;
use warnings;
$| = 1;  # fast pipes

use Math::Prime::Util qw/chinese lcm/;
use Math::ModInt qw(mod);
use Math::ModInt::ChineseRemainder qw(cr_combine);

my $limit = shift || 10_000;
my $printmod = int(($limit+77) / 78);

print "Running $limit random tests with 2 small inputs...\n";
for my $n (1 .. $limit) {
  print '.' unless $n % $printmod;

  my (@rm) = map { my $m = 1+int(rand(2**16)); [int(rand($m)),$m] } 1..2;
  #print join(" ", map { "[$_->[0] $_->[1]]" } @rm), "\n";

  my $mic = cr_combine( map { mod($_->[0],$_->[1]) } @rm );
  if ($mic->is_undefined) {
    my $mpu_res = chinese(@rm);
    die "MIC: undef  MPU: $mpu_res\n" if defined $mpu_res;
    next;
  }
  my $mic_res = $mic->residue;
  my $mic_mod = $mic->modulus;

  my $mpu_res = chinese(@rm);
  die "MIC: $mic_res $mic_mod  MPU: undef\n" if !defined $mpu_res;
  my $mpu_mod = lcm(map { $_->[1] } @rm);
  die "MIC: $mic_res $mic_mod   MPU: $mpu_res  $mpu_mod" if $mpu_res != $mic_res || $mpu_mod != $mic_mod;
}

print "\nRunning $limit random tests with 2 large inputs...\n";
for my $n (1 .. $limit) {
  print '.' unless $n % $printmod;

  my (@rm) = map { my $m = 1+int(rand(2**32)); [int(rand($m)),$m] } 1..2;
  #print join(" ", map { "[$_->[0] $_->[1]]" } @rm), "\n";

  my $mic = cr_combine( map { mod($_->[0],$_->[1]) } @rm );
  if ($mic->is_undefined) {
    my $mpu_res = chinese(@rm);
    die "MIC: undef  MPU: $mpu_res\n" if defined $mpu_res;
    next;
  }
  my $mic_res = $mic->residue;
  my $mic_mod = $mic->modulus;

  my $mpu_res = chinese(@rm);
  die "MIC: $mic_res $mic_mod  MPU: undef\n" if !defined $mpu_res;
  my $mpu_mod = lcm(map { $_->[1] } @rm);
  die "MIC: $mic_res $mic_mod   MPU: $mpu_res  $mpu_mod" if $mpu_res != $mic_res || $mpu_mod != $mic_mod;
}

print "\nRunning $limit random tests with 4 small inputs...\n";
for my $n (1 .. $limit) {
  print '.' unless $n % $printmod;

  my (@rm) = map { my $m = 1+int(rand(2**16)); [int(rand($m)),$m] } 1..4;
  #print join(" ", map { "[$_->[0] $_->[1]]" } @rm), "\n";

  my $mic = cr_combine( map { mod($_->[0],$_->[1]) } @rm );
  if ($mic->is_undefined) {
    my $mpu_res = chinese(@rm);
    die "MIC: undef  MPU: $mpu_res\n" if defined $mpu_res;
    next;
  }
  my $mic_res = $mic->residue;
  my $mic_mod = $mic->modulus;

  my $mpu_res = chinese(@rm);
  die "MIC: $mic_res $mic_mod  MPU: undef\n" if !defined $mpu_res;
  my $mpu_mod = lcm(map { $_->[1] } @rm);
  die "MIC: $mic_res $mic_mod   MPU: $mpu_res  $mpu_mod" if $mpu_res != $mic_res || $mpu_mod != $mic_mod;
}
print "\nRunning $limit random tests with 3 large inputs...\n";
for my $n (1 .. $limit) {
  print '.' unless $n % $printmod;

  my (@rm) = map { my $m = 1+int(rand(2**40)); [int(rand($m)),$m] } 1..3;
  #print join(" ", map { "[$_->[0] $_->[1]]" } @rm), "\n";

  my $mic = cr_combine( map { mod($_->[0],$_->[1]) } @rm );
  if ($mic->is_undefined) {
    my $mpu_res = chinese(@rm);
    die "MIC: undef  MPU: $mpu_res\n" if defined $mpu_res;
    next;
  }
  my $mic_res = $mic->residue;
  my $mic_mod = $mic->modulus;

  my $mpu_res = chinese(@rm);
  die "MIC: $mic_res $mic_mod  MPU: undef\n" if !defined $mpu_res;
  my $mpu_mod = lcm(map { $_->[1] } @rm);
  die "MIC: $mic_res $mic_mod   MPU: $mpu_res  $mpu_mod" if $mpu_res != $mic_res || $mpu_mod != $mic_mod;
}
print "\nRunning $limit random tests with 13 large inputs...\n";
for my $n (1 .. $limit) {
  print '.' unless $n % $printmod;

  my (@rm) = map { my $m = 1+int(rand(2**31)); [int(rand($m)),$m] } 1..13;
  #print join(" ", map { "[$_->[0] $_->[1]]" } @rm), "\n";

  my $mic = cr_combine( map { mod($_->[0],$_->[1]) } @rm );
  if ($mic->is_undefined) {
    my $mpu_res = chinese(@rm);
    die "MIC: undef  MPU: $mpu_res\n" if defined $mpu_res;
    next;
  }
  my $mic_res = $mic->residue;
  my $mic_mod = $mic->modulus;

  my $mpu_res = chinese(@rm);
  die "MIC: $mic_res $mic_mod  MPU: undef\n" if !defined $mpu_res;
  my $mpu_mod = lcm(map { $_->[1] } @rm);
  die "MIC: $mic_res $mic_mod   MPU: $mpu_res  $mpu_mod" if $mpu_res != $mic_res || $mpu_mod != $mic_mod;
}
print "\nDone\n";
