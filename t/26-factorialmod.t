#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/factorialmod factorial/;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;

plan tests => 9;

{
  my @result = map { my $m=$_; map { factorialmod($_,$m) } 0..$m-1; } 1 .. 40;
  my @expect = map { my $m=$_; map { "".(factorial($_) % $m); } 0..$m-1; } 1 .. 40;
  is_deeply( \@result, \@expect, "factorialmod n! mod m for m 1 to 50, n 0 to m" );
}

SKIP: {
  skip "medium size factorialmods in PP",2 unless $usexs || $extra;
  is( factorialmod(1000000000,1000000008), 0, "1000000000! mod 1000000008 is zero" );
  is( factorialmod(50000,10000019), 8482159, "50000! mod 10000019" );
}

SKIP: {
  skip "large value without EXTENDED_TESTING on 64-bit",1 unless $extra && $use64;
  is( "".factorialmod(5000001,"8000036000054000027"), "4179720539133404343", "factorialmod with large n and large composite non-square-free m" );
}

{
  my $m = 1000003;
  my @n = ($m-60 .. $m-1);
  my @got = map { factorialmod($_, $m) } @n;
  my @exp;
  $exp[0] = $got[0];
  for my $i (1 .. $#n) {
    $exp[$i] = ($exp[$i-1] * $n[$i]) % $m;
  }
  is_deeply(\@got, \@exp, "factorialmod recurrence near prime modulus (backward path)");
}

{
  my $m = 10000020;
  my @n = (1000 .. 1120);
  my @got = map { factorialmod($_, $m) } @n;
  my @exp;
  $exp[0] = $got[0];
  for my $i (1 .. $#n) {
    $exp[$i] = ($exp[$i-1] * $n[$i]) % $m;
  }
  is_deeply(\@got, \@exp, "factorialmod recurrence in segmented non-mont path");
}

{
  my $m = 10000019;
  my @n = (1000 .. 1120);
  my @got = map { factorialmod($_, $m) } @n;
  my @exp;
  $exp[0] = $got[0];
  for my $i (1 .. $#n) {
    $exp[$i] = ($exp[$i-1] * $n[$i]) % $m;
  }
  is_deeply(\@got, \@exp, "factorialmod recurrence in segmented mont path");
}

{
  my $m = 10403;  # 101 * 103
  is(factorialmod(102, $m), 5252,
     "factorialmod below composite-zero cutoff remains non-trivial");
  is(factorialmod(103, $m), 0,
     "factorialmod at composite-zero cutoff returns zero");
}
