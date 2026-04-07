#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/integer_complexity/;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};

plan tests => 4;

subtest 'small values', sub {
  # f(n) = min ones to represent n using + and *
  # f(6)=5: (1+1)*(1+1+1)     f(8)=6: (1+1)*(1+1)*(1+1)
  # f(9)=6: (1+1+1)*(1+1+1)   f(12)=7: (1+1+1)*(1+1)*(1+1)
  # Data from OEIS, 1..100
  my @A = (1,2,3,4,5,5,6,6,6,7,8,7,8,8,8,8,9,8,9,9,9,10,11,9,10,10,9,10,11,10,11,10,11,11,11,10,11,11,11,11,12,11,12,12,11,12,13,11,12,12,12,12,13,11,12,12,12,13,14,12,13,13,12,12,13,13,14,13,14,13,14,12,13,13,13,13,14,13,14,13,12,13,14,13,14,14,14,14,15,13,14,14,14,15,14,13,14,14,14,14);
  is_deeply([map { integer_complexity($_) } 1..100],\@A,"integer_complexity(1..100)");

  is(integer_complexity(0),undef,"integer_complexity(0) returns undef");
};

subtest 'powers of 2: f(2^k) = 2k', sub {
  # 2 = 1+1 uses 2 ones; each doubling adds f(2)=2
  is_deeply([map { integer_complexity(2**$_) } 1..16],
            [map { 2*$_ } 1..16],
            "f(2^k) = 2k for k=1..16");
};

subtest 'powers of 3: f(3^k) = 3k', sub {
  # 3 = 1+1+1 uses 3 ones; each tripling adds f(3)=3
  is_deeply([map { integer_complexity(3**$_) } 1..10],
            [map { 3*$_ } 1..10],
            "f(3^k) = 3k for k=1..10");
};

subtest 'selected values', sub {
  # f(2^a * 3^b) = 2a + 3b: optimal factorization tree uses only 2s and 3s
  is(integer_complexity(4 * 9),   2*2 + 3*2, "f(4*9) = f(2^2*3^2) = 10");
  is(integer_complexity(8 * 27),  2*3 + 3*3, "f(8*27) = f(2^3*3^3) = 15");
  is(integer_complexity(16 * 81), 2*4 + 3*4, "f(16*81) = f(2^4*3^4) = 20");

  if ($extra) {
    # Just above the 1M static cache boundary
    is(integer_complexity(2**20), 40, "f(2^20) = 40 (above cache)");
    is(integer_complexity(3**12), 36, "f(3^12) = 36 (above cache)");
  }
};
