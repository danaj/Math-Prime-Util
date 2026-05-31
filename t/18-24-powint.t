#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/powint negint/;

plan tests => 5;

subtest 'small values', sub {
  for my $a (-3 .. 3) {
    my @got = map { powint($a, $_) } 0..3;
    my @exp = (1, $a, $a*$a, $a*$a*$a);
    is_deeply( \@got, \@exp, "powint($a,0..3) = [@got] expect [@exp]" );
  }
};

subtest 'selected values', sub {
  my @vals = (
   [5, 6, 15625],
   [2, 16, 65536],
   [-544, 7, "-14099129446552305664"],
   [4294967295, 2, "18446744065119617025"],
   [2097152, 3, "9223372036854775808"],
   [2097153, 3, "9223385231000600577"],
   [65535, 4, "18445618199572250625"],
   [65536, 4, "18446744073709551616"],
   [4096, 5, "1152921504606846976"],
   [4097, 5, "1154329566852960257"],
   [16, 13, "4503599627370496"],
   [17, 13, "9904578032905937"],
   [8, 17, "2251799813685248"],
   [9, 17, "16677181699666569"],
  );

  is_deeply( [map { "".powint($_->[0],$_->[1]) } @vals],
             [map { $_->[2] } @vals],
             "powint a**b=c" );
};

subtest 'bigints', sub {
  is("".powint(powint(2,32),3),"79228162514264337593543950336","(2^32)^3");
  is("".powint(3,powint(2,7)),"11790184577738583171520872861412518665678211592275841109096961","3^(2^7)");
  ok(ref(powint(46,22)), "powint returns a bigint for 46,22");
  ok(ref(powint(-544,7)), "powint returns a bigint for -544,7");
};

subtest 'edge cases: 0^0, 0^n, 1^n, (-1)^n', sub {
  is(powint(0,0),  1, "0^0 = 1");
  is(powint(0,1),  0, "0^1 = 0");
  is(powint(0,5),  0, "0^5 = 0");
  is(powint(1,0),  1, "1^0 = 1");
  is(powint(1,100),1, "1^100 = 1");
  is_deeply([map { powint(-1,$_) } 0..7],
            [1,-1,1,-1,1,-1,1,-1],
            "(-1)^n alternates 1,-1");
  my $huge_even = "1000000000000000000000000000000";
  my $huge_odd  = "1000000000000000000000000000001";
  is(powint(0,  $huge_even),  0, "0^huge = 0");
  is(powint(1,  $huge_even),  1, "1^huge = 1");
  is(powint(-1, $huge_even),  1, "(-1)^huge_even = 1");
  is(powint(-1, $huge_odd),  -1, "(-1)^huge_odd = -1");
};

subtest 'negative bases parity check', sub {
  is( powint(-7,6),  117649, "powint(-7,6) ==  powint(7,6) ==  117649" );
  is( powint(-7,7), -823543, "powint(-7,7) == -powint(7,7) == -823543" );
};
