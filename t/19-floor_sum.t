#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/floor_sum/;

my $use64 = ~0 > 4294967295 && ~0 != 18446744073709550592;

my @tests = (
  [ 0, 7, 3, 4, 0 ],
  [ 1, 7, 3, 4, 0 ],
  [ 3, 3, 2, 1, 2 ],
  [ 3, 3, 4, 1, 4 ], # testing a >= m path
  [ 5, 7, 3, 4, 5 ],
  [ 4, 10, 6, 3, 3 ],
  [ 10, 6, 0, 5, 0 ],
  [ 10, 6, 0, 6, 10 ],
  [ 10, 6, 0, 13, 20 ],
  [ 1000, 1, 123, 456, "61894500" ],
  [ "100000000000000000000", 1234567, 1234567, 2469134,
    "5000000000000000000150000000000000000000" ],
);

my @overflow_tests = $use64 ?
  ( [6074001001, 1, 1, 0, "18446744077037500500", "triangular term"],
    [6074001000, 1, 2, 0, "36893488141926999000", "scaled triangular term"],
    [2, 1, 0, ~0, "36893488147419103230", "constant term"],
    [~0, ~0, 2, 0, "9223372036854775807", "a*n+b intermediate"],
    [1, 1, 0, ~0, "18446744073709551615", "exact UV_MAX result"],
  ) :
  (
    [92683, 1, 1, 0, "4295022903", "triangular term"],
    [92682, 1, 2, 0, "8589860442", "scaled triangular term"],
    [2, 1, 0, ~0, "8589934590", "constant term"],
    [~0, ~0, 2, 0, "2147483647", "a*n+b intermediate"],
    [1, 1, 0, ~0, "4294967295", "exact UV_MAX result"],
  );

plan tests => scalar(@tests) + scalar(@overflow_tests) + 4;

foreach my $t (@tests) {
  my($n,$m,$a,$b,$sum) = @$t;
  is("".floor_sum($n,$m,$a,$b), "$sum", "floor_sum($n,$m,$a,$b) = $sum");
}

foreach my $t (@overflow_tests) {
  my($n,$m,$a,$b,$sum,$desc) = @$t;
  is("".floor_sum($n,$m,$a,$b), "$sum", "floor_sum overflow: $desc");
}

ok(!eval { floor_sum(-1,7,3,4); 1 }, "floor_sum rejects negative n");
ok(!eval { floor_sum(5,0,3,4); 1 }, "floor_sum rejects zero m");
ok(!eval { floor_sum(5,7,-3,4); 1 }, "floor_sum rejects negative a");
ok(!eval { floor_sum(5,7,3,-4); 1 }, "floor_sum rejects negative b");
