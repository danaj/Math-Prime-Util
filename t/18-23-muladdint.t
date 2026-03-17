#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/muladdint mulsubint mulint/;

plan tests => 1   # muladdint small grid
            + 1   # mulsubint small grid
            + 1   # muladdint(n,m,0) == mulint(n,m)
            + 1   # mulsubint(n,m,0) == mulint(n,m)
            + 1   # selected exact values
            + 1   # product overflows UV but result fits
            + 1   # bigint
            ;

###### small grid

{
  my(@got, @exp);
  for my $a (-3 .. 3) {
    for my $b (-3 .. 3) {
      for my $c (-3 .. 3) {
        push @got, muladdint($a, $b, $c);
        push @exp, $a*$b + $c;
      }
    }
  }
  is_deeply(\@got, \@exp, "muladdint(-3..3, -3..3, -3..3)");
}

{
  my(@got, @exp);
  for my $a (-3 .. 3) {
    for my $b (-3 .. 3) {
      for my $c (-3 .. 3) {
        push @got, mulsubint($a, $b, $c);
        push @exp, $a*$b - $c;
      }
    }
  }
  is_deeply(\@got, \@exp, "mulsubint(-3..3, -3..3, -3..3)");
}

###### muladdint/mulsubint with c=0 matches mulint

{
  my @vals = (
    [qw/13282407956253574712 14991082624209354397/],
    [qw/4294967295 4294967296/],
    [qw/18446744073709551615 1/],
    [qw/9223372036854775807 -1/],
    [qw/-9223372036854775808 1/],
    [qw/-9223372036854775807 -1/],
    [qw/18446744073709551615 2/],
  );
  is_deeply(
    [map { "".muladdint($_->[0], $_->[1], 0) } @vals],
    [map { "".mulint(   $_->[0], $_->[1]   ) } @vals],
    "muladdint(a,b,0) == mulint(a,b)"
  );
  is_deeply(
    [map { "".mulsubint($_->[0], $_->[1], 0) } @vals],
    [map { "".mulint(   $_->[0], $_->[1]   ) } @vals],
    "mulsubint(a,b,0) == mulint(a,b)"
  );
}

###### selected exact values

subtest 'selected exact values', sub {
  # [a, b, c, muladdint=a*b+c, mulsubint=a*b-c]
  my @data = (
    # small signed
    [  3,  5,  7,   22,    8 ],
    [ -3,  5,  7,   -8,  -22 ],
    [  3, -5,  7,   -8,  -22 ],
    [  3,  5, -7,    8,   22 ],
    [ -3, -5,  7,   22,    8 ],
    [ -3, -5, -7,    8,   22 ],
    # mulsubint: a*b - (-c) = a*b + c
    [  6,  7, -1,   41,   43 ],
    # identity: a*b +/- 0
    [ 1000000, 1000000, 0, 1000000000000, 1000000000000 ],

    # IV_MAX boundary (product fits in UV as unsigned)
    # muladdint(IV_MAX, 1, 1) = IV_MAX + 1 = 9223372036854775808
    [qw/9223372036854775807  1  1  9223372036854775808  9223372036854775806/],
    # muladdint(IV_MAX, 2, 1) = 2*IV_MAX + 1 = UV_MAX
    [qw/9223372036854775807  2  1  18446744073709551615  18446744073709551613/],
    [qw/9223372036854775807  2 -1  18446744073709551613  18446744073709551615/],

    # muladdint(-IV_MAX, 1, -1) = -IV_MAX - 1 = IV_MIN = -9223372036854775808
    [qw/-9223372036854775807  1  -1  -9223372036854775808  -9223372036854775806/],

    # UV_MAX * 1 stays UV_MAX or shifts by one
    [qw/18446744073709551615  1   0  18446744073709551615  18446744073709551615/],
    [qw/18446744073709551615  1  -1  18446744073709551614  18446744073709551616/],
    [qw/18446744073709551615  1   1  18446744073709551616  18446744073709551614/],

    # muladdint(-1, UV_MAX, 1) = -UV_MAX + 1
    [qw/-1  18446744073709551615  1  -18446744073709551614  -18446744073709551616/],
  );

  for my $t (@data) {
    my($a, $b, $c, $exp_add, $exp_sub) = @$t;
    is("".muladdint($a,$b,$c), $exp_add, "muladdint($a,$b,$c)=$exp_add");
    is("".mulsubint($a,$b,$c), $exp_sub, "mulsubint($a,$b,$c)=$exp_sub");
  }
};

###### product overflows UV but final result fits (dispatches to PP/GMP)

subtest 'product overflow but result fits', sub {
  # 9223372036854775808 * 2 = 18446744073709551616 (> UV_MAX)
  # + (-9223372036854775808) = 9223372036854775808
  is("".muladdint("9223372036854775808", 2, "-9223372036854775808"),
     "9223372036854775808",
     "muladdint: product overflows, positive result");

  # mulsubint: same numeric result via subtraction
  is("".mulsubint("9223372036854775808", 2, "9223372036854775808"),
     "9223372036854775808",
     "mulsubint: product overflows, positive result");

  # 4294967296^2 = 18446744073709551616; minus that = 0
  is("".mulsubint("4294967296", "4294967296", "18446744073709551616"),
     "0",
     "mulsubint: product overflows, result 0");

  # UV_MAX * 2 - UV_MAX = UV_MAX
  is("".mulsubint("18446744073709551615", 2, "18446744073709551615"),
     "18446744073709551615",
     "mulsubint: 2*UV_MAX - UV_MAX = UV_MAX");

  # negative: -(UV_MAX)*2 + UV_MAX = -UV_MAX
  is("".muladdint("-18446744073709551615", 2, "18446744073709551615"),
     "-18446744073709551615",
     "muladdint: product overflows, negative result");
};

###### bigint inputs and results

subtest 'bigint', sub {
  my($a, $b) = qw/13282407956253574712 14991082624209354397/;
  my $prod    = "199117675120653046511338473800925208664";
  my $twoprod = "398235350241306093022676947601850417328";

  # muladdint(a, b, a*b) = 2*a*b
  is("".muladdint($a, $b, $prod),  $twoprod, "muladdint(a,b,a*b) = 2*a*b");
  # mulsubint(a, b, a*b) = 0
  is("".mulsubint($a, $b, $prod),  "0",      "mulsubint(a,b,a*b) = 0");

  # c dominates; product = 1
  is("".muladdint(1, 1, "999999999999999999999999999999"),
     "1000000000000000000000000000000",
     "muladdint(1,1,bigint_c)");
  is("".mulsubint("999999999999999999999999999999", 1, 1),
     "999999999999999999999999999998",
     "mulsubint(bigint,1,1)");

  # negative bigint result
  is("".mulsubint(1, 1, "999999999999999999999999999999"),
     "-999999999999999999999999999998",
     "mulsubint(1,1,bigint_c) negative");
};
