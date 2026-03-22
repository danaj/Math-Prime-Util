#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/toint powint/;
use Math::BigInt;
use Math::BigFloat;

my $use64 = (~0 > 4294967295);

sub is_native { !ref($_[0]) }
sub is_bigint {  ref($_[0]) }

plan tests => 1  # native ints
            + 1  # small floats
            + 1  # integer strings
            + 1  # float strings
            + 1  # scientific notation strings
            + 1  # bigints
            + 1  # bigfloat
            + 1  # undef and empty strings
            + 1  # invalid
            + 1  # check return types
            ;

###############################################################################

subtest 'native integers pass through unchanged' => sub {
  is(toint(0),   0,   "zero");
  is(toint(1),   1,   "one");
  is(toint(-1), -1,   "negative one");
  is(toint(42),  42,  "small positive");
  is(toint(-42), -42, "small negative");
  is(toint(4294967295), 4294967295, "UV32_MAX");
  is(toint(4294967296), 4294967296, "UV32_MAX + 1");
  SKIP: {
    skip "64-bit only", 2 unless $use64;
    is(toint(9223372036854775807),  9223372036854775807,  "IV64_MAX");
    is(toint(18446744073709551615), 18446744073709551615, "UV64_MAX");
  }
  ok( is_native(toint(42)),  "native int result is not a ref");
  ok( is_native(toint(-42)), "native neg result is not a ref");
};

###############################################################################

subtest 'NV floats truncate toward zero (not floor)' => sub {
  is(toint(3.9),    3,  "3.9 → 3");
  is(toint(-3.9),  -3,  "-3.9 → -3");
  is(toint(0.9),    0,  "0.9 → 0");
  is(toint(-0.9),   0,  "-0.9 → 0  (truncate not floor; floor would be -1)");
  is(toint(1.0),    1,  "1.0 → 1");
  is(toint(-1.0),  -1,  "-1.0 → -1");
  is(toint(0.0),    0,  "0.0 → 0");
  is(toint(0.5),    0,  "0.5 → 0");
  is(toint(-0.5),   0,  "-0.5 → 0  (truncate not floor; floor would be -1)");
  is(toint(99.99), 99,  "99.99 → 99");
  ok( is_native(toint(3.9)), "float result is native");
};

###############################################################################

subtest 'integer strings' => sub {
  is(toint("0"),     0,   "'0'");
  is(toint("42"),   42,   "'42'");
  is(toint("-42"), -42,   "'-42'");
  is(toint("+42"),  42,   "'+42' explicit plus");
  is(toint("007"),   7,   "'007' leading zeros stripped");
  is(toint("-007"), -7,   "'-007' leading zeros stripped");
  is(toint("00"),    0,   "'00' is zero");
  ok( is_native(toint("42")),  "string int → native");
  ok( is_native(toint("-42")), "string neg int → native");
  # 20-digit integer > UV64_MAX → bigint
  { my $r = toint("99999999999999999999");
    ok( is_bigint($r), "20-digit integer string → bigint");
    is("$r", "99999999999999999999", "20-digit value correct");
  }
};

###############################################################################

subtest 'float strings — truncate toward zero' => sub {
  is(toint("3.14"),    3,  "'3.14' → 3");
  is(toint("-3.14"),  -3,  "'-3.14' → -3");
  is(toint("3.9"),     3,  "'3.9' → 3");
  is(toint("-3.9"),   -3,  "'-3.9' → -3");
  is(toint("3."),      3,  "'3.' trailing dot → 3");
  is(toint("-3."),    -3,  "'-3.' trailing dot → -3");
  is(toint(".5"),      0,  "'.5' no leading digit → 0");
  is(toint("-.5"),     0,  "'-.5' → 0  (truncate not floor; floor would be -1)");
  is(toint("+.5"),     0,  "'+.5' → 0");
  is(toint("0.9"),     0,  "'0.9' → 0");
  is(toint("-0.9"),    0,  "'-0.9' → 0  (truncate not floor)");
  is(toint("0.0"),     0,  "'0.0' → 0");
  is(toint("99.99"), 99,   "'99.99' → 99");
};

###############################################################################

subtest 'scientific notation strings' => sub {
  is(toint("1e3"),      1000,  "'1e3'");
  is(toint("1E3"),      1000,  "'1E3' uppercase E");
  is(toint("1.5e2"),    150,   "'1.5e2'");
  is(toint("1.5e+2"),   150,   "'1.5e+2' explicit positive exponent");
  is(toint("3e0"),        3,   "'3e0'");
  is(toint("0e5"),        0,   "'0e5' = 0");
  is(toint("1.5e-1"),     0,   "'1.5e-1' = 0.15 → 0 (truncate)");
  is(toint("-1.5e-1"),    0,   "'-1.5e-1' = -0.15 → 0 (truncate, not -1)");
  SKIP: { skip "64-bit only", 1 unless $use64;
    is(toint("1e10"),  10000000000, "'1e10'");
  }
  # large sci notation → bigint
  { my $r = toint("1e50");
    ok( is_bigint($r), "'1e50' → bigint");
    is( length("$r"), 51, "'1e50': 51 chars (1 followed by 50 zeros)");
    is( substr("$r",0,1), "1", "'1e50': starts with 1");
  }
};

###############################################################################

subtest 'BigInt objects' => sub {
  is(toint(Math::BigInt->new(0)),    0,    "BigInt(0)");
  is(toint(Math::BigInt->new(42)),  42,    "BigInt(42)");
  is(toint(Math::BigInt->new(-42)), -42,   "BigInt(-42)");
  ok( is_native(toint(Math::BigInt->new(99))),  "BigInt in native range → native");
  ok( is_native(toint(Math::BigInt->new(-99))), "negative BigInt in native range → native");

  {
    my $big = powint(2, 64);
    my $r = toint($big);
    ok( is_bigint($r), "2^64 → bigint");
    is("$r", "18446744073709551616", "2^64 value correct");

    my $neg = Math::BigInt->new("-9223372036854775809");
    my $rn = toint($neg);
    ok( is_bigint($rn), "-(2^63+1) → bigint");
    is("$rn", "-9223372036854775809", "-(2^63+1) value correct");
  }

  # Already an integer — no truncation needed even for large values
  { my $big = powint(2, 128);
    my $r = toint($big);
    ok( is_bigint($r), "BigInt(2^128) result is a ref");
    is("$r", "340282366920938463463374607431768211456", "BigInt(2^128) value correct");
  }
  { my $r = toint(Math::BigInt->new("123456789012345678901234567890"));
    is("$r", "123456789012345678901234567890", "large BigInt preserved exactly");
  }
};

###############################################################################

subtest 'BigFloat objects — truncate toward zero' => sub {
  is(toint(Math::BigFloat->new("3.14")),    3,  "BigFloat(3.14) → 3");
  is(toint(Math::BigFloat->new("-3.14")),  -3,  "BigFloat(-3.14) → -3");
  is(toint(Math::BigFloat->new("0.9")),     0,  "BigFloat(0.9) → 0");
  is(toint(Math::BigFloat->new("-0.9")),    0,  "BigFloat(-0.9) → 0 (truncate not floor)");
  is(toint(Math::BigFloat->new("3.0")),     3,  "BigFloat(3.0)");
  is(toint(Math::BigFloat->new("0")),       0,  "BigFloat(0)");
  ok( is_native(toint(Math::BigFloat->new("42"))),  "BigFloat(42) → native");
  ok( is_native(toint(Math::BigFloat->new("-42"))), "BigFloat(-42) → native");
  # Large BigFloat — truncate and return bigint
  { my $r = toint(Math::BigFloat->new("1.5e30"));
    ok( is_bigint($r), "BigFloat(1.5e30) → bigint");
    is("$r", "1500000000000000000000000000000", "BigFloat(1.5e30) truncated correctly");
  }
};

###############################################################################

subtest 'undef and empty string return 0' => sub {
  is(toint(undef), 0, "toint(undef) = 0");
  is(toint(""),    0, "toint('') = 0");
  ok( is_native(toint(undef)), "undef result is native");
  ok( is_native(toint("")),    "empty string result is native");
};

###############################################################################

subtest 'invalid inputs croak' => sub {
  my @bad = (
    "abc",       # pure alphabetic
    "hello",     # pure alphabetic
    "3..14",     # double decimal point
    "e5",        # exponent without mantissa
    "3e",        # incomplete exponent
    ".",         # lone decimal point
    "+-3",       # double sign
    "1 2",       # embedded space
    #"3 ",        # trailing space
    #" 3",        # leading space
    #"0xFF",      # hex — use hex() first
    #"0b101",     # binary — use oct() first
    "inf",       # infinity string
    "Inf",
    "nan",       # not-a-number string
    "NaN",
  );
  for my $bad (@bad) {
    local $@;
    eval { toint($bad) };
    ok($@, "toint('$bad') croaks") or diag "  toint('$bad') returned without error";
  }
};

###############################################################################

# Error at compile time
#subtest 'toint() with no arguments croaks' => sub {
#  local $@;
#  eval { toint() };
#  ok($@, "toint() no args croaks");
#};

###############################################################################

subtest 'return type: native when fits, bigint when large' => sub {
  my @native_cases = (
    [0,                         "zero NV"],
    [42,                        "native int"],
    [-42,                       "native neg int"],
    [3.9,                       "float NV"],
    ["99",                      "integer string"],
    ["3.7",                     "float string"],
    [Math::BigInt->new(42),     "BigInt in range"],
    [Math::BigFloat->new("3.14"), "BigFloat in range"],
  );
  for my $case (@native_cases) {
    my($v, $label) = @$case;
    ok( is_native(toint($v)), "toint($label) is native");
  }
  ok( is_native(toint(undef)), "toint(undef) is native");

  my @bigint_cases = (
    [powint(2,128),                              "2^128"],
    [Math::BigInt->new("99999999999999999999999"), "large BigInt"],
    [Math::BigFloat->new("1.5e30"),              "large BigFloat"],
  );
  for my $case (@bigint_cases) {
    my($v, $label) = @$case;
    ok( is_bigint(toint($v)), "toint($label) is bigint");
  }
};
