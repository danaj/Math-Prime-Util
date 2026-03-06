#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/
  addint subint mulint negint absint cmpint
  divint modint cdivint divrem fdivrem cdivrem tdivrem
  powint
  addmod submod mulmod powmod divmod muladdmod mulsubmod
/;

my $use64 = (~0 > 4294967295);
my $bits  = $use64 ? 64 : 32;
my $uvmax = $use64 ? "18446744073709551615" : "4294967295";
my $ivmax = $use64 ? "9223372036854775807" : "2147483647";
my $ivmin = $use64 ? "-9223372036854775808" : "-2147483648";
# just past UV boundary
my $uvmax1 = $use64 ? "18446744073709551616" : "4294967296";
my $ivmax1 = $use64 ? "9223372036854775808" : "2147483648";

plan tests =>
            + 1                # addint associativity
            + 1                # subint(IV_MIN, 1) crosses to bigint
            + 1                # addint(UV_MAX, 1) crosses to bigint
            + 1                # addint(IV_MIN, -1) crosses to bigint
            + 1                # subint(0, UV_MAX) = -UV_MAX
            + 1                # subint(IV_MIN, IV_MAX) extreme spread
            + 1                # mulint near-boundary products
            + 1                # mulint distributive
            + 1                # mulint(-1, IV_MIN)
            + 1                # powint special cases
            + 1                # powint bigint overflow
            + 1                # cmpint systematic
            + 1                # cmpint equal values
            + 1                # cmpint transitivity
            + 1                # cmpint sign dominance
            + 1                # divint/modint by 1 and -1
            + 1                # divint exact division
            + 1                # divint small / large
            + 1                # divrem invariant on boundary values
            + 1                # modops: addmod algebraic identity
            + 1                # modops: mulmod algebraic identity
            + 1                # modops: powmod(a,0,m) == 1 mod m
            + 1                # modops: powmod(a,1,m) == a mod m
            + 1                # modops: mulmod near UV_MAX
            + 1                # modops: muladdmod/mulsubmod consistency
            + 1                # modops: addmod/submod inverse
            + 0;


##### addint: associativity (a+b)+c == a+(b+c)
{
  my @triples = (
    [1, 2, 3],
    [$ivmax, 1, 1],
    [-1, $ivmin, -1],
    [$uvmax, $uvmax, $uvmax],
    ["$ivmin", "$ivmin", "$ivmin"],
    [$uvmax, "-$ivmax", 1],
    ["123456789012345678901234", "-987654321098765432109876", "864197532086419753208642"],
  );
  my $ok = 1;
  for my $t (@triples) {
    my($a,$b,$c) = @$t;
    $ok = 0 if "".addint(addint($a,$b),$c) ne "".addint($a,addint($b,$c));
  }
  ok($ok, "addint is associative");
}

##### addint/subint: boundary crossings
is("".subint($ivmin, 1), "".addint($ivmin, -1), "subint(IV_MIN, 1) == addint(IV_MIN, -1)");
is("".addint($uvmax, 1), $uvmax1, "addint(UV_MAX, 1) crosses into bigint");
{
  my $expect_below_ivmin = subint($ivmin, 1);
  # IV_MIN - 1 should be one less
  is("".addint($expect_below_ivmin, 1), $ivmin, "addint(IV_MIN-1, 1) returns to IV_MIN");
}

##### subint: extreme negation and spread
is("".subint(0, $uvmax), "-$uvmax", "subint(0, UV_MAX) = -UV_MAX");
{
  # IV_MIN - IV_MAX should be a very negative number (2*IV_MIN + 1 in magnitude)
  my $result = subint($ivmin, $ivmax);
  my $check  = addint($result, $ivmax);
  is("$check", $ivmin, "subint(IV_MIN, IV_MAX) + IV_MAX == IV_MIN");
}

##### mulint: products near UV boundary
{
  # (2^32-1)*(2^32-1) = 2^64 - 2^33 + 1  (overflows UV on 64-bit)
  # (2^16-1)*(2^16-1) = 2^32 - 2^17 + 1  (overflows UV on 32-bit)
  my @cases;
  if ($use64) {
    push @cases,
      ["4294967295", "4294967295", "18446744065119617025"],
      ["4294967296", "4294967295", "18446744069414584320"],
      [$uvmax, $uvmax, "340282366920938463426481119284349108225"];
  } else {
    push @cases,
      ["65535", "65535", "4294836225"],
      [$uvmax, $uvmax, "18446744065119617025"];
  }
  my $ok = 1;
  for my $c (@cases) {
    $ok = 0 if "".mulint($c->[0],$c->[1]) ne $c->[2];
  }
  ok($ok, "mulint products near UV boundary");
}

##### mulint: distributive property a*(b+c) == a*b + a*c
{
  my @cases = (
    [7, 11, 13],
    [-5, $ivmax, 1],
    [3, $uvmax, "-$ivmax"],
    ["123456789", "-987654321", "864197532"],
  );
  my $ok = 1;
  for my $c (@cases) {
    my($a,$b,$cc) = @$c;
    my $lhs = mulint($a, addint($b,$cc));
    my $rhs = addint(mulint($a,$b), mulint($a,$cc));
    $ok = 0 if "$lhs" ne "$rhs";
  }
  ok($ok, "mulint distributive: a*(b+c) == a*b + a*c");
}

##### mulint: -1 * IV_MIN (tricky: result is positive and > IV_MAX)
{
  is("".mulint(-1, $ivmin), $ivmax1, "mulint(-1, IV_MIN) = IV_MAX+1");
}

##### powint: special cases
{
  my $ok = 1;
  $ok = 0 if powint(0, 0) != 1;
  $ok = 0 if powint(0, 1) != 0;
  $ok = 0 if powint(0, 100) != 0;
  $ok = 0 if powint(1, 0) != 1;
  $ok = 0 if powint(1, 9999999) != 1;
  $ok = 0 if powint(-1, 0) != 1;
  $ok = 0 if powint(-1, 1) != -1;
  $ok = 0 if powint(-1, 2) != 1;
  $ok = 0 if powint(-1, 99) != -1;
  $ok = 0 if powint(-1, 100) != 1;
  $ok = 0 if powint(2, 0) != 1;
  $ok = 0 if powint(2, 1) != 2;
  $ok = 0 if powint(2, 10) != 1024;
  $ok = 0 if "".powint(-2, 3) ne "-8";
  $ok = 0 if "".powint(-2, 4) ne "16";
  ok($ok, "powint special cases: 0^k, 1^k, (-1)^k, 2^k, (-2)^k");
}

##### powint: overflow to bigint
{
  # 2^BITS should produce UV_MAX+1 as a bigint
  is("".powint(2, $bits), $uvmax1, "powint(2, BITS) produces bigint UV_MAX+1");
}

##### cmpint: systematic coverage
{
  my @vals = (0, 1, -1, $ivmax, $ivmin, $uvmax, "-$uvmax",
              $uvmax1, "-$uvmax1");
  my $ok = 1;
  for my $a (@vals) {
    for my $b (@vals) {
      my $cmp = cmpint($a, $b);
      # Verify: cmpint agrees with subint sign
      my $diff = subint($a, $b);
      my $expected = ($diff eq "0") ? 0 : (substr("$diff",0,1) eq '-') ? -1 : 1;
      if ($cmp != $expected) {
        $ok = 0;
      }
    }
  }
  ok($ok, "cmpint agrees with sign(subint(a,b)) for boundary values");
}

##### cmpint: equal values
{
  my @vals = (0, 1, -1, $ivmax, $ivmin, $uvmax, $uvmax1, "-$uvmax1",
              "340282366920938463463374607431768211456");
  my $ok = 1;
  for my $v (@vals) {
    $ok = 0 if cmpint($v, $v) != 0;
  }
  ok($ok, "cmpint(n,n) == 0 for all test values");
}

##### cmpint: transitivity
{
  my @sorted = ("-340282366920938463463374607431768211456",
                "-$uvmax1", "-$uvmax", $ivmin, -1, 0, 1, $ivmax,
                $uvmax, $uvmax1,
                "340282366920938463463374607431768211456");
  my $ok = 1;
  for my $i (0 .. $#sorted-1) {
    $ok = 0 if cmpint($sorted[$i], $sorted[$i+1]) >= 0;
    # Also verify reverse
    $ok = 0 if cmpint($sorted[$i+1], $sorted[$i]) <= 0;
  }
  ok($ok, "cmpint: sorted list has correct ordering (transitivity)");
}

##### cmpint: sign dominance - any negative < any non-negative
{
  my @neg = (-1, $ivmin, "-$uvmax", "-$uvmax1");
  my @pos = (0, 1, $ivmax, $uvmax, $uvmax1);
  my $ok = 1;
  for my $n (@neg) {
    for my $p (@pos) {
      $ok = 0 if cmpint($n, $p) >= 0;
    }
  }
  ok($ok, "cmpint: all negative values < all non-negative values");
}

##### divint/modint: divide by 1 and -1
{
  my @vals = (0, 1, -1, 7, -7, $ivmax, $ivmin, $uvmax, "-$uvmax",
              "39458349850349850394853049583049");
  my $ok = 1;
  for my $n (@vals) {
    # divint(n,1) == n, modint(n,1) == 0
    $ok = 0 if "".divint($n, 1) ne "$n";
    $ok = 0 if modint($n, 1) != 0;
    # divint(n,-1) == -n, modint(n,-1) == 0
    $ok = 0 if "".divint($n, -1) ne "".negint($n);
    $ok = 0 if modint($n, -1) != 0;
  }
  ok($ok, "divint(n,1)==n, modint(n,1)==0, divint(n,-1)==negint(n)");
}

##### divint: exact division (zero remainder)
{
  my @cases = (
    [0, 7],
    [42, 7],
    [-42, 7],
    [42, -7],
    [$uvmax, $uvmax],
    ["$ivmin", 1],
    ["39458349850349850394853049583049", "85889"],
  );
  my $ok = 1;
  for my $c (@cases) {
    my($prod, $d) = @$c;
    # Multiply first to guarantee exact division
    my $n = mulint($prod, $d);
    # Now divint(n, d) should give prod back (floor division)
    my($q, $r) = fdivrem($n, $d);
    $ok = 0 if "$r" ne "0";
    $ok = 0 if "$q" ne "".divint($n, $d);
  }
  ok($ok, "fdivrem gives zero remainder for exact multiples");
}

##### divint: small numerator / large denominator
{
  my @cases = (
    [0, $uvmax, 0],
    [1, $uvmax, 0],
    [-1, $uvmax, -1],
    [1, $uvmax1, 0],
    [$ivmax, $uvmax, 0],
  );
  my $ok = 1;
  for my $c (@cases) {
    $ok = 0 if "".divint($c->[0], $c->[1]) ne "$c->[2]";
  }
  ok($ok, "divint: small / large gives correct quotient");
}

##### divrem invariant on boundary values
{
  my @pairs = (
    [$uvmax, 2], [$uvmax, 3], [$uvmax, $ivmax],
    ["$ivmin", 3], ["$ivmin", -3], ["$ivmin", $ivmax],
    [$uvmax1, 7], ["-$uvmax1", 7],
    ["340282366920938463463374607431768211456", $uvmax],
  );
  my $ok = 1;
  for my $p (@pairs) {
    my($a, $b) = @$p;
    for my $fn (\&divrem, \&tdivrem, \&fdivrem, \&cdivrem) {
      my($q, $r) = $fn->($a, $b);
      # invariant: a == b*q + r
      $ok = 0 if "".addint(mulint($b, $q), $r) ne "$a";
    }
  }
  ok($ok, "a == b*q + r for all division modes on boundary values");
}


##### modops: addmod algebraic identity
#   addmod(a,b,m) == (a+b) mod m  when m > 0
{
  my @cases = (
    [0, 0, 7],
    [$uvmax, 1, $uvmax],
    [$uvmax, $uvmax, $uvmax],
    [$ivmax, $ivmax, $uvmax],
    ["340282366920938463463374607431768211456", $uvmax, "999999999999999989"],
  );
  my $ok = 1;
  for my $c (@cases) {
    my($a,$b,$m) = @$c;
    my $got = addmod($a, $b, $m);
    my $exp = modint(addint($a, $b), $m);
    $ok = 0 if "$got" ne "$exp";
  }
  ok($ok, "addmod(a,b,m) == modint(addint(a,b), m)");
}

##### modops: mulmod algebraic identity
#   mulmod(a,b,m) == (a*b) mod m  when m > 0
{
  my @cases = (
    [0, 0, 7],
    [0, $uvmax, 13],
    [$uvmax, 1, $uvmax],
    [$uvmax, 2, $uvmax],
    [$ivmax, $ivmax, $uvmax],
    [$uvmax, $uvmax, "999999999999999989"],
  );
  my $ok = 1;
  for my $c (@cases) {
    my($a,$b,$m) = @$c;
    my $got = mulmod($a, $b, $m);
    my $exp = modint(mulint($a, $b), $m);
    $ok = 0 if "$got" ne "$exp";
  }
  ok($ok, "mulmod(a,b,m) == modint(mulint(a,b), m)");
}

##### modops: powmod(a, 0, m) == 1 mod m (for m > 1)
{
  my @bases = (0, 1, 2, -1, $ivmax, $uvmax);
  my @mods  = (2, 3, 7, $uvmax, "999999999999999989");
  my $ok = 1;
  for my $a (@bases) {
    for my $m (@mods) {
      my $got = powmod($a, 0, $m);
      # a^0 mod m should be 1 mod m, which is 1 when m > 1
      $ok = 0 if "$got" ne "1";
    }
  }
  ok($ok, "powmod(a, 0, m) == 1 for m > 1");
}

##### modops: powmod(a, 1, m) == a mod m
{
  my @cases = (
    [0, 7],
    [1, 7],
    [6, 7],
    [7, 7],
    [-1, 7],
    [$uvmax, "999999999999999989"],
    [$ivmax, $uvmax],
  );
  my $ok = 1;
  for my $c (@cases) {
    my($a, $m) = @$c;
    my $got = powmod($a, 1, $m);
    my $exp = modint($a, $m);
    $ok = 0 if "$got" ne "$exp";
  }
  ok($ok, "powmod(a, 1, m) == a mod m");
}

##### modops: mulmod near UV_MAX (exercises 128-bit intermediate)
{
  my @cases;
  if ($use64) {
    @cases = (
      # (UV_MAX-1) * (UV_MAX-1) mod UV_MAX should be 1
      ["18446744073709551614", "18446744073709551614", $uvmax, 1],
      # UV_MAX * UV_MAX mod (UV_MAX - 1) should be 1
      [$uvmax, $uvmax, "18446744073709551614", 1],
      # large prime modulus
      [$uvmax, $uvmax, "999999999999999989", "587155414672247084"],
    );
  } else {
    @cases = (
      ["4294967294", "4294967294", $uvmax, 1],
      [$uvmax, $uvmax, "4294967294", 1],
      [$uvmax, $uvmax, "999999937", "264566326"],
    );
  }
  my $ok = 1;
  for my $c (@cases) {
    my($a, $b, $m, $exp) = @$c;
    my $got = mulmod($a, $b, $m);
    $ok = 0 if "$got" ne "$exp";
  }
  ok($ok, "mulmod near UV_MAX boundary");
}

##### modops: muladdmod and mulsubmod consistency
#   muladdmod(a,b,c,m) == addmod(mulmod(a,b,m), c, m)
#   mulsubmod(a,b,c,m) == submod(mulmod(a,b,m), c, m)
{
  my @cases = (
    [0, 0, 0, 7],
    [3, 5, 7, 11],
    [$ivmax, 2, $ivmax, $uvmax],
    [$uvmax, $uvmax, $uvmax, "999999999999999989"],
    ["123456789012345678901234", "987654321098765432109876", "111111111111111111111111", "314159265358979323846263"],
  );
  my $ok = 1;
  for my $c (@cases) {
    my($a,$b,$cc,$m) = @$c;
    my $mam = muladdmod($a, $b, $cc, $m);
    my $msm = mulsubmod($a, $b, $cc, $m);
    my $mam_exp = addmod(mulmod($a, $b, $m), $cc, $m);
    my $msm_exp = submod(mulmod($a, $b, $m), $cc, $m);
    $ok = 0 if "$mam" ne "$mam_exp";
    $ok = 0 if "$msm" ne "$msm_exp";
  }
  ok($ok, "muladdmod/mulsubmod == addmod/submod(mulmod(...))");
}

##### modops: addmod and submod are inverses
#   submod(addmod(a, b, m), b, m) == a mod m
{
  my @cases = (
    [0, 0, 7],
    [3, 5, 7],
    [6, 6, 7],
    [$ivmax, 1, $uvmax],
    [$uvmax, $uvmax, "999999999999999989"],
    [0, $uvmax, $uvmax],
  );
  my $ok = 1;
  for my $c (@cases) {
    my($a, $b, $m) = @$c;
    my $a_mod_m = modint($a, $m);
    my $roundtrip = submod(addmod($a, $b, $m), $b, $m);
    $ok = 0 if "$roundtrip" ne "$a_mod_m";
  }
  ok($ok, "submod(addmod(a,b,m), b, m) == a mod m");
}
