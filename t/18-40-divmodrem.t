#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/divint modint cdivint divrem fdivrem cdivrem tdivrem
                         addint mulint divint modint/;
use Math::BigInt;


my @quotients = (  # trunc, floor, ceiling, euclidian
  ["S + +", qw/ 9949242744253247 64  155456917878956  155456917878956  155456917878957  155456917878956/],
  ["S - +", qw/-9949242744253247 64 -155456917878956 -155456917878957 -155456917878956 -155456917878957/],

  ["L + +", qw/ 39458349850349850394853049583049  85889  459410982202026457344398579  459410982202026457344398579  459410982202026457344398580  459410982202026457344398579/],
  ["L + -", qw/ 39458349850349850394853049583049 -85889 -459410982202026457344398579 -459410982202026457344398580 -459410982202026457344398579 -459410982202026457344398579/],
  ["L - +", qw/-39458349850349850394853049583049  85889 -459410982202026457344398579 -459410982202026457344398580 -459410982202026457344398579 -459410982202026457344398580/],
  ["L - -", qw/-39458349850349850394853049583049 -85889  459410982202026457344398579  459410982202026457344398579  459410982202026457344398580  459410982202026457344398580/],
);

plan tests => 0
            + 1                    # divbyzero
            + 2                    # divint
            + 2 + 1                # modint
            + 12                   # table 1.3 from Leijen 2001
            + 11                   # divint with large neg returns
            + 3                    # cdivrem extra
            + scalar(@quotients)   # signed bigint divint+
            + scalar(@quotients)   # signed bigint divrem+
            + 1                    # check b*q+r=a
            + 4                    # verify a=b*q+r
            + 1                    # Euclidean remainder in range
            + 1                    # cdivrem remainder sign
            + 1                    # fdivrem = divint/modint
            + 0;

###### divide by 0 should error
{
  my $ok = 0;
  $ok++ if !defined eval {  divint(0,0); } && $@ =~ /divide by zero/;
  $ok++ if !defined eval {  divint(1,0); } && $@ =~ /divide by zero/;
  $ok++ if !defined eval { cdivint(0,0); } && $@ =~ /divide by zero/;
  $ok++ if !defined eval { cdivint(1,0); } && $@ =~ /divide by zero/;
  $ok++ if !defined eval {  modint(0,0); } && $@ =~ /divide by zero/;
  $ok++ if !defined eval {  modint(1,0); } && $@ =~ /divide by zero/;
  $ok++ if !defined eval {  divrem(0,0); } && $@ =~ /divide by zero/;
  $ok++ if !defined eval {  divrem(1,0); } && $@ =~ /divide by zero/;
  $ok++ if !defined eval { fdivrem(0,0); } && $@ =~ /divide by zero/;
  $ok++ if !defined eval { fdivrem(1,0); } && $@ =~ /divide by zero/;
  $ok++ if !defined eval { cdivrem(0,0); } && $@ =~ /divide by zero/;
  $ok++ if !defined eval { cdivrem(1,0); } && $@ =~ /divide by zero/;
  $ok++ if !defined eval { tdivrem(0,0); } && $@ =~ /divide by zero/;
  $ok++ if !defined eval { tdivrem(1,0); } && $@ =~ /divide by zero/;

  is($ok, 14, "divide by zero correctly trapped");
}

# For negative inputs, the div and mod operations might be different than Perl's builtins.
# It matches Math::BigInt bdiv / bmod (post 1.997 Sep 2015).

my @qpos1024 = map { int(1024/$_) } 1 .. 1025;
my @qneg1024 = map { my $d=-1024/$_; my $i = int($d);  ($d==$i) ? $i : $i-1; } 1 .. 1025;

my @rpos1024 = map {  1024 - $_ * $qpos1024[$_-1] } 1 .. 1025;
my @rneg1024 = map { -1024 - $_ * $qneg1024[$_-1] } 1 .. 1025;

is_deeply( [map { divint(1024,$_) } 1..1025], \@qpos1024, "divint(1024,x) for 1 .. 1025" );
is_deeply( [map { divint(-1024,$_) } 1..1025], \@qneg1024, "divint(-1024,x) for 1 .. 1025" );

###### modint
is_deeply( [map { modint(1024,$_) } 1..1025], \@rpos1024, "modint(1024,x) for 1 .. 1025" );
is_deeply( [map { modint(-1024,$_) } 1..1025], \@rneg1024, "modint(-1024,x) for 1 .. 1025" );

is(modint("-1117091728166568014",59), 4, "modint(-1117091728166568014,59) = 4");

###### small values through divint, modint, divrem, fdivrem, cdivrem, tdivrem
is( join(" ", tdivrem(8,3), tdivrem(8,-3), tdivrem(-8,3), tdivrem(-8,-3)),
    "2 2 -2 2 -2 -2 2 -2",
    "tdivrem with +/- 8,3" );
is( join(" ", divrem(8,3), divrem(8,-3), divrem(-8,3), divrem(-8,-3)),
    "2 2 -2 2 -3 1 3 1",
    "divrem with +/- 8,3" );
is( join(" ", fdivrem(8,3), fdivrem(8,-3), fdivrem(-8,3), fdivrem(-8,-3)),
    "2 2 -3 -1 -3 1 2 -2",
    "fdivrem with +/- 8,3" );
is( join(" ", cdivrem(8,3), cdivrem(8,-3), cdivrem(-8,3), cdivrem(-8,-3)),
    "3 -1 -2 2 -2 -2 3 1",
    "cdivrem with +/- 8,3" );
is( join(" ", divint(8,3), modint(8,3), divint(8,-3), modint(8,-3), divint(-8,3), modint(-8,3), divint(-8,-3), modint(-8,-3)),
    "2 2 -3 -1 -3 1 2 -2",
    "divint+modint with +/- 8,3" );
is( join(" ", cdivint(8,3),cdivint(8,-3),cdivint(-8,3),cdivint(-8,-3)),
    "3 -2 -2 3",
    "cdivint with +/- 8,3" );

is( join(" ", tdivrem(1,2), tdivrem(1,-2), tdivrem(-1,2), tdivrem(-1,-2)),
    "0 1 0 1 0 -1 0 -1",
    "tdivrem with +/- 1,2" );
is( join(" ", divrem(1,2), divrem(1,-2), divrem(-1,2), divrem(-1,-2)),
    "0 1 0 1 -1 1 1 1",
    "divrem with +/- 1,2" );
is( join(" ", fdivrem(1,2), fdivrem(1,-2), fdivrem(-1,2), fdivrem(-1,-2)),
    "0 1 -1 -1 -1 1 0 -1",
    "fdivrem with +/- 1,2" );
is( join(" ", cdivrem(1,2), cdivrem(1,-2), cdivrem(-1,2), cdivrem(-1,-2)),
    "1 -1 0 1 0 -1 1 1",
    "cdivrem with +/- 1,2" );
is( join(" ", divint(1,2), modint(1,2), divint(1,-2), modint(1,-2), divint(-1,2), modint(-1,2), divint(-1,-2), modint(-1,-2)),
    "0 1 -1 -1 -1 1 0 -1",
    "divint+modint with +/- 1,2" );
is( join(" ", cdivint(1,2),cdivint(1,-2),cdivint(-1,2),cdivint(-1,-2)),
    "1 0 0 1",
    "cdivint with +/- 1,2" );

###### divint and modint with interesting values
is("".divint("1895315831", -1), "-1895315831", "Divide 31-bit input by -1");
is("".divint("3483637757", -1), "-3483637757", "Divide 32-bit input by -1");
is("".cdivint("3483637757", -1), "-3483637757", "Divide 32-bit input by -1 (ceiling)");
is("".divint("6127303089832103323", -1), "-6127303089832103323", "Divide 63-bit input by -1");
is("".divint("13026328650942325963", -1), "-13026328650942325963", "Divide 64-bit input by -1");
is("".divint("14123555781055773270", 2), "7061777890527886635", "Divide 64-bit input by 2");
is("".divint("12844039487317506779", "12844039487317506779"), 1, "Divide 64-bit input by itself");
is(divint(3, "12844039487317506779"), 0, "Divide small int by 64-bit input");
# Note this is floor division:
is(divint(-3, "12844039487317506779"), -1, "Divide negative small int by 64-bit input");
# Now ceiling
is(cdivint(3, "12844039487317506779"), 1, "Divide (ceil) small int by 64-bit input");
is(cdivint(-3, "12844039487317506779"), 0, "Divide (ceil) negative small int by 64-bit input");

###### cdivrem special test
is( join(" ",cdivrem(3, "12844039487317506779")), "1 -12844039487317506776", "cdivrem with small quotient and 64-bit denominator shouldn't overflow IV" );
{
 my $x248 = "281474976710656";
 is_deeply([cdivint("9223372036854775807",$x248),   # 2^63-1 / 2^48
            cdivint("9223372036854775808",$x248),   # 2^63   / 2^48
            cdivint("9223372036854775809",$x248)],  # 2^63+1 / 2^48
           [32768, 32768, 32769],
           "cdivint (2^63 +/- 1) / 2^48");
 is_deeply([cdivint("18446744073709551615",$x248),  # 2^64-1 / 2^48
            cdivint("18446744073709551616",$x248),  # 2^64   / 2^48
            cdivint("18446744073709551617",$x248)], # 2^64+1 / 2^48
           [65536, 65536, 65537],
           "cdivint (2^64 +/- 1) / 2^48");
}

###### large values through divint, cdivint, modint,
######                      divrem, tdivrem, fdivrem, cdivrem
for my $s (@quotients) {
  my($signs, $n, $m, $qt, $qf, $qc, $qe) = @$s;
  my($bn,$bm) = map { Math::BigInt->new($_) } ($n,$m);
  my($rt, $rf, $rc, $re) = map{"$_"}map { $bn - $bm * $_ } ($qt, $qf, $qc, $qe);

  #is( "".divint($n, $m), $qf, "large divint  $signs" );
  #is( "".modint($n, $m), $rf, "large modint  $signs" );
  #is( "".cdivint($n, $m), $qc, "large cdivint  $signs" );
  #is_deeply( [map{"$_"}divrem($n, $m)], [$qe, $re], "large divrem  $signs" );
  #is_deeply( [map{"$_"}tdivrem($n, $m)], [$qt, $rt], "large tdivrem $signs" );
  #is_deeply( [map{"$_"}fdivrem($n, $m)], [$qf, $rf], "large fdivrem $signs" );
  #is_deeply( [map{"$_"}cdivrem($n, $m)], [$qc, $rc], "large cdivrem $signs" );

  is_deeply( ["".divint($n,$m), "".cdivint($n,$m), "".modint($n,$m)],
             [   $qf,              $qc,               $rf,         ],
             "$signs   divint, cdivint, modint" );

  is_deeply( [[map{"$_"} divrem($n,$m)], [map{"$_"}tdivrem($n,$m)],
              [map{"$_"}fdivrem($n,$m)], [map{"$_"}cdivrem($n,$m)]],
             [[$qe,$re], [$qt,$rt], [$qf,$rf], [$qc,$rc]],
             "$signs   divrem, tdivrem, fdivrem, cdivrem" );
}

# --- divint/modint invariant:  b*q + r == a ---
{
  my @pairs = (
    [13, 4], [-13, 4], [13, -4], [-13, -4],
    ["9949242744253247", 64], ["-9949242744253247", 64],
    ["18446744073709551615", 7], ["-18446744073709551617", 13],
    ["39458349850349850394853049583049", 85889],
    ["-39458349850349850394853049583049", 85889],
  );
  my $ok = 1;
  for my $p (@pairs) {
    my($a,$b) = @$p;
    $ok = 0 unless "".addint(mulint($b,divint($a,$b)),modint($a,$b)) eq "$a";
  }
  ok($ok, "b*divint(a,b) + modint(a,b) == a for all test pairs");
}

{
  my @pairs = ([8,3],[-8,3],[8,-3],[-8,-3],[1,2],[-1,2],[7,1],[-7,1]);
  my @names = qw/divrem tdivrem fdivrem cdivrem/;
  my @fns   = (\&divrem, \&tdivrem, \&fdivrem, \&cdivrem);
  for my $i (0..$#fns) {
    my $ok = 1;
    for my $p (@pairs) {
      my($a,$b) = @$p;
      my($q,$r) = $fns[$i]->($a,$b);
      $ok = 0 unless $b*$q + $r == $a;
    }
    ok($ok, "$names[$i]: a == b*q + r for all sign combinations");
  }
}

# --- Euclidean remainder is always >= 0 ---
{
  my @pairs = ([8,3],[-8,3],[8,-3],[-8,-3],
               ["39458349850349850394853049583049",85889],
               ["-39458349850349850394853049583049",85889]);
  my $ok = 1;
  for my $p (@pairs) {
    my($q,$r) = divrem($p->[0], $p->[1]);
    $ok = 0 if $r < 0;
  }
  ok($ok, "divrem: remainder is always >= 0");
}

# --- cdivrem: remainder has opposite sign from divisor ---
{
  my @cases = ([8,3],[8,-3],[-8,3],[-8,-3]);
  my $ok = 1;
  for my $c (@cases) {
    my($a,$b) = @$c;
    my($q,$r) = cdivrem($a,$b);
    # r == 0 or sign(r) != sign(b)
    $ok = 0 if $r != 0 && (($r > 0) == ($b > 0));
  }
  ok($ok, "cdivrem: remainder has opposite sign from divisor (or is zero)");
}

# --- fdivrem results match divint + modint ---
{
  my @pairs = ([13,4],[-13,4],[13,-4],[-13,-4],
               ["18446744073709551617",7],["-18446744073709551617",13]);
  my $ok = 1;
  for my $p (@pairs) {
    my($a,$b) = @$p;
    my($q,$r) = fdivrem($a,$b);
    $ok = 0 unless "$q" eq "".divint($a,$b) && "$r" eq "".modint($a,$b);
  }
  ok($ok, "fdivrem results match divint+modint");
}
