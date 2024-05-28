#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/ powint mulint addint subint
                          add1int sub1int
                          divint modint cdivint
                          divrem fdivrem cdivrem tdivrem
                          lshiftint rshiftint rashiftint
                          absint negint cmpint signint/;
use Math::BigInt;


my $use64 = (~0 > 4294967296 && 18446744073709550592 != ~0);
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
#my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
#my $usegmp= Math::Prime::Util::prime_get_config->{'gmp'};

my @powints = (
 [5, 6, 15625],
 [2, 16, 65536],
);
my @mulints = (
  [qw/13282407956253574712 14991082624209354397 199117675120653046511338473800925208664/],
  [qw/65536 65536 4294967296/],
  [qw/4294967296 4294967296 18446744073709551616/],
);
my @addints = (
  [qw/1178630961471601951655862 827639478068904540012 1179458600949670856195874/],
  [qw/-2555488174170453670799 1726145541361106236340 -829342632809347434459/],
  [qw/9223372036854775808 9223372036854775808 18446744073709551616/],
);
my @subints = (
  [qw/1178630961471601951655862 827639478068904540012 1177803321993533047115850/],
  [qw/-2555488174170453670799 1726145541361106236340 -4281633715531559907139/],
  [qw/9223372036854775808 -9223372036854775808 18446744073709551616/],
);
my @anvals = (
  "1178630961471601951655862","827639478068904540012",
  "-2555488174170453670799","1726145541361106236340",
  "-829342632809347434459","-4281633715531559907139",
  "-230948092384903284908329048239084023984092384",
  "982349082340982348502392937523840234029384908325098234",
  "32767","-32767","32768","-32768",
  "4294967295","-4294967295","4294967296","-4294967296",
  "9223372036854775807","-9223372036854775807",
  "9223372036854775808","-9223372036854775808",
  "18446744073709551614","-18446744073709551614",
  "18446744073709551615","-18446744073709551615",
  "18446744073709551616","-18446744073709551616",
  "18446744073709551617","-18446744073709551617",
);
my @quotients = (  # trunc, floor, ceiling, euclidian
  ["S + +", "9949242744253247", "64", "155456917878956", "155456917878956", "155456917878957", "155456917878956"],
  ["S - +", "-9949242744253247", "64", "-155456917878956", "-155456917878957", "-155456917878956", "-155456917878957"],
  ["L + +",  "39458349850349850394853049583049",  "85889",  "459410982202026457344398579",  "459410982202026457344398579",  "459410982202026457344398580",  "459410982202026457344398579"],
  ["L + -",  "39458349850349850394853049583049", "-85889", "-459410982202026457344398579", "-459410982202026457344398580", "-459410982202026457344398579", "-459410982202026457344398579"],
  ["L - +", "-39458349850349850394853049583049",  "85889", "-459410982202026457344398579", "-459410982202026457344398580", "-459410982202026457344398579", "-459410982202026457344398580"],
  ["L - -", "-39458349850349850394853049583049", "-85889",  "459410982202026457344398579",  "459410982202026457344398579",  "459410982202026457344398580",  "459410982202026457344398580"],
);
my @negshifts = (
  # n, k,  >>, >>arith
  [ 0, 1,  0, 0],
  [-1, 1,  0, -1],
  [-5, 1,  -2, -3],
  [-8, 2,  -2, -2],
  ["-307385513", 6, -4802898, -4802899],
  ["-637526413", 6, -9961350, -9961351],
  ["-2045651239", 6, -31963300, -31963301],
  ["-3675663743", 6, -57432245, -57432246],
  ["-2332267979728172537", 6, "-36441687183252695", "-36441687183252696"],
  ["-8408654401686460807", 6, "-131385225026350950", "-131385225026350951"],
  ["-17640827963513397449", 6, "-275637936929896835", "-275637936929896836"],
  ["-32659506018295865747", 6, "-510304781535872902", "-510304781535872903"],
  ["-79231600218559026832557301750107210001", 6, "-1237993753414984794258707839845425156", "-1237993753414984794258707839845425157"],
  ["-131954888069700539887213633881194728277", 6, "-2061795126089070935737713029393667629", "-2061795126089070935737713029393667630"],
  ["-254262665582332530470619504253273698569", 6, "-3972854149723945788603429753957401540", "-3972854149723945788603429753957401541"],
  ["-416649423645764932216789232242651032187", 6, "-6510147244465077065887331753791422377", "-6510147244465077065887331753791422378"],
);


plan tests => 0
            + 7*4 + scalar(@powints) + 5     # powint
            + 1 + scalar(@mulints)           # mulint
            + 1 + scalar(@addints)           # addint
            + 1 + scalar(@subints)           # subint
            + 1                              # add1int
            + 1                              # sub1int
            + 2 + 2                          # divint
            + 2 + 2 + 1                      # modint
            + 2                              # divrem
            + 2                              # fdivrem
            + 2                              # cdivrem
            + 2                              # tdivrem
            + 7 * scalar(@quotients)         # signed bigint division
            + 12                             # table 1.3 from Leijen 2001
            + 11                             # divint with large neg returns
            + 3                              # cdivrem extra
            + 4 + 3*scalar(@negshifts)       # shiftint
            + 4                              # absint
            + 6                              # negint
            + 5                              # cmpint
            + 3                              # signint
            + 0;

###### powint
for my $a (-3 .. 3) {
  is(powint($a, 0), 1, "powint($a,0) = 1");
  is(powint($a, 1), $a, "powint($a,1) = $a");
  is(powint($a, 2), $a*$a, "powint($a,2) = " . $a*$a);
  is(powint($a, 3), $a*$a*$a, "powint($a,3) = " . $a*$a*$a);
}
foreach my $r (@powints) {
  my($a, $b, $exp) = @$r;
  is( powint($a,$b), $exp, "powint($a,$b) = ".((defined $exp)?$exp:"<undef>") );
}
is(powint(powint(2,32),3),"79228162514264337593543950336","(2^32)^3");
is(powint(3,powint(2,7)),"11790184577738583171520872861412518665678211592275841109096961","3^(2^7)");
is(powint(46,22)+1, "3807783932766699862493193563344470017", "powint returns a bigint if needed");
is(powint(-544,7),"-14099129446552305664","powint(-544,7) = -14099129446552305664");
is(powint(-544,7)-1,"-14099129446552305665","powint(-544,7)-1 = -14099129446552305665");

###### mulint
{ my(@got,@exp);
  for my $a (-3 .. 3) {
    for my $b (-3 .. 3) {
      push @got, mulint($a,$b);
      push @exp, ($a == 0 || $b == 0) ? 0 : $a*$b;  # Perl 5.6: -1*0 = -0
    }
  }
  is_deeply( \@got, \@exp, "mulint( -3 .. 3, -3 .. 3)" );
}
foreach my $r (@mulints) {
  my($a, $b, $exp) = @$r;
  is( mulint($a,$b), $exp, "mulint($a,$b) = ".((defined $exp)?$exp:"<undef>") );
}

###### addint
{ my(@got,@exp);
  for my $a (-3 .. 3) {
    for my $b (-3 .. 3) {
      push @got, addint($a,$b);
      push @exp, $a+$b;
    }
  }
  is_deeply( \@got, \@exp, "addint( -3 .. 3, -3 .. 3)" );
}
foreach my $r (@addints) {
  my($a, $b, $exp) = @$r;
  is( addint($a,$b), $exp, "addint($a,$b) = ".((defined $exp)?$exp:"<undef>") );
}

###### subint
{ my(@got,@exp);
  for my $a (-3 .. 3) {
    for my $b (-3 .. 3) {
      push @got, subint($a,$b);
      push @exp, $a-$b;
    }
  }
  is_deeply( \@got, \@exp, "subint( -3 .. 3, -3 .. 3)" );
}
foreach my $r (@subints) {
  my($a, $b, $exp) = @$r;
  is( subint($a,$b), $exp, "subint($a,$b) = ".((defined $exp)?$exp:"<undef>") );
}

###### add1int / sub1int
{
  my @N = (-17 .. 17,
           "4294967295", "4294967296", "4294967297",
           "9223372036854775807", "9223372036854775808", "9223372036854775809",
           "18446744073709551615", "18446744073709551616", "18446744073709551617",
           "158456325028528675187087900671");
  is_deeply([map { add1int($_) } @N], [map { addint($_,1) } @N], "add1int");
  is_deeply([map { sub1int($_) } @N], [map { subint($_,1) } @N], "sub1int");
}

###### divint
ok(!eval { divint(0,0); }, "divint(0,0)");
ok(!eval { divint(1,0); }, "divint(1,0)");

# For negative inputs, the div and mod operations might be different than Perl's builtins.
# It matches Math::BigInt bdiv / bmod (post 1.997 Sep 2015).

my @qpos1024 = map { int(1024/$_) } 1 .. 1025;
my @qneg1024 = map { my $d=-1024/$_; my $i = int($d);  ($d==$i) ? $i : $i-1; } 1 .. 1025;

my @rpos1024 = map {  1024 - $_ * $qpos1024[$_-1] } 1 .. 1025;
my @rneg1024 = map { -1024 - $_ * $qneg1024[$_-1] } 1 .. 1025;

is_deeply( [map { divint(1024,$_) } 1..1025], \@qpos1024, "divint(1024,x) for 1 .. 1025" );
is_deeply( [map { divint(-1024,$_) } 1..1025], \@qneg1024, "divint(-1024,x) for 1 .. 1025" );

###### modint
ok(!eval { modint(0,0); }, "modint(0,0)");
ok(!eval { modint(1,0); }, "modint(1,0)");

is_deeply( [map { modint(1024,$_) } 1..1025], \@rpos1024, "modint(1024,x) for 1 .. 1025" );
is_deeply( [map { modint(-1024,$_) } 1..1025], \@rneg1024, "modint(-1024,x) for 1 .. 1025" );

is(modint("-1117091728166568014",59), 4, "modint(-1117091728166568014,59) = 4");

###### divrem
ok(!eval { divrem(0,0); }, "divrem(0,0)");
ok(!eval { divrem(1,0); }, "divrem(1,0)");

###### fdivrem
ok(!eval { fdivrem(0,0); }, "fdivrem(0,0)");
ok(!eval { fdivrem(1,0); }, "fdivrem(1,0)");

###### cdivrem
ok(!eval { cdivrem(0,0); }, "cdivrem(0,0)");
ok(!eval { cdivrem(1,0); }, "cdivrem(1,0)");

###### tdivrem
ok(!eval { tdivrem(0,0); }, "tdivrem(0,0)");
ok(!eval { tdivrem(1,0); }, "tdivrem(1,0)");


###### large values through divint, cdivint, modint,
######                      divrem, tdivrem, fdivrem, cdivrem
for my $s (@quotients) {
  my($signs, $n, $m, $qt, $qf, $qc, $qe) = @$s;
  my($bn,$bm) = map { Math::BigInt->new($_) } ($n,$m);
  my($rt, $rf, $rc, $re) = map { $bn - $bm * $_ } ($qt, $qf, $qc, $qe);
  is( divint($n, $m), $qf, "large divint  $signs" );
  is( modint($n, $m), $rf, "large modint  $signs" );
  is( cdivint($n, $m), $qc, "large cdivint  $signs" );
  is_deeply( [divrem($n, $m)], [$qe, $re], "large divrem  $signs" );
  is_deeply( [tdivrem($n, $m)], [$qt, $rt], "large tdivrem $signs" );
  is_deeply( [fdivrem($n, $m)], [$qf, $rf], "large fdivrem $signs" );
  is_deeply( [cdivrem($n, $m)], [$qc, $rc], "large cdivrem $signs" );
}

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
is(divint("1895315831", -1), "-1895315831", "Divide 31-bit input by -1");
is(divint("3483637757", -1), "-3483637757", "Divide 32-bit input by -1");
is(cdivint("3483637757", -1), "-3483637757", "Divide 32-bit input by -1 (ceiling)");
is(divint("6127303089832103323", -1), "-6127303089832103323", "Divide 63-bit input by -1");
is(divint("13026328650942325963", -1), "-13026328650942325963", "Divide 64-bit input by -1");
is(divint("14123555781055773270", 2), "7061777890527886635", "Divide 64-bit input by 2");
is(divint("12844039487317506779", "12844039487317506779"), 1, "Divide 64-bit input by itself");
is(divint(3, "12844039487317506779"), 0, "Divide small int by 64-bit input");
# Note this is floor division:
is(divint(-3, "12844039487317506779"), -1, "Divide negative small int by 64-bit input");
# Now ceiling
is(cdivint(3, "12844039487317506779"), 1, "Divide (ceil) small int by 64-bit input");
is(cdivint(-3, "12844039487317506779"), 0, "Divide (ceil) negative small int by 64-bit input");

###### cdivrem special test
is( join(" ",cdivrem(3, "12844039487317506779")), "1 -12844039487317506776", "cdivrem with small quotient and 64-bit denominator shouldn't overflow IV" );
{
 my ($x248, $x263, $x264) = (powint(2,48), powint(2,63), powint(2,64));
 is_deeply([cdivint($x263-1,$x248),cdivint($x263,$x248),cdivint($x263+1,$x248)],
           [32768, 32768, 32769],
           "cdivint (2^63 +/- eps) / 2^48");
 is_deeply([cdivint(subint($x264,1),$x248),cdivint($x264,$x248),cdivint(addint($x264,1),$x248)],
           [65536, 65536, 65537],
           "cdivint (2^64 +/- eps) / 2^48");
}

###### lshiftint
is_deeply([map { lshiftint($_) } 0..50], [map { $_ << 1 } 0..50], "lshiftint(0..50)");
is_deeply([map { rshiftint($_) } 0..50], [map { $_ >> 1 } 0..50], "rshiftint(0..50)");
is_deeply([map { rashiftint($_) } 0..50], [map { $_ >> 1 } 0..50], "rashiftint(0..50)");
is_deeply([map { lshiftint($_,5) } -65 .. 65], [map { $_ * 32 } -65 .. 65], "lshiftint(-65 .. 65, 5)");

for my $d (@negshifts) {
  my($n, $k, $rs, $ras) = @$d;
  my $ls = mulint($n, powint(2,$k));
  is( lshiftint($n,$k), $ls, "lshiftint($n,$k) = $ls" );
  is( rshiftint($n,$k), $rs, "rshiftint($n,$k) = $rs" );
  is( rashiftint($n,$k), $ras, "rashiftint($n,$k) = $ras" );
}

###### absint
{ my(@got,@exp);
  for my $n (-100 .. 100) {
    push @got, absint($n);
    push @exp, abs($n);
  }
  is_deeply( \@got, \@exp, "absint( -100 .. 100)" );
  is( absint("0"), 0, "absint(0) = 0" );
  is( absint("-0"), 0, "absint(-0) = 0" );
}
{ my(@got,@exp);
  for my $n (@anvals) {
    my $av = $n;  $av =~ s/^-//;
    push @got, absint($n);
    push @exp, $av;
  }
  is_deeply( \@got, \@exp, "absint on large values" );
}

###### negint
{ my(@got,@exp);
  for my $n (-100 .. 100) {
    push @got, negint($n);
    push @exp, ($n == 0) ? 0 : -$n;
  }
  is_deeply( \@got, \@exp, "negint( -100 .. 100)" );
  is( negint("0"), 0, "negint(0) = 0" );
  is( negint("-0"), 0, "negint(-0) = 0" );
}
{ my(@got,@exp);
  for my $n (@anvals) {
    my $nv = $n;
    if ($nv =~ /^-/) {
      $nv =~ s/^-//;
    } else {
      $nv = "-$nv";
    }
    push @got, negint($n);
    push @exp, $nv;
  }
  is_deeply( \@got, \@exp, "negint on large values" );
}
{
  my @pos = (qw/1073741823 1073741824 1073741825 2147483647 2147483648 2147483649 4294967295 4294967296 4294967297 8589934591 8589934592 8589934593 4611686018427387903 4611686018427387904 4611686018427387905 9223372036854775807 9223372036854775808 9223372036854775809 18446744073709551615 18446744073709551616 18446744073709551617 36893488147419103231 36893488147419103232 36893488147419103233 170141183460469231731687303715884105727 170141183460469231731687303715884105728 170141183460469231731687303715884105729 340282366920938463463374607431768211455 340282366920938463463374607431768211456 340282366920938463463374607431768211457/);
  my @neg = map { '-' . $_ } @pos;
  is_deeply( [map { negint($_) } @pos], \@neg, "negint on positive powers of 2 transitions");
  is_deeply( [map { negint($_) } @neg], \@pos, "negint on negative powers of 2 transitions");
}

###### cmpint

is(cmpint(1,2),-1,"1 < 2");
is(cmpint(2,1), 1,"2 > 1");
is(cmpint(2,2), 0,"2 == 2");
is(cmpint("18446744073709553664","18446744073709551615"),1,"2^64+2048 > 2^64-1");
is(cmpint("18446744073709551664","18446744073709551615"),1,"2^64+1048 > 2^64-1");

###### signint
is(signint(-13), -1, "signint(-13) = -1");
is(signint(0), 0, "signint(0) = 0");
is(signint(13), 1, "signint(13) = 1");
