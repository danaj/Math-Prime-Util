#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/ powint mulint addint divint modint divrem tdivrem /;
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
  ["13282407956253574712","14991082624209354397","199117675120653046511338473800925208664"],
);
my @addints = (
  ["1178630961471601951655862","827639478068904540012","1179458600949670856195874"],
  ["-2555488174170453670799","1726145541361106236340","-829342632809347434459"],
);
my @quotients = (  # trunc, floor, euclidian
  ["+ +",  "39458349850349850394853049583049",  "85889",  "459410982202026457344398579",  "459410982202026457344398579",  "459410982202026457344398579"],
  ["+ -",  "39458349850349850394853049583049", "-85889", "-459410982202026457344398579", "-459410982202026457344398580", "-459410982202026457344398579"],
  ["- +", "-39458349850349850394853049583049",  "85889", "-459410982202026457344398579", "-459410982202026457344398580", "-459410982202026457344398580"],
  ["- -", "-39458349850349850394853049583049", "-85889",  "459410982202026457344398579",  "459410982202026457344398579",  "459410982202026457344398580"],
);

plan tests => 0
            + 7*4 + scalar(@powints) + 3     # powint
            + 1 + scalar(@mulints)           # mulint
            + 1 + scalar(@addints)           # addint
            + 2 + 2                          # divint
            + 2 + 2                          # modint
            + 2                              # divrem
            + 2                              # tdivrem
            + 4 * scalar(@quotients)         # signed bigint division
            + 6                              # table 1.3 from Leijen 2001
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

###### mulint
{ my(@got,@exp);
  for my $a (-3 .. 3) {
    for my $b (-3 .. 3) {
      push @got, mulint($a,$b);
      push @exp, $a*$b;
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

###### divint
ok(!eval { divint(0,0); }, "divint(1,0)");
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
ok(!eval { modint(0,0); }, "modint(1,0)");
ok(!eval { modint(1,0); }, "modint(1,0)");

is_deeply( [map { modint(1024,$_) } 1..1025], \@rpos1024, "modint(1024,x) for 1 .. 1025" );
is_deeply( [map { modint(-1024,$_) } 1..1025], \@rneg1024, "modint(-1024,x) for 1 .. 1025" );

###### divrem
ok(!eval { divrem(0,0); }, "divrem(1,0)");
ok(!eval { divrem(1,0); }, "divrem(1,0)");

###### tdivrem
ok(!eval { tdivrem(0,0); }, "tdivrem(1,0)");
ok(!eval { tdivrem(1,0); }, "tdivrem(1,0)");


###### large values through divint, modint, divrem, tdivrem
for my $s (@quotients) {
  my($signs, $n, $m, $qt, $qf, $qe) = @$s;
  my($bn,$bm) = map { Math::BigInt->new($_) } ($n,$m);
  my($rt, $rf, $re) = map { $bn - $bm * $_ } ($qt, $qf, $qe);
  is( divint($n, $m), $qf, "large divint  $signs" );
  is( modint($n, $m), $rf, "large modint  $signs" );
  is_deeply( [divrem($n, $m)], [$qe, $re], "large divrem  $signs" );
  is_deeply( [tdivrem($n, $m)], [$qt, $rt], "large tdivrem $signs" );
}

###### small values through divint, modint, divrem, tdivrem
is( join(" ", tdivrem(8,3), tdivrem(8,-3), tdivrem(-8,3), tdivrem(-8,-3)),
    "2 2 -2 2 -2 -2 2 -2",
    "tdivrem with +/- 8,3" );
is( join(" ", divrem(8,3), divrem(8,-3), divrem(-8,3), divrem(-8,-3)),
    "2 2 -2 2 -3 1 3 1",
    "divrem with +/- 8,3" );
is( join(" ", divint(8,3), modint(8,3), divint(8,-3), modint(8,-3), divint(-8,3), modint(-8,3), divint(-8,-3), modint(-8,-3)),
    "2 2 -3 -1 -3 1 2 -2",
    "divint+modint with +/- 8,3" );

is( join(" ", tdivrem(1,2), tdivrem(1,-2), tdivrem(-1,2), tdivrem(-1,-2)),
    "0 1 0 1 0 -1 0 -1",
    "tdivrem with +/- 1,2" );
is( join(" ", divrem(1,2), divrem(1,-2), divrem(-1,2), divrem(-1,-2)),
    "0 1 0 1 -1 1 1 1",
    "divrem with +/- 1,2" );
is( join(" ", divint(1,2), modint(1,2), divint(1,-2), modint(1,-2), divint(-1,2), modint(-1,2), divint(-1,-2), modint(-1,-2)),
    "0 1 -1 -1 -1 1 0 -1",
    "divint+modint with +/- 1,2" );
