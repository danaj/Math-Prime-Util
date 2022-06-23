#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/invmod sqrtmod allsqrtmod rootmod allrootmod addmod submod mulmod divmod powmod/;
use Math::BigInt try=>"GMP,Pari";

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
$use64 = 0 if $use64 && 18446744073709550592 == ~0;
my $usexs = Math::Prime::Util::prime_get_config->{'xs'};

my @invmods = (
 [ 0, 0, undef],
 [ 1, 0, undef],
 [ 1, 1, 0],
 [ 45, 59, 21],
 [  42,  2017, 1969],
 [  42, -2017, 1969],
 [ -42,  2017, 48],
 [ -42, -2017, 48],
 [ 14, 28474, undef],
);
if ($use64) {
 push @invmods, [ 13, 9223372036854775808, 5675921253449092805 ];
 push @invmods, [ 14, 18446744073709551615, 17129119497016012214 ];
}
my @sqrtmods = (
 [ 0, 0, undef ],
 [ 1, 0, undef ],
 [ 0, 1, 0 ],
 [ 1, 1, 0 ],

 # prime moduli
 [ -1, 17, [4,13], 1 ],  # github #51
 [ 58, 101, [19,82], 1 ],
 [ 111, 113, [26,87], 1 ],
 [ 160, 461, undef ],
 [ 37, 999221, [9946,989275], 1 ],
 [ 30, 1000969, [89676,911293], 1 ],
 [ "9223372036854775808", "5675921253449092823", ["22172359690642254","5653748893758450569"], 1 ],
 [ "18446744073709551625", "340282366920938463463374607431768211507", ["57825146747270203522128844001742059051","282457220173668259941245763430026152456"], 1 ],
 [ 2, 72388801, [20312446,52076355], 1 ],

 # composite moduli
 [ 30, 74, [20,54] ],
 [ 56, 1018, [458,560] ],
 [ 42, 979986, [356034,623952] ],
 [ 5, 301, undef ],
 [ 5, 302, [55,247] ],
 [ 5, 404, [45,157,247,359] ],
 [ 5, 400, undef ],
 [ 9, 400, [3,53,147,197,203,253,347,397] ],
 [ 15, 402, [45,357] ],
 [ 1242, 1849, [851, 998] ],  # prime power
 [ 0, 4, [0,2] ],
 [ 1, 4, [1,3] ],
 [ 4, 8, [2,6] ],        # github #52
 [ 4, 16, [2,6,10,14] ],
 [ 0, 9, [0,3,6] ],
 [ 3, 9, undef ],
 [ 0, 27, [0,9,18] ],
 [ 9, 27, [3,6,12,15,21,24] ],
 [ 0, 36, [0,6,12,18,24,30] ],
 [ 4, 36, [2,16,20,34] ],
 [ 13556, 26076, undef ],
 [ 15347, 38565, undef ],
 [ 588,2912, undef ],
 [ 24684, 69944, [2138,17126,17846,32834,37110,52098,52818,67806] ],
);

# Check:
#   p {prime, prime power, square-free composite, non-SF composite}
#   k {prime, prime power, square-free composite, non-SF composite}
my @rootmods = (
 # prime moduli
 [14,-3, 101, 17],
 [13, 6, 107, [24,83]],
 [13,-6, 107, [49,58]],
 [64, 6, 101, [2, 99]],
 [9, -2, 101, [34, 67]],
 [ 2, 3, 3, 2],
 [ 2, 3, 7, undef],
 [17, 29, 19, 6],
 [ 5, 3, 13, [7,8,11]],
 [53, 3, 151, [15,27,109]],
 [3,3,73, [25,54,67]],
 [7,3,73, [13,29,31]],
 [49,3,73, [12,23,38]],
 [44082,4,100003, [2003,98000]],
 [90594,6,100019, [37071,62948]],
 [6,5,31, [11,13,21,22,26]],
 [0,2,2, 0],
 [2,4,5, undef],
 [51,12,10009,[64,1203,3183,3247,3999,4807,5202,6010,6762,6826,8806,9945]],
 [15,3,"1000000000000000000117",[qw/72574612502199260377 361680004182786118804 565745383315014620936/]],
 [1,0,13, [0,1,2,3,4,5,6,7,8,9,10,11,12]],
 [2,0,13, undef],
 [0,5,0,  undef],
 [0,-1,3, undef],

 # composite moduli.
 # Pari will usually give a *wrong* answer for these if using Mod(a,p).
 # The right way with Pari is to use p-adic.
 [ 4, 2, 10, [2,8]],
 [ 4, 2, 18, [2,16]],
 [ 2, 3, 21, undef],  # Pari says 2
 [ 8, 3, 27, [2,11,20]],      # Pari says 26
 [22, 3, 1505, [148,578,673,793,813,1103,1243,1318,1458] ],  # Pari says 1408
 [58787, 3, 100035, [3773,8633,10793,13763,19163,24293,26183,26588,31313,37118,41978,44138,47108,52508,57638,59528,59933,64658,70463,75323,77483,80453,85853,90983,92873,93278,98003]],
 [3748, 2, 4992, [154,262,314,518,730,934,986,1094,1402,1510,1562,1766,1978,2182,2234,2342,2650,2758,2810,3014,3226,3430,3482,3590,3898,4006,4058,4262,4474,4678,4730,4838]],
 [68,2,2048, [46,466,558,978,1070,1490,1582,2002]],
 [96,5,128, [6,14,22,30,38,46,54,62,70,78,86,94,102,110,118,126]],
 [2912,5,4992, [182,494,806,1118,1430,1742,2054,2366,2678,2990,3302,3614,3926,4238,4550,4862]],
 [ 2, 3, 4, undef],
 [ 3, 2, 4, undef],
 [ 3, 4, 19, undef],
 [ 1, 4, 20, [1,3,7,9,11,13,17,19]],
 [ 9, 2, 24, [3,9,15,21]],
 [ 6, 6, 35, undef],
 [ 36,2, 40, [6,14,26,34]],
 [ 16, 12, 48, [2,4,8,10,14,16,20,22,26,28,32,34,38,40,44,46]],
 [ 13, 6, 112, undef],
 [ 52, 6, 117, undef],
 [ 48, 3, 128, undef],
 [ 382, 3, 1000, undef],
 [ 10, 3, 81, [13,40,67]],
 [ 26, 5, 625, [81,206,331,456,581]],
 [ 51, 5, 625, [61,186,311,436,561]],
 ["9833625071",3,"10000000071", [qw/3333332807 6666666164 9999999521/]],
 #[2131968,5,10000000000, [...]],   # Far too many
 [198,-1,519, undef],
);

plan tests => 0
            + 5 + scalar(@invmods)
            + scalar(@sqrtmods)*2    # sqrtmod / allsqrtmod
            + 5*2
            + 1                      # addmod
            + 2                      # submod / addmod
            + 2                      # mulmod
            + 2 + 1                  # divmod
            + 2                      # powmod
            + 6                      # large negative args
            + 5                      # rootmod
            + scalar(@rootmods)*2    # allrootmod
            + 1                      # more rootmod
            + 0;

###### invmod
ok(!eval { invmod(undef,11); }, "invmod(undef,11)");
ok(!eval { invmod(11,undef); }, "invmod(11,undef)");
ok(!eval { invmod('nan',11); }, "invmod('nan',11)");

foreach my $r (@invmods) {
  my($a, $n, $exp) = @$r;
  is( invmod($a,$n), $exp, "invmod($a,$n) = ".((defined $exp)?$exp:"<undef>") );
}
# Pari, Mathematica, SAGE, Math::BigInt  all return 0 for this case.
is( invmod(0,1), 0, "invmod(0,1) = 0");
is( invmod(0,-1), 0, "invmod(0,-1) = 0");
# my $res = invmod(0,1);   $res = "<undef>" if !defined $res;
# ok($res eq '0' || $res eq '<undef>', "invmod(0,1) = $res");

###### sqrtmod
foreach my $r (@sqrtmods) {
  my($a, $n, $exp, $prime) = @$r;
  if (!defined $exp) {
    is( sqrtmod($a,$n), $exp, "sqrtmod($a,$n) = <undef>");
    is_deeply( [allsqrtmod($a,$n)], [], "allsqrtmod($a,$n) = ()");
  } elsif (!ref($exp)) {
    is( sqrtmod($a,$n), $exp, "sqrtmod($a,$n) = $exp");
    is_deeply( [allsqrtmod($a,$n)], [$exp], "allsqrtmod($a,$n) = ($exp)");
  } else {
    my $val = sqrtmod($a,$n);
    if ($prime) {
        # sqrtmod() must return least root for prime modulus
        is($val, $exp->[0], "sqrtmod($a,$n) = $exp->[0]");
    } else {
        ok( is_one_of($val, @$exp), "sqrtmod($a,$n) = $val, roots [@$exp]" );
    }
    is_deeply([allsqrtmod($a,$n)], $exp, "allsqrtmod($a,$n) = (@$exp)");
  }
}


my $num = 99;
$num = 29 if Math::BigInt->config()->{lib} !~ /(GMP|Pari)/;

my @i1 = map { nrand() } 0 .. $num;
my @i2 = map { nrand() } 0 .. $num;
my @i2t= map { $i2[$_] >> 1 } 0 .. $num;
my @i3 = map { nrand() } 0 .. $num;
my(@exp,@res);


###### add/mul/div/pow with small arguments
@exp = map { undef } 0..27;
is_deeply(\@exp, [map { addmod($_ & 3, ($_>>2)-3, 0) } 0..27], "addmod(..,0)");
is_deeply(\@exp, [map { submod($_ & 3, ($_>>2)-3, 0) } 0..27], "submod(..,0)");
is_deeply(\@exp, [map { mulmod($_ & 3, ($_>>2)-3, 0) } 0..27], "mulmod(..,0)");
is_deeply(\@exp, [map { divmod($_ & 3, ($_>>2)-3, 0) } 0..27], "divmod(..,0)");
is_deeply(\@exp, [map { powmod($_ & 3, ($_>>2)-3, 0) } 0..27], "powmod(..,0)");

@exp = map { 0 } 0..27;
is_deeply(\@exp, [map { addmod($_ & 3, ($_>>2)-3, 1) } 0..27], "addmod(..,1)");
is_deeply(\@exp, [map { submod($_ & 3, ($_>>2)-3, 1) } 0..27], "submod(..,1)");
is_deeply(\@exp, [map { mulmod($_ & 3, ($_>>2)-3, 1) } 0..27], "mulmod(..,1)");
is_deeply(\@exp, [map { divmod($_ & 3, ($_>>2)-3, 1) } 0..27], "divmod(..,1)");
is_deeply(\@exp, [map { powmod($_ & 3, ($_>>2)-3, 1) } 0..27], "powmod(..,1)");


###### addmod
@exp = (); @res = ();
for (0 .. $num) {
  push @exp, Math::BigInt->new("$i1[$_]")->badd("$i2[$_]")->bmod("$i3[$_]");
  push @res, addmod($i1[$_], $i2[$_], $i3[$_]);
}
is_deeply( \@res, \@exp, "addmod on ".($num+1)." random inputs" );

###### submod
@exp = (); @res = ();
for (0 .. $num) {
  push @exp, Math::BigInt->new("$i1[$_]")->bsub("$i2t[$_]")->bmod("$i3[$_]");
  push @res, submod($i1[$_], $i2t[$_], $i3[$_]);
}
is_deeply( \@res, \@exp, "submod on ".($num+1)." random inputs" );
##### addmod with negative
@res = ();
for (0 .. $num) {
  push @res, addmod($i1[$_], -$i2t[$_], $i3[$_]);
}
is_deeply( \@res, \@exp, "addmod with negative second input on ".($num+1)." random inputs" );

###### mulmod
@exp = (); @res = ();
for (0 .. $num) {
  push @exp, Math::BigInt->new("$i1[$_]")->bmul("$i2[$_]")->bmod("$i3[$_]");
  push @res, mulmod($i1[$_], $i2[$_], $i3[$_]);
}
is_deeply( \@res, \@exp, "mulmod on ".($num+1)." random inputs" );

###### mulmod (neg)
@exp = (); @res = ();
for (0 .. $num) {
  push @exp, Math::BigInt->new("$i1[$_]")->bmul("-$i2t[$_]")->bmod("$i3[$_]");
  push @res, mulmod($i1[$_], -$i2t[$_], $i3[$_]);
}
is_deeply( \@res, \@exp, "mulmod with negative second input on ".($num+1)." random inputs" );

###### divmod
is(divmod(0,14,53), 0, "divmod(0,14,53) = mulmod(0,invmod(14,53),53) = mulmod(0,19,53) = 0");

@exp = (); @res = ();
for (0 .. $num) {
  push @exp, Math::BigInt->new("$i2[$_]")->bmodinv("$i3[$_]")->bmul("$i1[$_]")->bmod("$i3[$_]");
  push @res, divmod($i1[$_], $i2[$_], $i3[$_]);
}
@exp = map { $_->is_nan() ? undef : $_ } @exp;
is_deeply( \@res, \@exp, "divmod on ".($num+1)." random inputs" );

###### divmod (neg)
@exp = (); @res = ();
# Old Math::BigInt will die with FP exception.  Work around.
#for (0 .. $num) {
#  push @exp, Math::BigInt->new("-$i2t[$_]")->bmodinv("$i3[$_]")->bmul("$i1[$_]")->bmod("$i3[$_]");
#  push @res, divmod($i1[$_], -$i2t[$_], $i3[$_]);
#}
#@exp = map { $_->is_nan() ? undef : $_ } @exp;
for (0 .. $num) {
  my $r = divmod($i1[$_], -$i2t[$_], $i3[$_]);
  push @res, $r;
  if (defined $r) {
    push @exp, Math::BigInt->new("-$i2t[$_]")->bmodinv("$i3[$_]")->bmul("$i1[$_]")->bmod("$i3[$_]");
  } else {
    push @exp, undef;
  }
}
is_deeply( \@res, \@exp, "divmod with negative second input on ".($num+1)." random inputs" );

###### powmod
@exp = (); @res = ();
for (0 .. $num) {
  push @exp, Math::BigInt->new("$i1[$_]")->bmodpow("$i2[$_]","$i3[$_]");
  push @res, powmod($i1[$_], $i2[$_], $i3[$_]);
}
is_deeply( \@res, \@exp, "powmod on ".($num+1)." random inputs" );

###### powmod (neg)
@exp = (); @res = ();
for (0 .. $num) {
  push @exp, Math::BigInt->new("$i1[$_]")->bmodpow("-$i2t[$_]","$i3[$_]");
  push @res, powmod($i1[$_], -$i2t[$_], $i3[$_]);
}
@exp = map { $_->is_nan() ? undef : $_ } @exp;
is_deeply( \@res, \@exp, "powmod with negative exponent on ".($num+1)." random inputs" );

###### large negative args (github issue 43)
{
  my($a, $b, $m) = (1363362182, "-26315271553053477373", 2000000011);
  is( addmod($a,$b,$m), 1043877553, "addmod with large negative arg" );
  is( submod($a,$b,$m), 1682846811, "submod with large negative arg" );
  is( mulmod($a,$b,$m), 1486752452, "mulmod with large negative arg" );
  is( divmod($a,$b,$m),  160625959, "mulmod with large negative arg" );
  is( powmod($a,$b,$m), 1550454861, "powmod with large negative arg" );
  is( powmod($b,$a,$m),   16491583, "powmod with large negative arg" );
}


###### rootmod
{
  my(@out0,@out1);
  for my $a (0..3) {
    for my $k (0..3) {
      push @out0, "$a,$k,0"  if defined rootmod($a,$k,0);
      push @out1, "$a,$k,1"  unless iseq(0,rootmod($a,$k,1));
    }
  }
  is(join('  ',@out0),'',"rootmod(a,k,0) should be undef");
  is(join('  ',@out1),'',"rootmod(a,k,1) should be 0");
}
{
  my(@out0,@out1,@out2);
  for my $a (0..19) {
    my $a17 = $a % 17;
    push @out0, "$a,0,17" if ($a17 == 1 && !iseq(1,rootmod($a,0,17)))
                          || ($a17 != 1 && defined rootmod($a,0,17));
    push @out1, "$a,1,17"  unless iseq($a17,rootmod($a,1,17));
    push @out2, "$a,2,17"  unless iseq(sqrtmod($a,17),  rootmod($a,2,17));
  }
  is(join('  ',@out0),'',"rootmod(a,0,17) should be 1 or undef");
  is(join('  ',@out1),'',"rootmod(a,1,17) should be a mod 17");
  is(join('  ',@out2),'',"rootmod(a,2,17) should be sqrtmod(a,17)");
}

foreach my $r (@rootmods) {
  my($a, $k, $n, $exp) = @$r;
  if (!defined $exp) {
    is( rootmod($a,$k,$n), $exp, "rootmod($a,$k,$n) = <undef>");
    is_deeply( [allrootmod($a,$k,$n)], [], "allrootmod($a,$k,$n) = ()");
  } elsif (!ref($exp)) {
    is( rootmod($a,$k,$n), $exp, "rootmod($a,$k,$n) = $exp");
    is_deeply( [allrootmod($a,$k,$n)], [$exp], "allrootmod($a,$k,$n) = ($exp)");
  } else {
    my $val = rootmod($a,$k,$n);
    ok( is_one_of($val, @$exp), "rootmod($a,$k,$n) = $val, roots [@$exp]" );
    is_deeply([allrootmod($a,$k,$n)], $exp, "allrootmod($a,$k,$n) = (@$exp)");
  }
}

is(powmod(rootmod(12,41,1147),41,1147), 12, "41st root of 12 mod 1147 is correct");




sub nrand {
  my $r = int(rand(4294967296));
  $r = ($r << 32) + int(rand(4294967296)) if $use64;
  $r;
}

sub is_one_of {
  my($n, @list) = @_;
  if (defined $n) {
    for (@list) {
      return 1 if defined $_ && $n eq $_;
    }
  } else {
    for (@list) {
      return 1 if !defined $_;
    }
  }
  0;
}

sub iseq {    # numerical comparison allowing undef = undef
  my($x,$y) = @_;
  return 1 if !defined $x && !defined $y;
  return 0 if !defined $x || !defined $y;
  $x == $y;
}
