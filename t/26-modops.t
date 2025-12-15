#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/negmod invmod addmod submod mulmod muladdmod mulsubmod divmod powmod/;
use Math::BigInt try=>"GMP,GMPz,Pari";

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
$use64 = 0 if $use64 && 18446744073709550592 == ~0;
my $usexs = Math::Prime::Util::prime_get_config->{'xs'};

my @invmods = (
 [ 0, 0, undef],
 [ 1, 0, undef],
 [ 0, 1, 0],
 [ 0, 2, undef],
 [ 1, 1, 0],
 [ 45, 59, 21],
 [  42,  2017, 1969],
 [  42, -2017, 1969],
 [ -42,  2017, 48],
 [ -42, -2017, 48],
 [ 14, 28474, undef],
 [ 13, "9223372036854775808", "5675921253449092805" ],
 [ 14, "18446744073709551615", "17129119497016012214" ],
);
plan tests => 0
            + 10                     # negmod
            + 5 + scalar(@invmods)   # invmod
            + 5*2
            + 1                      # addmod
            + 2                      # submod / addmod
            + 2                      # mulmod
            + 2 + 1                  # divmod
            + 2                      # powmod
            + 6                      # large negative args
            + 1                      # muladdmod
            + 1                      # mulsubmod
            + 4                      # muladdmod and mulsubmod large inputs
            + 1                      # large negative modulus
            + 0;

###### negmod

# For n != 0, negmod(a,n) = modint(-a,n).
# For all inputs, negmod(a,n) = $n ? submod(0, modint(a,n), n) : undef;
is(negmod(0,0), undef, "negmod(0,0) = undef");
is(negmod(1,0), undef, "negmod(1,0) = undef");
is(negmod(0,1), 0, "negmod(0,1) = 0");
is(negmod(100,1), 0, "negmod(100,1) = 0");
is(negmod( 100, 123), 23, "negmod(100, 123) = 23");
is(negmod( 100,-123), 23, "negmod(100,-123) = 23");
is(negmod(-100, 123), 100, "negmod(-100, 123) = 100");
is(negmod( 10000, 123), 86, "negmod(10000, 123) = 86");
is(negmod( 10000,-123), 86, "negmod(10000,-123) = 86");
is(negmod(-10000, 123), 37, "negmod(-10000, 123) = 37");

###### invmod
ok(!eval { invmod(undef,11); }, "invmod(undef,11)");
ok(!eval { invmod(11,undef); }, "invmod(11,undef)");
ok(!eval { invmod('nan',11); }, "invmod('nan',11)");

foreach my $r (@invmods) {
  my($a, $n, $exp) = @$r;
  my $got = invmod($a,$n);
  if (!defined $exp) {
    is($got, $exp, "invmod($a,$n) = <undef>");
  } elsif (!defined $got) {
    is($got, $exp, "invmod($a,$n) = $exp");
  } else {
    is("$got", $exp, "invmod($a,$n) = $exp");
  }
}
# Pari, Mathematica, SAGE, Math::BigInt  all return 0 for this case.
is( invmod(0,1), 0, "invmod(0,1) = 0");
is( invmod(0,-1), 0, "invmod(0,-1) = 0");
# my $res = invmod(0,1);   $res = "<undef>" if !defined $res;
# ok($res eq '0' || $res eq '<undef>', "invmod(0,1) = $res");


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
  is( divmod($a,$b,$m),  160625959, "divmod with large negative arg" );
  is( powmod($a,$b,$m), 1550454861, "powmod with large negative arg" );
  is( powmod($b,$a,$m),   16491583, "powmod with large negative arg" );
}


my @ic = map { nrand() } 0 .. $num;

###### muladdmod
@exp = (); @res = ();
for (0 .. $num) {
  push @exp, Math::BigInt->new($i1[$_])->bmul(-$i2t[$_])->badd($ic[$_])->bmod($i3[$_]);
  push @res, muladdmod($i1[$_], -$i2t[$_], $ic[$_], $i3[$_]);
}
is_deeply( \@res, \@exp, "muladdmod on ".($num+1)." random inputs" );

###### mulsubmod
@exp = (); @res = ();
for (0 .. $num) {
  push @exp, Math::BigInt->new($i1[$_])->bmul(-$i2t[$_])->bsub($ic[$_])->bmod($i3[$_]);
  push @res, mulsubmod($i1[$_], -$i2t[$_], $ic[$_], $i3[$_]);
}
is_deeply( \@res, \@exp, "mulsubmod on ".($num+1)." random inputs" );

# Arbitrary non-tiny values
is("".muladdmod("293482938498234","982498230923490234234","982349823092355","87777777777757"), "20728855000562", "muladdmod with medium size inputs");
is("".mulsubmod("293482938498234","982498230923490234234","982349823092355","87777777777757"), "74918097704263", "mulsubmod with medium size inputs");
# 128-bit inputs mod a 126-bit prime
is("".muladdmod("175109911729618543589989257539043768012","21887412602962542281538131483385626868","263466159656861646486075450888763957942","83494980727347746728226137271418851789"), "28511529282241296497665677199750506129", "muladdmod with 128-bit inputs mod a 126-bit prime");
is("".mulsubmod("171821502870939196679518625154011220409","182569474286058024586486590841354369890","329619958784558749434516006469339236320","49684876044205406960769234394385141897"), "43334010019275970236275282850802256285", "mulsubmod with 128-bit inputs mod a 126-bit prime");

subtest 'big raw negative mod ', sub {
  is("".addmod("18446744073709551615","18446744073709551615","-19446744073709551616"),"17446744073709551614");
  is("".submod("17446744073709551614",0,"-19446744073709551616"),"17446744073709551614");
  is("".mulmod("18446744073709551615",2,"-19446744073709551616"),"17446744073709551614");
};


sub nrand {
  my $r = int(rand(4294967296));
  $r = ($r << 32) + int(rand(4294967296)) if $use64;
  $r;
}
