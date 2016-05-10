#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/invmod sqrtmod addmod mulmod divmod powmod/;
use Math::BigInt try=>"GMP,Pari";

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
$use64 = 0 if $use64 && 18446744073709550592 == ~0;

my @invmods = (
 [ 0, 0, undef],
 [ 1, 0, undef],
 [ 0, 1, undef],
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
 [ 0, 0, undef],
 [ 1, 0, undef],
 [ 0, 1, 0],
 [ 1, 1, 0],
 [ 58, 101, 19],
 [ 111, 113, 26],
 [ 37, 999221, 9946],
 [ 30, 1000969, 89676],
 [ "9223372036854775808", "5675921253449092823", "22172359690642254" ],
 [ "18446744073709551625", "340282366920938463463374607431768211507", "57825146747270203522128844001742059051" ],
 [ 30, 74, 20 ],
 [ 56, 1018, 458 ],
 [ 42, 979986, 356034 ],
);

plan tests => 0
            + 3 + scalar(@invmods)
            + scalar(@sqrtmods)
            + 4*2
            + 1                      # addmod
            + 1                      # submod
            + 2                      # mulmod
            + 2 + 1                  # divmod
            + 2                      # powmod
            + 0;

###### invmod
ok(!eval { invmod(undef,11); }, "invmod(undef,11)");
ok(!eval { invmod(11,undef); }, "invmod(11,undef)");
ok(!eval { invmod('nan',11); }, "invmod('nan',11)");

foreach my $r (@invmods) {
  my($a, $n, $exp) = @$r;
  is( invmod($a,$n), $exp, "invmod($a,$n) = ".((defined $exp)?$exp:"<undef>") );
}

###### sqrtmod
foreach my $r (@sqrtmods) {
  my($a, $n, $exp) = @$r;
  is( sqrtmod($a,$n), $exp, "sqrtmod($a,$n) = ".((defined $exp)?$exp:"<undef>") );
}


my $num = 99;
my @i1 = map { nrand() } 0 .. $num;
my @i2 = map { nrand() } 0 .. $num;
my @i2t= map { $i2[$_] >> 1 } 0 .. $num;
my @i3 = map { nrand() } 0 .. $num;
my(@exp,@res);


###### add/mul/div/pow with small arguments
@exp = map { 0 } 0..27;
is_deeply(\@exp, [map { addmod($_ & 3, ($_>>2)-3, 0) } 0..27], "addmod(..,0)");
is_deeply(\@exp, [map { mulmod($_ & 3, ($_>>2)-3, 0) } 0..27], "mulmod(..,0)");
is_deeply(\@exp, [map { divmod($_ & 3, ($_>>2)-3, 0) } 0..27], "divmod(..,0)");
is_deeply(\@exp, [map { powmod($_ & 3, ($_>>2)-3, 0) } 0..27], "powmod(..,0)");

is_deeply(\@exp, [map { addmod($_ & 3, ($_>>2)-3, 1) } 0..27], "addmod(..,1)");
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



sub nrand {
  my $r = int(rand(4294967296));
  $r = ($r << 32) + int(rand(4294967296)) if $use64;
  $r;
}
