#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/ znprimroot is_primitive_root qnr is_qr/;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
my $usegmp= Math::Prime::Util::prime_get_config->{'gmp'};
$use64 = 0 if $use64 && 18446744073709550592 == ~0;

my %primroots = (
   -11 => 2,
     1 => 0,
     2 => 1,
     3 => 2,
     4 => 3,
     5 => 2,
     6 => 5,
     7 => 3,
     8 => undef,
     9 => 2,
    10 => 3,          # 3 is the smallest root.  Pari gives the other root 7.
      1729 => undef,  # Pari goes into an infinite loop.
   5109721 =>  94,
  17551561 =>  97,
  90441961 => 113,
1407827621 =>   2,
1520874431 =>  17,
1685283601 => 164,
 100000001 => undef,  # Without an early exit, this will essentially hang.
);
if ($use64) {
  $primroots{2232881419280027} = 6;         # factor divide goes to FP
  $primroots{14123555781055773271} = 6;     # bmodpow hits RT 71548
  $primroots{89637484042681} = 335;         # smallest root is large
  $primroots{9223372036854775837} = 5;      # Pari #905
  $primroots{36002292036481} = 13;
  $primroots{72004584072962} = 13;
  $primroots{2067900233973681742} = 17;
  $primroots{8000468009126059319} = 13;
}
if ($usegmp || $extra) {  # in each case, p-1 is very easy to factor
  # p^2
  $primroots{"474264225821700214950222988868518911801235024731324721"} = 7;
  # 2p^2
  $primroots{"1580603145023079446166874838636458851122"} = 7;
  # p^3
  $primroots{"11154774760949852441478897023837868805975434161260919037124141673071282481903446814549"} = 2;
  # 2p^3
  $primroots{"44434394326141300867665315903406029736550298166159399085858"} = 23;
}

plan tests => 0
                + scalar(keys %primroots) + 1  # znprimroot
                + scalar(keys %primroots) + 4  # is_primitive_root
                + 7                            # qnr
                + 9                            # is_qr
                ;

###### znprimroot
while (my($n, $root) = each (%primroots)) {
  is( znprimroot($n), $root, "znprimroot($n) == " . ((defined $root) ? $root : "<undef>") );
}
is( znprimroot("-100000898"), 31, "znprimroot(\"-100000898\") == 31" );
# I don't think we should rely on this parsing correctly.
#is( znprimroot("+100000898"), 31, "znprimroot(\"+100000898\") == 31" );

###### is_primitive_root
while (my($n, $root) = each (%primroots)) {
  if (defined $root) {
    is( is_primitive_root(0+$root,$n), 1, "$root is a primitive root mod $n" );
  } else {
    is( is_primitive_root(2,$n), 0, "2 is not a primitive root mod $n" );
  }
}
is(is_primitive_root(2,0), undef, "is_primitive_root(2,0) => undef");
is(is_primitive_root(19,191), 1, "19 is a primitive root mod 191");
is(is_primitive_root(13,191), 0, "13 is not a primitive root mod 191");
is(is_primitive_root(35,982), 0, "35 is not a primitive root mod 982");

###### qnr
is(qnr(0), undef, "qnr(0) returns undef");
is_deeply([map{qnr($_)}1..15], [1,2,2,2,2,2,3,2,2,2,2,2,2,3,2], "qnr(1..15)");
is_deeply([map{qnr(2**$_)}1..16], [map{2}1..16], "qnr(2^k) = 2 for k>=1");
is(qnr(5711), 19, "The least quadratice non-residue of 5711 is 19");
is(qnr(366791), 43, "The least quadratice non-residue of 366791 is 43");
is(qnr(2737), 3, "qnr(7*17*23) = 2");
is(qnr(9257330), 2, "qnr(2*5*925733) = 2");

###### is_qr
is(is_qr(0,0), undef, "is_qr(x,0) returns undef");
is_deeply([map{is_qr($_,1)}0..5], [1,1,1,1,1,1], "is_qr(a,1) = 1");
is_deeply([map{is_qr($_,2)}0..5], [1,1,1,1,1,1], "is_qr(a,2) = 1");
is_deeply([map{is_qr($_,3)}0..10], [1,1,0,1,1,0,1,1,0,1,1], "is_qr(0..10,3)");
is_deeply([map{is_qr($_,4)}0..10], [1,1,0,0,1,1,0,0,1,1,0], "is_qr(0..10,4)");
is_deeply([map{is_qr($_,6)}0..10], [1,1,0,1,1,0,1,1,0,1,1], "is_qr(0..10,6)");
is_deeply([map{is_qr($_,9)}0..20], [1,1,0,0,1,0,0,1,0,1,1,0,0,1,0,0,1,0,1,1,0], "is_qr(0..20,9)");
is_deeply([map{is_qr($_,15)}0..32], [1,1,0,0,1,0,1,0,0,1,1,0,0,0,0,1,1,0,0,1,0,1,0,0,1,1,0,0,0,0,1,1,0], "is_qr(0..32,15)");
is(is_qr("2636542937688", "3409243234243"), 1, "2636542937688 is a qr mod 3409243234243");
