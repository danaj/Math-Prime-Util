#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/ znprimroot is_primitive_root /;

#my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
#my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
#my $usegmp= Math::Prime::Util::prime_get_config->{'gmp'};
$use64 = 0 if $use64 && 18446744073709550592 == ~0;

my %primroots = (
   -11 => 2,
     0 => undef,
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
}

plan tests => 0
                + scalar(keys %primroots) + 1  # znprimroot
                + scalar(keys %primroots) + 3  # is_primitive_root
                ;

###### znprimroot
while (my($n, $root) = each (%primroots)) {
  is( znprimroot(0+$n), $root, "znprimroot($n) == " . ((defined $root) ? $root : "<undef>") );
}
is( znprimroot("-100000898"), 31, "znprimroot(\"-100000898\") == 31" );
# I don't think we should rely on this parsing correctly.
#is( znprimroot("+100000898"), 31, "znprimroot(\"+100000898\") == 31" );

###### is_primitive_root
while (my($n, $root) = each (%primroots)) {
  if (defined $root) {
    is( is_primitive_root(0+$root,0+$n), 1, "$root is a primitive root mod $n" );
  } else {
    is( is_primitive_root(2,0+$n), 0, "2 is not a primitive root mod $n" );
  }
}
is(is_primitive_root(19,191), 1, "19 is a primitive root mod 191");
is(is_primitive_root(13,191), 0, "13 is not a primitive root mod 191");
is(is_primitive_root(35,982), 0, "35 is not a primitive root mod 982");
