#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/legendre_phi/;

#my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
#my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
#my $usegmp= Math::Prime::Util::prime_get_config->{'gmp'};
#my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
#$use64 = 0 if $use64 && 18446744073709550592 == ~0;

my @legendre_sums = (
  [ 0,  92372, 0],
  [ 5,  15, 1],
  [ 89, 4, 21 ],
  [ 46, 4, 11 ],
  [ 47, 4, 12 ],
  [ 48, 4, 12 ],
  [ 52, 4, 12 ],
  [ 53, 4, 13 ],
  [10000, 5, 2077],
  [526, 7, 95],
  [588, 6, 111],
  [100000, 5, 20779],
  [5882, 6, 1128],
  [100000, 7, 18053],
  [10000, 8, 1711],
  [1000000, 168, 78331],
  [800000, 213, 63739],
  [4000,255,296],
);

plan tests => scalar(@legendre_sums);

###### Legendre phi
foreach my $r (@legendre_sums) {
  my($x, $a, $exp) = @$r;
  is( legendre_phi($x, $a), $exp, "legendre_phi($x,$a) = $exp" );
}
