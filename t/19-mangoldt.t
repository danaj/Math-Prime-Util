#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/exp_mangoldt/;

#my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
#my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
#my $usegmp= Math::Prime::Util::prime_get_config->{'gmp'};
#my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
#$use64 = 0 if $use64 && 18446744073709550592 == ~0;

my %mangoldt = (
#-13 => 1,    # Non-negative inputs
  0 => 1,
  1 => 1,
  2 => 2,
  3 => 3,
  4 => 2,
  5 => 5,
  6 => 1,
  7 => 7,
  8 => 2,
  9 => 3,
 10 => 1,
 11 => 11,
 25 => 5,
 27 => 3,
 399981 => 1,
 399982 => 1,
 399983 => 399983,
 823543 => 7,
 83521 => 17,
 130321 => 19,
);

plan tests => scalar(keys %mangoldt);

###### Exponential of von Mangoldt
while (my($n, $em) = each (%mangoldt)) {
  is( exp_mangoldt(0+$n), $em, "exp_mangoldt($n) == $em" );
}
