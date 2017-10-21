#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/valuation/;

#my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
#my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
#my $usegmp= Math::Prime::Util::prime_get_config->{'gmp'};
#my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
#$use64 = 0 if $use64 && 18446744073709550592 == ~0;

my @valuations = (
  [-4,2, 2],
  [0,0, 0],
  [1,0, 0],
  [96552,6, 3],
  [1879048192,2, 28],
);

plan tests => scalar(@valuations);

###### valuation
foreach my $r (@valuations) {
  my($n, $k, $exp) = @$r;
  is( valuation($n, $k), $exp, "valuation($n,$k) = $exp" );
}
