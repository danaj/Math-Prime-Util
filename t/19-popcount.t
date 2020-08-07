#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/hammingweight prime_set_config/;

#my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
#my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
my $usegmp= Math::Prime::Util::prime_get_config->{'gmp'};
#my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
#$use64 = 0 if $use64 && 18446744073709550592 == ~0;

my @popcounts = (
  [0, 0],
  [1, 1],
  [2, 1],
  [3, 2],
  [452398, 12],
  [-452398, 12],
  [4294967295, 32],
  ["777777777777777714523989234823498234098249108234236", 83],
  ["65520150907877741108803406077280119039314703968014509493068998974809747144933", 118],
);

plan tests => scalar(@popcounts) * 2;

###### hammingweight
foreach my $r (@popcounts) {
  my($n, $exp) = @$r;
  is( hammingweight($n), $exp, "hammingweight($n) = $exp" );
}

###### Turn off gmp and try
SKIP: {
  skip "No need to test non-GMP when not using GMP", scalar(@popcounts) unless $usegmp;
  prime_set_config(gmp=>0);
  foreach my $r (@popcounts) {
    my($n, $exp) = @$r;
    is( hammingweight($n), $exp, "non-GMP hammingweight($n) = $exp" );
  }
}
