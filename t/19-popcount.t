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

my @inputs = (qw/0 1 2 3 15 -15 452398 -452398 4294967295 506198202011006650616143 777777777777777714523989234823498234098249108234236 65520150907877741108803406077280119039314703968014509493068998974809747144933/);
my @popcnt = (0,1,1,2,4,4,12,12,32,45,83,118);

plan tests => 2;

###### hammingweight
is_deeply( [map { hammingweight($_) } @inputs], \@popcnt, "hammingweight for various inputs" );

###### Turn off gmp and try
SKIP: {
  skip "No need to test non-GMP when not using GMP", 1 unless $usegmp;
  prime_set_config(gmp=>0);
  is_deeply( [map { hammingweight($_) } @inputs], \@popcnt, "non-GMP hammingweight for various inputs" );
}
