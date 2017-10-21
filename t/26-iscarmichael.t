#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/ is_carmichael is_quasi_carmichael /;

#my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
#my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
#my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
my $usegmp= Math::Prime::Util::prime_get_config->{'gmp'};
#$use64 = 0 if $use64 && 18446744073709550592 == ~0;

plan tests => 2 + 4;

###### is_carmichael
is_deeply( [grep { is_carmichael($_) } 1 .. 20000],
           [561,1105,1729,2465,2821,6601,8911,10585,15841],
           "Carmichael numbers to 20000" );
SKIP: {
  skip "Skipping large Carmichael", 1 unless $usegmp;
  ok( is_carmichael("341627175004511735787409078802107169251"), "Large Carmichael" );
}

###### is_quasi_carmichael
is_deeply( [grep { is_quasi_carmichael($_) } 1 .. 400],
           [35,77,143,165,187,209,221,231,247,273,299,323,357,391,399],
           "Quasi-Carmichael numbers to 400" );
is( scalar(grep { is_quasi_carmichael($_) } 1 .. 5000),
           95,
           "95 Quasi-Carmichael numbers under 5000" );
is(is_quasi_carmichael(5092583), 1, "5092583 is a Quasi-Carmichael number with 1 base");
is(is_quasi_carmichael(777923), 7, "777923 is a Quasi-Carmichael number with 7 bases");
