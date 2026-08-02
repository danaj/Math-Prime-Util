#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/random_factored_integer irand factor vecprod/;

my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
#my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
#my $maxbits = $use64 ? 64 : 32;

plan tests => 7;

########

my($n, $factors) = random_factored_integer(1000000);
ok($n > 0, "random_factored_integer did not return 0");
ok($n <= 1000000, "random_factored_integer in requested range");
# random_factored_integer returns factors in sorted order
is_deeply( $factors, [factor($n)], "factors match factor routine");
is( vecprod(@$factors), $n, "product of factors = n");

SKIP: {
  skip "large native random_factored_integer requires 64-bit", 3 unless $use64;
  my $limit = "1000000000000";
  my($large, $large_factors) = random_factored_integer($limit);
  ok($large > 0 && $large <= $limit, "large random_factored_integer in requested range");
  is_deeply($large_factors, [factor($large)], "large result factors match factor routine");
  is(vecprod(@$large_factors), $large, "large result is the product of its factors");
}
