#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/prime_count prime_count_lower prime_count_upper prime_count_approx/;

my $use64 = Math::Prime::Util::_maxbits > 32;
my $extra = defined $ENV{RELEASE_TESTING} && $ENV{RELEASE_TESTING};

plan tests => 10*3 + ($extra ? 10 : 7) + ($use64 ? 2*9 : 0);

my %pivals32 = (
                  1 => 0,
                 10 => 4,
                100 => 25,
               1000 => 168,
              10000 => 1229,
             100000 => 9592,
            1000000 => 78498,
           10000000 => 664579,
          100000000 => 5761455,
         1000000000 => 50847534,
);
my %pivals64 = (
        10000000000 => 455052511,
       100000000000 => 4118054813,
      1000000000000 => 37607912018,
     10000000000000 => 346065536839,
    100000000000000 => 3204941750802,
   1000000000000000 => 29844570422669,
  10000000000000000 => 279238341033925,
 100000000000000000 => 2623557157654233,
1000000000000000000 => 24739954287740860,
);
while (my($n, $pin) = each (%pivals32)) {
  cmp_ok( prime_count_upper($n), '>=', $pin, "Pi($n) <= upper estimate" );
  cmp_ok( prime_count_lower($n), '<=', $pin, "Pi($n) >= lower estimate" );
  if ( ($n <= 2000000) || $extra ) {
    is( prime_count($n), $pin, "Pi($n) = $pin" );
  }
  my $approx_range = abs($pin - prime_count_approx($n));
  my $range_limit = 1100;
  cmp_ok( $approx_range, '<=', $range_limit, "prime_count_approx($n) within $range_limit");
}
if ($use64) {
  while (my($n, $pin) = each (%pivals64)) {
    cmp_ok( prime_count_upper($n), '>=', $pin, "Pi($n) <= upper estimate" );
    cmp_ok( prime_count_lower($n), '<=', $pin, "Pi($n) >= lower estimate" );
  }
}

