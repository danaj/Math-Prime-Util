#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_prime miller_rabin/;

my $use64 = Math::Prime::Util::_maxbits > 32;
my $extra = defined $ENV{RELEASE_TESTING} && $ENV{RELEASE_TESTING};

#plan tests => 6 + 19 + 3573 + (5 + 29 + 22 + 23 + 16) + 15 + 27
#              + ($use64 ? 5+1 : 0)
#              + ($extra ? 6 : 0)
#              + (($extra && $use64) ? 19 : 0);

map { ok(miller_rabin($_,2), "Pseudoprime (base 2) $_ passes MR" ) }
  qw/2047 3277 4033 4681 8321 15841 29341 42799 49141 52633 65281 74665
     80581 85489 88357 90751/;

my $mr2fail = 0;
for (2 .. 4032) {
  next if $_ == 2047 || $_ == 3277;
  if (is_prime($_)) {
    if (!miller_rabin($_,2)) { $mr2fail = $_; last; }
  } else {
    if (miller_rabin($_,2))  { $mr2fail = $_; last; }
  }
}
is($mr2fail, 0, "miller_rabin matched is_prime for 2-4032 (skipping 2047 & 3277)");

map { ok(miller_rabin($_,3), "Pseudoprime (base 3) $_ passes MR" ) }
  qw/121 703 1891 3281 8401 8911 10585 12403 16531 18721 19345 23521 31621
     44287 47197 55969 63139 74593 79003 82513 87913 88573 97567/;

map { ok(miller_rabin($_,5), "Pseudoprime (base 5) $_ passes MR" ) }
  qw/781 1541 5461 5611 7813 13021 14981 15751 24211 25351 29539 38081
     40501 44801 53971 79381/;

done_testing();
