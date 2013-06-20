#!/usr/bin/env perl
use strict;
use warnings;

# I found these issues when doing some testing of is_provable_prime.  When
# bignum is loaded, we get some strange behavior.  There are two fixes for
# it in the code:
#  1) make sure every divide and bdiv is coerced back to an integer.
#  2) turn off upgrade in input validation.
# The second method in theory is all that is needed.

use Math::Prime::Util qw/:all/;
use bignum;

use Test::More tests => 2;

if ($] < 5.008) {
  diag "A prototype warning was expected with old, old Perl";
}

my $n = 100199294509778143137521762187425301691197073534078445671945250753109628678272;
# 2 2 2 2 2 2 2 3 7 509 277772399 263650456338779643073784729209358382310353002641378210462709359

my @partial_factor = Math::Prime::Util::PP::prho_factor(100199294509778143137521762187425301691197073534078445671945250753109628678272, 5);

my @expected_factors =
   map { ($_ <= 4294967295 && ref($_)) ? int($_->bstr) : $_ }
   (2,2,2,2,2,2,2,3,7,37276523255125797298185179385202865212498911284999421752955822452793760669);

is_deeply( \@partial_factor, \@expected_factors,
           "PP prho factors correctly with 'use bignum'" );

# The same thing happens in random primes, PP holf factoring,
# PP is_provable_primes, and possibly elsewhere

################################################################################

# Here is another test case that has to do with reference counting
# in the XS subroutine callback code.
SKIP: {
  skip "No MPU::GMP, skipping callback test",1 unless prime_get_config->{'gmp'};
  my $n = 10**1200+5226;
  my $expect = $n+1;
  is(next_prime($n), $expect, "next_prime(10^1200+5226) = 10^1200+5227");
}
