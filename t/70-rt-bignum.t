#!/usr/bin/env perl
use strict;
use warnings;
use Test::More;

# I found these issues when doing some testing of is_provable_prime.  When
# bignum is loaded, we get some strange behavior.  There are two fixes for
# it in the code:
#  1) make sure every divide and bdiv is coerced back to an integer.
#  2) turn off upgrade in input validation.
# The second method in theory is all that is needed.  We do both.

# Limit to AUTHOR_TESTING only.
# There is too much danger of effectively failing by hang.
#
# We try to discourage use of bignum because of the way it changes
# behavior of integers to be more like floats.
#
# This test *should* work, and we'd prefer to find it in testing than prod.
#
# At least one knowledgeable person has indicated that it's perfectly
# natural the 'bignum' module will wreak havoc with us and this test is
# going beyond what's needed.

BEGIN {
  unless ($ENV{AUTHOR_TESTING}) {
    plan( skip_all => 'these tests are for author testing' );
  }
}

use Math::Prime::Util qw/:all/;
use Math::Prime::Util::PP;
use bignum;

plan tests => 3+1;

if ($] < 5.008) {
  diag "A prototype warning was expected with old, old Perl";
}

my $n = 100199294509778143137521762187425301691197073534078445671945250753109628678272;
# 2 2 2 2 2 2 2 3 7 509 277772399 263650456338779643073784729209358382310353002641378210462709359

my @partial_factor = Math::Prime::Util::PP::prho_factor(100199294509778143137521762187425301691197073534078445671945250753109628678272, 5);

# Don't assume a fixed set of factors beyond 2/3/5 will be found.
{
  my @S = @partial_factor[0..7];
  my @L = @partial_factor[8..$#partial_factor];
  is_deeply(\@S, [2,2,2,2,2,2,2,3], "PP prho with 'use bignum' small factors are correct");
  my $expf = [qw/7 509 3563 277772399 1944406793 141386151091 989703057637/];
  ok(setcontainsany(\@L, $expf), "We found one of the expected small factors or products");
  is("".vecprod(@L), "260935662785880581087296255696420056487492378994995952270690757169556324683", "product of factors is the input");
}

# The same thing happens in random primes, PP holf factoring,
# PP is_provable_primes, and possibly elsewhere

################################################################################

# Here is another test case that has to do with reference counting
# in the XS subroutine callback code.  Perl 5.8.x [x < 8] will get lost
# and just exit with no message.
SKIP: {
  skip "No MPU::GMP, skipping callback test",1 unless prime_get_config->{'gmp'};
  my $n = 10**1200+5226;
  my $expect = $n+1;
  is(next_prime($n), $expect, "next_prime(10^1200+5226) = 10^1200+5227");
}
