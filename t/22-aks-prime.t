#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_aks_prime/;

my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $ispp = !Math::Prime::Util::prime_get_config->{xs};

plan tests =>   6   # range
              + 1   # small number
              + 2   # medium numbers
              + 1*$extra
              + 0;

ok(!eval { is_aks_prime(undef); }, "is_prime(undef)");
ok( is_aks_prime(2),  '2 is prime');
ok(!is_aks_prime(1),  '1 is not prime');
ok(!is_aks_prime(0),  '0 is not prime');
ok(!is_aks_prime(-1), '-1 is not prime');
ok(!is_aks_prime(-2), '-2 is not prime');

# Simple number (cought by sqrt test)
is( is_aks_prime(877), 1, "is_aks_prime(877) is true" );

# Perhaps let them know this is probably not a hung test?
# This runs in milliseconds on an i3930K, but 1.5 minutes on an UltraSPARC.
# These are the smallest numbers that actually run the code, so I don't know
# how to make them run any faster.  On the 32-bit UltraSPARC, it's the mulmod
# that is painfully slow.

#diag "Unfortunately these tests are very slow.";

SKIP: {
  # If we're pure Perl, then this is definitely too slow.
  # Arguably we should check to see if they have the GMP code.
  skip "Skipping PP AKS on PP -- just too slow.", 1 if $ispp;
  # If we have 64-bit available in the compiler (e.g. uint64_t), this can
  # still be quite fast.  However for pretty much everyone else, this is
  # just far too slow for running in a test suite.
  skip "Skipping PP AKS on 32-bit -- just too slow.", 1 if !$use64;
  # The first number that makes it past the sqrt test to actually run.
  is( is_aks_prime(69197), 1, "is_aks_prime(69197) is true" );
}

# A small composite that runs the real primality test.
is( is_aks_prime(69199), 0, "is_aks_prime(69199) is false" );

if ($extra) {
  # A composite (product of two 5-digit primes)
  is( is_aks_prime(1262952907), 0, "is_aks_prime(1262952907) is false" );
}
