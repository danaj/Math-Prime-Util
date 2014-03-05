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

# This test can take a very long time if mulmods are very slow (e.g. on
# UltraSPARC).  With the B+V improvements this should be fast enough for
# the little example that we are ok.

#diag "Unfortunately these tests are very slow.";

SKIP: {
  # If we're pure Perl, then this is definitely too slow.
  # Arguably we should check to see if they have the GMP code.
  skip "Skipping PP AKS on PP -- just too slow.", 1 if $ispp;
  # The least number that performs the full test with either implementation.
  is( is_aks_prime(69197), 1, "is_aks_prime(69197) is true" );
}

# A small composite that runs the real primality test.
is( is_aks_prime(69199), 0, "is_aks_prime(69199) is false" );

if ($extra) {
  # A composite (product of two 5-digit primes)
  is( is_aks_prime(1262952907), 0, "is_aks_prime(1262952907) is false" );
}
