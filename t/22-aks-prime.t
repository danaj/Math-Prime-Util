#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_aks_prime/;

my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
my $usegmp = Math::Prime::Util::prime_get_config->{'gmp'};

plan tests =>   6   # range
              + 1   # small number
              + 2   # medium numbers
              + 1*$extra
              + 4   # Numbers for BERN41 version
              + 0;

# Note: AKS testing is *extremely* sparse due to its lack of speed.
#       This does almost nothing to test whether AKS is working properly.
#
# If you are concerned about AKS correctness, you really need to use
# the xt/primality-aks.pl test.

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

# These run the polynomial test with the V6 code
is( is_aks_prime(69197), 1, "is_aks_prime(69197) is true" );
is( is_aks_prime(69199), 0, "is_aks_prime(69199) is false" );

if ($extra) {
  # A composite (product of two 3-digit primes)
  #is( is_aks_prime(370781), 0, "is_aks_prime(370781) is false" );
  # A composite (product of two 4-digit primes)
  #is( is_aks_prime(37809463), 0, "is_aks_prime(37809463) is false" );
  # A composite (product of two 5-digit primes)
  is( is_aks_prime(1262952907), 0, "is_aks_prime(1262952907) is false" );
}

# With Bernstein 4.1, we need to use different numbers.
# It is highly unlikely that a composite will even run the polynomial test.
{
  is(is_aks_prime(101), 1, "is_aks_prime(101)=1");
  is(is_aks_prime(15481), 0, "is_aks_prime(15481)=0");
  is(is_aks_prime(12109), 1, "is_aks_prime(12109)=1");
  SKIP: {
    skip "Skipping 29-bit prime in PP",1 unless $usexs;
    # 0.006 s   XS 64-bit
    # 0.08  s   XS 32-bit
    # 0.02  s   GMP
    # 0.8   s   64-bit PP
    # 7m 50 s   32-bit PP (ouch)
    is(is_aks_prime(536891893), 1, "is_aks_prime(536891893)=1");
  }
}
