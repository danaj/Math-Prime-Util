#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/ chebyshev_theta chebyshev_psi /;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
my $usegmp= Math::Prime::Util::prime_get_config->{'gmp'};

my %chebyshev1 = (
         0 =>         0,
         1 =>         0,
         2 =>         0.693147180559945,
         3 =>         1.79175946922805,
         4 =>         1.79175946922805,
         5 =>         3.40119738166216,
       243 =>       226.59350713646702688959187294159019343,
      2010 =>      1947.44160163822684808718124934809069276,
      2011 =>      1955.04798902799950013394060154173889009,
    123456 =>    123034.09173991353285470434178654771360,
  69201234 =>  69192285.907794531021524894279820294952,
);
my %chebyshev2 = (
         0 =>       0,
         1 =>       0,
         2 =>       0.693147180559945,
         3 =>       1.79175946922805,
         4 =>       2.484906649788,
         5 =>       4.0943445622221,
       243 =>     245.274469978683,
    123456 =>  123435.148054491
);
if ($extra) {
  $chebyshev1{1234567} = 1233272.80087825;
  $chebyshev2{1234567} = 1234515.17962833;
}
if (!$usexs && !$extra) {
  delete $chebyshev1{$_} for grep { $_ > 50000 } keys %chebyshev1;
  delete $chebyshev2{$_} for grep { $_ > 50000 } keys %chebyshev2;
}

plan tests => scalar(keys %chebyshev1) + scalar(keys %chebyshev2) + 3;

###### first Chebyshev function
while (my($n, $c1) = each (%chebyshev1)) {
  cmp_closeto( chebyshev_theta(0+$n), $c1, 1e-9*abs($n), "chebyshev_theta($n)" );
}
SKIP: {
  skip "Chebyshev checkpoint tests require 64-bit XS", 3 unless $usexs && $use64;
  my $n = 10000000000000;
  my $theta = 9999996988293.03419965318214160283774190;
  cmp_closeto(chebyshev_theta($n), $theta, 0.004, "theta at a checkpoint");
  cmp_closeto(chebyshev_theta($n+36), $theta, 0.004, "theta before the next prime");
  cmp_closeto(chebyshev_theta($n+37), 9999996988322.96780586210843549507162, 0.004,
             "theta includes the first prime after the checkpoint");
}
###### second Chebyshev function
while (my($n, $c2) = each (%chebyshev2)) {
  cmp_closeto( chebyshev_psi(0+$n), $c2, 1e-9*abs($n), "chebyshev_psi($n)" );
}


sub cmp_closeto {
  my $got = shift;
  my $expect = shift;
  my $tolerance = shift;
  my $message = shift;
  cmp_ok( abs($got - $expect), '<=', $tolerance, $message );
}
