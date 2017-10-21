#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/ chebyshev_theta chebyshev_psi /;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
my $usegmp= Math::Prime::Util::prime_get_config->{'gmp'};
$use64 = 0 if $use64 && 18446744073709550592 == ~0;

my %chebyshev1 = (
       0 =>       0,
       1 =>       0,
       2 =>       0.693147180559945,
       3 =>       1.79175946922805,
       4 =>       1.79175946922805,
       5 =>       3.40119738166216,
     243 =>     226.593507136467,
  123456 =>  123034.091739914,
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

plan tests => scalar(keys %chebyshev1) + scalar(keys %chebyshev2);

###### first Chebyshev function
while (my($n, $c1) = each (%chebyshev1)) {
  cmp_closeto( chebyshev_theta(0+$n), $c1, 1e-9*abs($n), "chebyshev_theta($n)" );
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
