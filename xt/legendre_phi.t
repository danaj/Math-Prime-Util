#!/usr/bin/env perl
use warnings;
use strict;

use Test::More;
use Math::Prime::Util qw/nth_prime prime_count/;

my $x = shift || 50000;

my $sqrtx = int(sqrt($x));
my $pcx = prime_count($x);
my $pcsqrtx = prime_count($sqrtx);

my @a = 1 .. $x;

foreach my $a (0 .. $sqrtx+1) {
  if ($a > 0) {
    my $p = nth_prime($a);
    @a = grep { $_ % $p } @a;
  }
  my $expect = scalar @a;
  if ($a > $pcsqrtx) {
    is ( $expect, $pcx - $a + 1, "sieved phi($x,$a) = Pi($x) - $a + 1" );
  }
  my $phixa = Math::Prime::Util::_XS_legendre_phi($x, $a);
  is( $phixa, $expect, "Legendre phi($x,$a) = $expect" );
}
done_testing();
