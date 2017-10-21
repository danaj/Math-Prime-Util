#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/znorder/;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
my $usegmp= Math::Prime::Util::prime_get_config->{'gmp'};
my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
$use64 = 0 if $use64 && 18446744073709550592 == ~0;

my @mult_orders = (
  [1, 35, 1],
  [2, 35, 12],
  [4, 35, 6],
  [6, 35, 2],
  [7, 35],
  #[2,1000000000000031,81788975100],
  [1, 1, 1],
  [0, 0],
  [1, 0],
  [25, 0],
  [1, 1, 1],
  [19, 1, 1],
  [1, 19, 1],
  [2, 19, 18],
  [3, 20, 4],
  [57,1000000003,189618],
  [67,999999749,30612237],
  [22,999991815,69844],
  [10,2147475467,31448382],
  [141,2147475467,1655178],
  [7410,2147475467,39409],
  [31407,2147475467,266],
);
if ($use64) {
  push @mult_orders, [2, 2405286912458753, 1073741824];  # Pari #1031
}

plan tests => scalar(@mult_orders);

###### znorder
foreach my $moarg (@mult_orders) {
  my ($a, $n, $exp) = @$moarg;
  my $zn = znorder($a, $n);
  is( $zn, $exp, "znorder($a, $n) = " . ((defined $exp) ? $exp : "<undef>") );
}
