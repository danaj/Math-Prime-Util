#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/kronecker/;

#my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
#my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
#my $usegmp= Math::Prime::Util::prime_get_config->{'gmp'};
my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
$use64 = 0 if $use64 && 18446744073709550592 == ~0;

my @kroneckers = (
  [ 109981, 737777,  1],
  [ 737779, 121080, -1],
  [-737779, 121080,  1],
  [ 737779,-121080, -1],
  [-737779,-121080, -1],
  [12345,331,-1],
  [1001,9907,-1],
  [19,45,1],
  [8,21,-1],
  [5,21,1],
  [5,1237,-1],
  [10, 49, 1],
  [123,4567,-1],
  [3,18,0], [3,-18,0],
  [-2, 0, 0],  [-1, 0, 1],  [ 0, 0, 0],  [ 1, 0, 1],  [ 2, 0, 0],
  [-2, 1, 1],  [-1, 1, 1],  [ 0, 1, 1],  [ 1, 1, 1],  [ 2, 1, 1],
  [-2,-1,-1],  [-1,-1,-1],  [ 0,-1, 1],  [ 1,-1, 1],  [ 2,-1, 1],
  # Some cases trying to make sure we're not turning UVs into IVs
  [ 3686556869,  428192857,  1],
  [-1453096827,  364435739, -1],
  [ 3527710253, -306243569, 1],
  [-1843526669, -332265377, 1],
  [  321781679, 4095783323, -1],
  [  454249403,  -79475159, -1],
);
if ($use64) {
  push @kroneckers, [17483840153492293897, 455592493, 1];
  push @kroneckers, [-1402663995299718225, 391125073, 1];
  push @kroneckers, [16715440823750591903, -534621209, -1];
  push @kroneckers, [13106964391619451641,16744199040925208803, 1];
  push @kroneckers, [11172354269896048081,10442187294190042188,-1];
  push @kroneckers, [-5694706465843977004,9365273357682496999,-1];
}

plan tests => scalar(@kroneckers);

###### kronecker
foreach my $karg (@kroneckers) {
  my($a, $n, $exp) = @$karg;
  my $k = kronecker($a, $n);
  is( $k, $exp, "kronecker($a, $n) = $exp" );
}
