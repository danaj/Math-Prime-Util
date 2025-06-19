#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/binomialmod/;
#use Math::Prime::Util qw/binomial modint/;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $usexs  = Math::Prime::Util::prime_get_config->{'xs'};

my @tests = (
  [0,0,7, 1],
  [0,1,7, 0],
  [0,2,7, 0],
  [3,0,7, 1],
  [7,5,11, 10],
  [950,100,123456, 24942],
  [950,100,7, 2],
  [8100,4000,1155, 924],
  [950,100,1000000007, 640644226],
  [189,34,877, 81],
  [189,34,253009, 47560],
  [189,34,36481, 14169],
  [1900,17,41, 0],
  [5000,654,101223721, 59171352],

  [-112,5,351, 313],
  [-189,34,877, 141],
  [-23,-29,377, 117],

  [189,-34,877, 0],
);

if ($usexs) {
  push @tests, [100000000,87654321,1005973, 937361],
               [100000000,7654321,1299709, 582708],
               [100000000,7654321,12345678, 4152168],
               [100000,7654,32768, 12288],
               [100000,7654,196608, 110592],
               [100000,7654,101223721, 5918452];
}

if ($usexs && $extra) {
  push @tests, [100000000,7654321,32768, 24576],
               [100000000,7654321,196608, 122880],
               [100000000,7654321,101223721, 5463123];
}

plan tests => 1
            + scalar(@tests)
            ;

for my $data (@tests) {
  my($n,$k,$m,$exp) = @$data;
  is( binomialmod($n,$k,$m), $exp, "binomialmod($n,$k,$m) = $exp" );
  #is( modint(binomial($n,$k),$m), $exp, "binomial($n,$k) mod $m = $exp" );
}

{
  my $sum = 0;
  for my $p (1 .. 27) {
    for my $n (1..40) {
      for my $k (1 .. $n) {
        $sum += binomialmod($n,$k,$p);
      }
    }
  }
  # s=0;for(p=1,27,for(n=1,40,for(k=1,n,s+=lift(binomod(n,k,p))))); s
  is( $sum, 99531, "Small binomialmod works" );
}
