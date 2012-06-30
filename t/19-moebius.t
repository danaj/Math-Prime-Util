#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/moebius euler_phi/;

my @moeb_vals = (qw/ 1 -1 -1 0 -1 1 -1 0 0 1 -1 0 -1 1 1 0 -1 0 -1 0 /);
my %mertens = (
        1 =>    1,
        2 =>    0,
        3 =>   -1,
        4 =>   -1,
        5 =>   -2,
       10 =>   -1,
      100 =>    1,
     1000 =>    2,
    10000 =>  -23,
   100000 =>  -48,
#  1000000 =>  212,
# 10000000 => 1037,
);

my %totients = (
          1 => 1,
          2 => 1,
         20 => 8,
         36 => 12,
     123456 => 41088, 
     123457 => 123456,
  123456789 => 82260072,
);

plan tests => 0 + 1
                + scalar @moeb_vals
                + scalar(keys %mertens)
                + scalar(keys %totients);

ok(!eval { moebius(0); }, "moebius(0)");

my $i = 1;
foreach my $m (@moeb_vals) {
  is( moebius($i), $m, "moebius($i) == $m" );
  $i++;
}

while (my($n, $mertens) = each (%mertens)) {
  my $M = 0;
  $M += moebius($_) for (1 .. $n);
  is( $M, $mertens, "Mertens($n) == $mertens" );
}

while (my($n, $phi) = each (%totients)) {
  is( euler_phi($n), $phi, "euler_phi($n) == $phi" );
}
