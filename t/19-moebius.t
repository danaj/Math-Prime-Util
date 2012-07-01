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
     123456 => 41088, 
     123457 => 123456,
  123456789 => 82260072,
);
my @A000010 = (0,1,1,2,2,4,2,6,4,6,4,10,4,12,6,8,8,16,6,18,8,12,10,22,8,20,12,18,12,28,8,30,16,20,16,24,12,36,18,24,16,40,12,42,20,24,22,46,16,42,20,32,24,52,18,40,24,36,28,58,16,60,30,36,32,48,20,66,32,44);
@totients{0..$#A000010} = @A000010;

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
