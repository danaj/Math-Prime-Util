#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/setbinop/;

plan tests => 8;

###### setbinop

is_deeply( [setbinop { $a+$b } [],[1,2,3]], [], "setbinop with an empty set" );

is_deeply( [setbinop { $a+$b } [2,4,6]],
           [4,6,8,10,12],
           "setbinop A+A" );

is_deeply( [setbinop { $a+$b } [2,4,6,8],[3,5,7]],
           [5,7,9,11,13,15],
           "setbinop A+B" );

is_deeply( [setbinop { $a-$b } [1,2,3],[5,7,9]],
           [-8,-7,-6,-5,-4,-3,-2],
           "setbinop A-B" );

is_deeply( [setbinop { $b-$a } [1,2,3],[5,7,9]],
           [2,3,4,5,6,7,8],
           "setbinop B-A" );

is_deeply( [setbinop { $a+2*$b } [1,2,3],[5,7,9]],
           [11, 12, 13, 15, 16, 17, 19, 20, 21],
           "setbinop A+2B" );

is_deeply( [setbinop { ($a+$b) % 4 } [3,7,14]],
           [0,1,2],
           "setbinop A+A mod 4" );

# from Granville and Roesler
{
  my @s = grep { /^[124]{2}$/ } 1..99;
  is(scalar(setbinop { $a-$b } \@s), 7**2, "[124]{2} has 3^2 elements, A-A has 7^2 elements");
}
