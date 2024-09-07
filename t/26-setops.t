#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/setunion setintersect setminus setdelta toset
                         powint addint/;
use Math::BigInt;

my $bi1 = Math::BigInt->new("98198086365677506205371483123156488634");
my $bi2 = Math::BigInt->new("59724578844314338843734830435499460367");
my $bi3 = Math::BigInt->new("606571739116749108206251582180042583662");
my $pr = "1844674407370955161";

my @sets = (
  [ [0,1,2], [2,3,4], "simple unsigned", [0,1,2,3,4], [2], [0,1], [0,1,3,4] ],
  [ [7,1,3,5,1], [3,7,8,3,9], "simple unsigned unsorted with dups", [1,3,5,7,8,9], [3,7], [1,5], [1,5,8,9] ],
  [ [], [1,2,3], "empty first list", [1,2,3], [], [], [1,2,3] ],
  [ [1,2,3], [], "empty second list", [1,2,3], [], [1,2,3], [1,2,3] ],
  [ [], [], "empty lists", [], [], [], [] ],
  [ [-5..1],[-1..5], "signed overlap", [-5..5],[-1..1],[-5..-2],[-5..-2,2..5] ],
  [ ["9223372036854775808","9223372036854775807"],["9223372036854775807","9223372036854775810"], "too big for IV", ["9223372036854775807","9223372036854775808","9223372036854775810"], ["9223372036854775807"], ["9223372036854775808"],["9223372036854775808","9223372036854775810"] ],
  [ ["10223372036854775808","-9223372036854775807"],["-9223372036854775807","10223372036854775810"], "range bigger than IV or UV", ["-9223372036854775807","10223372036854775808","10223372036854775810"], ["-9223372036854775807"], ["10223372036854775808"],["10223372036854775808","10223372036854775810"] ],
  [ [$bi1,$bi2],[$bi3,$bi2], "bigints", [$bi2,$bi1,$bi3],[$bi2],[$bi1],[$bi1,$bi3] ],
  # This tests that we sort correctly even when given strings that Perl
  # doesn't compare properly.  This will fail if we try to use sort {$a<=>$b}.
  [ [$pr.6,$pr.5,$pr.7,$pr.4,$pr.8], [$pr.6,$pr.5,$pr.3], "mix 64-bit and 65-bit as strings", [map { $pr.$_ } 3..8], [$pr.5,$pr.6], [$pr.4,$pr.7,$pr.8], [$pr.3,$pr.4,$pr.7,$pr.8] ],
);

plan tests => 2 + 4*scalar(@sets) + 6;

###### some specific tests

is_deeply([setunion([1,2,3],[-5,10,-11])],[-11,-5,1,2,3,10],"union signed properly sorted");
is( scalar setdelta([7,1,3,5,1], [3,7,8,3,9]), 4, "scalar setdelta with dups yields same as array" );

###### union, intersect, difference, symmetric difference

for my $info (@sets) {
  my($x,$y,$str,$aunion,$ainter,$aminus,$adelta) = @$info;
  is_deeply( [setunion($x,$y)], $aunion, "$str   union" );
  is_deeply( [setintersect($x,$y)], $ainter, "$str   intersect" );
  is_deeply( [setminus($x,$y)], $aminus, "$str   minus" );
  is_deeply( [setdelta($x,$y)], $adelta, "$str   delta" );
}

###### toset
is_deeply( [toset([])],[],"toset: empty list" );
is_deeply( [toset([1])],[1],"toset: one value" );
is_deeply( [toset([3,-2,1,3,3,-14])],[-14,-2,1,3],"toset: simple" );
is_deeply( [toset([qw/1 -2147483647 3 2 2147483648/])],
           [qw/-2147483647 1 2 3 2147483648/],
           "toset: 32-bit mix of sign and unsigned" );
is_deeply( [toset([qw/1 -9223372036854775807 3 2 9223372036854775808/])],
           [qw/-9223372036854775807 1 2 3 9223372036854775808/],
           "toset: 64-bit mix of sign and unsigned" );
{
  my $b = powint(2,129);
  my @L = map { addint($b,$_) } (-2,3,0,0,-2,0,3);
  my @R = map { addint($b,$_) } (-2,0,3);
  is_deeply( [toset(\@L)], \@R, "toset: 129-bit unsigned inputs" );
}


###### setcomplement

#is_deeply([setcomplement([1,3,5],0,8)], [0,2,4,6,7,8], "setcomplement 0..8");
#is_deeply([setcomplement([1,3,5],-2,8)], [-2,-1,0,2,4,6,7,8], "setcomplement -2..8");
#is_deeply([setcomplement([],0,5)], [0..5], "setcomplement 0..5 with empty set");
#is_deeply([setcomplement([],0,5)], [0..5], "setcomplement 0..5 with empty set");
#is_deeply([setcomplement([0..5],6,8)], [6..8], "setcomplement 0..5 with 6..8");
#is_deeply([setcomplement(["-9223372036854775810","-9223372036854775806"],"-9223372036854775810","-9223372036854775806")], ["-9223372036854775809","-9223372036854775808","-9223372036854775807"], "setcomplement crossing IV");
#is_deeply([setcomplement(["18446744073709551614","18446744073709551617"],"18446744073709551614","18446744073709551617")], ["18446744073709551615","18446744073709551616"], "setcomplement crossing UV");
