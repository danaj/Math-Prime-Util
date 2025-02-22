#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/setunion setintersect setminus setdelta
                         setcontains
                         toset is_subset is_sidon_set is_sumfree_set
                         powint addint subint negint/;
use Math::BigInt;

my $bi1 = Math::BigInt->new("59724578844314338843734830435499460367");
my $bi2 = Math::BigInt->new("98198086365677506205371483123156488634");
my $bi3 = Math::BigInt->new("606571739116749108206251582180042583662");
my $pr = "1844674407370955161";

# Inputs, description, results:  union, intersection, minus, delta

# Generic numeric lists (not sorted, maybe duplicates).
my @vecs = (
  [ [7,1,3,5,1], [3,7,8,3,9], "simple unsigned unsorted with dups", [1,3,5,7,8,9], [3,7], [1,5], [1,5,8,9] ],
  [ ["9223372036854775808","9223372036854775807"],["9223372036854775807","9223372036854775810"], "too big for IV", ["9223372036854775807","9223372036854775808","9223372036854775810"], ["9223372036854775807"], ["9223372036854775808"],["9223372036854775808","9223372036854775810"] ],
  [ ["10223372036854775808","-9223372036854775807"],["-9223372036854775807","10223372036854775810"], "range bigger than IV or UV", ["-9223372036854775807","10223372036854775808","10223372036854775810"], ["-9223372036854775807"], ["10223372036854775808"],["10223372036854775808","10223372036854775810"] ],
  [ [$bi1,$bi2],[$bi3,$bi2], "bigints", [$bi2,$bi1,$bi3],[$bi2],[$bi1],[$bi1,$bi3] ],
  # This tests that we sort correctly even when given strings that Perl
  # doesn't compare properly.  This will fail if we try to use sort {$a<=>$b}.
  [ [$pr.6,$pr.5,$pr.7,$pr.4,$pr.8], [$pr.6,$pr.5,$pr.3], "mix 64-bit and 65-bit as strings", [map { $pr.$_ } 3..8], [$pr.5,$pr.6], [$pr.4,$pr.7,$pr.8], [$pr.3,$pr.4,$pr.7,$pr.8] ],
);

# Inputs in proper set form (numerically sorted with no duplicates).
my @sets = (
  [ [0,1,2], [2,3,4], "simple unsigned", [0,1,2,3,4], [2], [0,1], [0,1,3,4] ],
  [ [0,1,2], [0,2,3,4], "simple unsigned", [0,1,2,3,4], [0,2], [1], [1,3,4] ],
  [ [], [1,2,3], "empty first list", [1,2,3], [], [], [1,2,3] ],
  [ [1,2,3], [], "empty second list", [1,2,3], [], [1,2,3], [1,2,3] ],
  [ [], [], "empty lists", [], [], [], [] ],
  [ [-5..1],[-1..5], "signed overlap", [-5..5],[-1..1],[-5..-2],[-5..-2,2..5] ],
  [ ["-9223372036854775807","10223372036854775808"],["-9223372036854775807","10223372036854775810"], "range bigger than IV or UV", ["-9223372036854775807","10223372036854775808","10223372036854775810"], ["-9223372036854775807"], ["10223372036854775808"],["10223372036854775808","10223372036854775810"] ],
  # More sign overlap
  [ [-20,-16,-14,-12,-10,0,12,14], [-30,-18,-14,-11,-10,-8,1,13,14], "sign overlap", [qw/-30 -20 -18 -16 -14 -12 -11 -10 -8 0 1 12 13 14/], [qw/-14 -10 14/], [qw/-20 -16 -12 0 12/], [qw/-30 -20 -18 -16 -12 -11 -8 0 1 12 13/] ],
);

plan tests => 2 + 4*scalar(@sets) + 4*scalar(@vecs) + 6 + 6 + 2 + 2
            + 3+2+5+4+6+6+7+7+1;  # setcontains

###### some specific tests

is_deeply([setunion([1,2,3],[-11,-5,10])],[-11,-5,1,2,3,10],"setunion signed properly sorted");
is( scalar setdelta([7,1,3,5,1], [3,7,8,3,9]), 4, "scalar setdelta with dups yields same as array" );

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

###### union, intersect, difference, symmetric difference

for my $info (@sets) {
  my($x,$y,$str,$aunion,$ainter,$aminus,$adelta) = @$info;

  is_deeply( [setunion($x,$y)], $aunion,     "setunion     $str" );
  is_deeply( [setintersect($x,$y)], $ainter, "setintersect $str" );
  is_deeply( [setminus($x,$y)], $aminus,     "setminus     $str" );
  is_deeply( [setdelta($x,$y)], $adelta,     "setdelta     $str" );
}

# With unsorted and/or non-unique lists, the set functions still work,
# though possibly slower.
for my $info (@vecs) {
  my($x,$y,$str,$aunion,$ainter,$aminus,$adelta) = @$info;
  my($vunion, $vinter, $vminus, $vdelta) = (
  [setunion($x,$y)], [setintersect($x,$y)], [setminus($x,$y)], [setdelta($x,$y)]
  );
  ($aunion,$ainter,$aminus,$adelta,$vunion,$vinter,$vminus,$vdelta) =
      map { [sort @$_] }
      ($aunion,$ainter,$aminus,$adelta,$vunion,$vinter,$vminus,$vdelta);
  is_deeply( $vunion, $aunion, "setunion     $str" );
  is_deeply( $vinter, $ainter, "setintersect $str" );
  is_deeply( $vminus, $aminus, "setminus     $str" );
  is_deeply( $vdelta, $adelta, "setdelta     $str" );
}

###### setcomplement

#is_deeply([setcomplement([1,3,5],0,8)], [0,2,4,6,7,8], "setcomplement 0..8");
#is_deeply([setcomplement([1,3,5],-2,8)], [-2,-1,0,2,4,6,7,8], "setcomplement -2..8");
#is_deeply([setcomplement([],0,5)], [0..5], "setcomplement 0..5 with empty set");
#is_deeply([setcomplement([],0,5)], [0..5], "setcomplement 0..5 with empty set");
#is_deeply([setcomplement([0..5],6,8)], [6..8], "setcomplement 0..5 with 6..8");
#is_deeply([setcomplement(["-9223372036854775810","-9223372036854775806"],"-9223372036854775810","-9223372036854775806")], ["-9223372036854775809","-9223372036854775808","-9223372036854775807"], "setcomplement crossing IV");
#is_deeply([setcomplement(["18446744073709551614","18446744073709551617"],"18446744073709551614","18446744073709551617")], ["18446744073709551615","18446744073709551616"], "setcomplement crossing UV");

###### is_subset
is( is_subset([3,1],[1,2,3]), 1, "{1,3} is a subset of {1,2,3}" );
is( is_subset([],[1,2,3]), 1, "The empty list is a subset of any list" );
is( is_subset([],[0]), 1, "The empty list is a subset of any list" );
is_deeply( [
            is_subset([],[]),
            is_subset([0],[0]),
            is_subset([1,2,3],[1,2,3]),
            is_subset([~0],[~0]),
            is_subset(["-9223372036854775807"],["-9223372036854775807"]),
            is_subset(["-9223372036854775808"],["-9223372036854775808"]),
            is_subset(["18446744073709551615"],["18446744073709551615"]),
            is_subset([0,"18446744073709551615"],[0,"18446744073709551615"]),
            is_subset(["-9223372036854775808",0,"18446744073709551615"],["-9223372036854775808",0,"18446744073709551615"]),
           ], [1,1,1,1,1,1,1,1,1], "Every list is a subset of itself" );
is_deeply( [
            is_subset([14,17,18],[qw/4 12 14 17 18 19 20/]),
           ], [1], "Test some subsets");
is_deeply( [
            is_subset([0],[1,2,3]),
            is_subset([~0],[1,2,3]),
            is_subset(["18446744073709551615"],["-9223372036854775808"]),
            is_subset(["18446744073709551615"],["-1"]),
            is_subset([2,3,5,7],[2,3,5,8]),
            is_subset([14,17,18],[qw/2 5 8 11 12 13 17 18 20/]),
           ], [0,0,0,0,0,0], "Test some non-subsets");

###### is_sidon_set
my @sidons = (
  [], [0], [0,1], [0,1,3], [0,1,4,6], [0,1,4,9,11],
  [0,1,4,10,12,17], [0,1,8,11,13,17],
  [qw/239 106 56 53 161/],[qw/9 10 1 3 14/],
  [qw/18446744073709551614 18446744073709551612 18446744073709551606/],
  [qw/0 5 20 51 57 83 136 169 196 292 425 434 544 586 786 910 1016 1187 1210 1228 1350 1369 1405 1453 1507 1760 1785 1850 1920 1964 2130 2223 2237 2318 2352 2390 2429 2439 2533 2601 2609 2622 2683 2808 2848 2870 2872 2917 2928 2945 3016 3045 3075 3229 3230 3321 3356 3576 3579 3677 3684 3727 3901 3905 3917/],
);
my @nonsidons = (
  [2,7,12], [1,10,11,12], [2,5,13,21,26], [0,1,4,10,12,16],
  [-1], [0,-1],[-9,17],
);
is_deeply( [map { is_sidon_set($_) } @sidons], [map { 1 } 0..$#sidons], "is_sidon_set with Sidon sets" );
is_deeply( [map { is_sidon_set($_) } @nonsidons], [map { 0 } 0..$#nonsidons], "is_sidon_set with non-Sidon sets" );

###### is_sumfree_set
my @sf = (
  [], [1], [~0], [10,25], [3,24,28], [1,7,16,31], [7,11,15,12,13],
  [2],[2,3],[2,3,7],[2,3,7,11],[2,3,7,11,15],
  [-2,-5,8],[-1,2,10,-5],
);
my @nsf = (
  [0], [8,16], [3,8,11], [4,5,25,30],[3,8,9,14,23],
  [15,-14,1],[-5,8,-13,5],
  [2,1],[2,3,5],[2,3,7,9],[2,3,7,11,13],
);
is_deeply( [map { is_sumfree_set($_) } @sf], [map { 1 } 0..$#sf], "is_sumfree_set with sumfree sets" );
is_deeply( [map { is_sumfree_set($_) } @nsf], [map { 0 } 0..$#nsf], "is_sumfree_set with non-sumfree sets" );

###### setcontains
is( setcontains([],[]), 1, "empty set contains empty set");
is( setcontains([1],[]), 1, "regular set contains empty set");
is( setcontains([],[1]), 0, "empty set does not contain regular set");

is( setcontains([1,3,5,8],[1,3,5]), 1, "setcontains basic true");
is( setcontains([1,3,5,8],[3,5,9]), 0, "setcontains basic false");

is( setcontains([1,8],[1,8,9]), 0, "setcontains with bigger subset");
is( setcontains([1..8],[-5..-1]), 0, "setcontains with smaller subset");
is( setcontains([1..8],[-5..-1]), 0, "setcontains with smaller subset");
is( setcontains([-5..5],[-7..-3]), 0, "setcontains with small bottom overlap");
is( setcontains([-5..5],[3..8]), 0, "setcontains with small top overlap");

is( setcontains([-5..5],[-3..3]), 1, "setcontains both signs subset true");
is( setcontains([-5..5],[-1..6]), 0, "setcontains both signs subset false");
is( setcontains([-5..-1],[-5,-3,-1]), 1, "setcontains neg true");
is( setcontains([-5..-1],[-6,-4,-2]), 0, "setcontains neg false");

my $ivpos = subint(powint(2,63),1000000);
my $ivneg = negint($ivpos);
my $uvpos = addint(powint(2,63),1000000);
my $uvneg = negint($uvpos);

# mix large IVs near min/max
is( setcontains([$ivneg,$ivpos],[$ivneg]), 1, "setcontains ivneg 1");
is( setcontains([$ivneg,$ivpos],[$ivpos]), 1, "setcontains ivneg 2");
is( setcontains([$ivpos],$ivneg), 0, "setcontains ivneg 3");
is( setcontains([$ivneg],$ivneg), 1, "setcontains ivneg 4");
is( setcontains([$ivpos],$ivpos), 1, "setcontains ivneg 5");
is( setcontains([$ivneg],$ivpos), 0, "setcontains ivneg 6");

# mix negative IV and positive UV
is( setcontains([$ivneg,$uvpos],[$ivneg]), 1, "setcontains ivneg 1");
is( setcontains([$ivneg,$uvpos],[$uvpos]), 1, "setcontains ivneg 2");
is( setcontains([$uvpos],$ivneg), 0, "setcontains ivneg 3");
is( setcontains([$ivneg],$ivneg), 1, "setcontains ivneg 4");
is( setcontains([$uvpos],$uvpos), 1, "setcontains ivneg 5");
is( setcontains([$ivneg],$uvpos), 0, "setcontains ivneg 6");

# mix negative and positive both near 64-bit
is( setcontains([$uvneg,$uvpos],[$uvneg]), 1, "setcontains uvneg 1");
is( setcontains([$uvneg,$uvpos],[$uvpos]), 1, "setcontains uvneg 2");
is( setcontains([$uvpos],$uvneg), 0, "setcontains uvneg 3");
is( setcontains([$uvneg],$uvneg), 1, "setcontains uvneg 4");
is( setcontains([$uvpos],$uvpos), 1, "setcontains uvneg 5");
is( setcontains([$uvneg],$uvpos), 0, "setcontains uvneg 6");
is( setcontains([$uvneg],-1000000), 0, "setcontains uvneg 7");

is( setcontains([1,3,5,8],$bi1), 0, "setcontains bigint false");
is( setcontains([$bi2,$bi3],$bi1), 0, "setcontains bigint false");
is( setcontains([$bi1,$bi3],$bi1), 1, "setcontains bigint true");
is( setcontains([$bi1,$bi3],[]), 1, "setcontains bigint empty set");
is( setcontains([$bi1,$bi3],[1]), 0, "setcontains bigint false");
is( setcontains([$bi1,$bi3],[$bi2]), 0, "setcontains bigint false");
is( setcontains([$bi1,$bi2],[$bi2]), 1, "setcontains bigint true");

# We think we work with unsorted second arg
is( setcontains([1..8],[5,4,5,1,3]), 1, "setcontains with non-set second arg");
