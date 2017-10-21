#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/ sqrtint rootint logint /;

#my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
#my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
#my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
#my $usegmp= Math::Prime::Util::prime_get_config->{'gmp'};
#$use64 = 0 if $use64 && 18446744073709550592 == ~0;

my @roots = (
  [25,  3, 15625],
  [13,  4, 28561],
  [13,  5, 371293],
  [25,  6, 244140625],
  [ 7,  7, 823543],
  [13,  8, 815730721],
  [ 7,  9, 40353607],
  [13, 10, "137858491849"],
  [21, 11, "350277500542221"],
  [25, 12, "59604644775390625"],
  [ 7, 13, "96889010407"],
  [ 7, 14, "678223072849"],
  [13, 16, "665416609183179841"],
  [13, 17, "8650415919381337933"],
  [ 7, 18, "1628413597910449"],
  [ 6, 19, "609359740010496"],
  [ 3, 21, "10460353203"],
  [ 3, 23, "94143178827"],
  [ 3, 25, "847288609443"],
  [ 3, 29, "68630377364883"],
);

plan tests => 0 + 4  # sqrtint
                + 6  # rootint
                + 5; # logint

###### sqrtint
is_deeply( [map { sqrtint($_) } 0..100], [map { int(sqrt($_)) } 0..100], "sqrtint 0 .. 100" );
is( sqrtint(1524155677489), 1234567, "sqrtint(1234567^2) = 1234567" );
is( sqrtint(1524158146623), 1234567, "sqrtint(1234568^2-1) = 1234567" );
is( sqrtint(1524155677488), 1234566, "sqrtint(1234567^2-1) = 1234566" );

###### rootint
# TODO: croak if n < 0 or k < 1
is(rootint(928342398,1), 928342398, "rootint(928342398,1) returns 928342398");
is(rootint(88875,3), 44, "rootint(88875,3) returns 44");
is(rootint("266667176579895999",3), 643659, "integer third root of 266667176579895999 is 643659");
{
  my(@got, @expected);
  for my $arr (@roots) {
    my($b, $k, $n) = @$arr;
    push @expected, [$b,$n];
    my $rk;
    my $r = rootint($n,$k,\$rk);
    push @got, [$r,$rk];
  }
  is_deeply( \@got, \@expected, "rootint on perfect powers where log fails" );
}
is( rootint("43091031920942300256108314560009772304748698124094750326895058640841523270081624169128280918534127523222564290447104831706207227117677890695945149868732770531628297914633063561406978145215542597509491443634033203125",23), 2147483645, "integer 23rd root of a large 23rd power" );
is( rootint("43091031920942300256108314560009772304748698124094750326895058640841523270081624169128280918534127523222564290447104831706207227117677890695945149868732770531628297914633063561406978145215542597509491443634033203124",23), 2147483644, "integer 23rd root of almost a large 23rd power" );

###### logint
is_deeply( [map { logint($_,2) } 1..200], [map { int(log($_)/log(2)+1e-10) } 1..200], "logint base 2: 0 .. 200" );
is_deeply( [map { logint($_,3) } 1..200], [map { int(log($_)/log(3)+1e-10) } 1..200], "logint base 3: 0 .. 200" );
is_deeply( [map { logint($_,5) } 1..200], [map { int(log($_)/log(5)+1e-10) } 1..200], "logint base 3: 0 .. 200" );
{
  my $be;
  is( logint(19284098234,16,\$be), 8, "logint(19284098234,16) = 8" );
  is( $be, 16**8, "power is 16^8" );
}
