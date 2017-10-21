#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_square_free/;

my %isf = map { $_ => 0 } 0 .. 16,758096738,434420340,870589313,695486396,602721315,418431087;
$isf{$_} = 1 for 1,2,3,5,6,7,10,11,13,14,15,752518565,723570005,506916483,617459403;

plan tests => scalar(keys %isf);

while (my($n, $isf) = each (%isf)) {
  is( is_square_free($n), $isf, "is_square_free($n)" );
}
