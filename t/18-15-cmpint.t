#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/cmpint/;

plan tests => 7 + 1;

###### cmpint
is(cmpint(1,2),-1,"1 < 2");
is(cmpint(2,1), 1,"2 > 1");
is(cmpint(2,2), 0,"2 == 2");
is(cmpint("18446744073709553664","18446744073709551615"),1,"2^64+2048 > 2^64-1");
is(cmpint("18446744073709551664","18446744073709551615"),1,"2^64+1048 > 2^64-1");
is(cmpint("18446744073709551615","18446744073709551616"),-1,"2^64-1 < 2^64");
is(cmpint("-18446744073709551615","18446744073709551615"),-1,"-2^64-1 < 2^64-1");

my @arr = (-36, 1, 17, 19, 2400);
is_deeply( [map { mediancmp($_,\@arr) } -50,-36,0,1,2,17,18,19,20,2400,3000],
           [-5,-4,-3,-2,-1,0,1,2,3,4,5],
           "Use cmpint as part of array median comparator" );

sub mediancmp {
  my($val,$ar) = @_;
  #return vecsum(map { cmpint($val,$_) } @$ar);
  my $sum = 0;
  $sum += cmpint($val,$_) for @$ar;
  $sum;
}
