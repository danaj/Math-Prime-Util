#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_powerfree powerfree_count powerfree_sum
                         is_square_free vecsum vecmax factor_exp/;

my @simple = (0 .. 16,
              758096738,434420340,870589313,695486396,602721315,418431087,
              752518565,723570005,506916483,617459403);

plan tests => 3     # simple is square free
            + 11*2  # powerfree_count, powerfree_sum
            + 6+2;  # ""

is_deeply( [map { is_powerfree($_)   } @simple],
           [map { is_square_free($_) } @simple],
           "is_powerfree(n) matches is_square_free(n)" );
is_deeply( [map { is_powerfree($_)   } @simple],
           [map { ipf($_)            } @simple],
           "is_powerfree(n) works for simple inputs" );
is_deeply( [map { is_powerfree($_,3) } @simple],
           [map { ipf($_,3)          } @simple],
           "is_powerfree(n,3) works for simple inputs" );

foreach my $k (0..10) {
    my $n = 100;
    is_deeply(
        [map { powerfree_count($_, $k) } 0..$n],
        [map { scalar grep { is_powerfree($_, $k) } 1..$_ } 0..$n],
        "powerfree_count(0..$n, $k)"
    );
    is_deeply(
        [map { powerfree_sum($_, $k) } 0..$n],
        [map { vecsum(grep { is_powerfree($_, $k) } 1..$_) } 0..$n],
        "powerfree_sum(0..$n, $k)"
    );
}

is( powerfree_count(12345,2),  7503, "powerfree_count(12345,2) = 7503");
is( powerfree_count(12345,3), 10272, "powerfree_count(12345,3) = 10272");
is( powerfree_count(12345,4), 11408, "powerfree_count(12345,4) = 11408");
is( powerfree_sum(12345,2), 46286859, "powerfree_sum(12345,2) = 46286859");
is( powerfree_sum(12345,3), 63404053, "powerfree_sum(12345,3) = 63404053");
is( powerfree_sum(12345,4), 70415676, "powerfree_sum(12345,4) = 70415676");

is( powerfree_count(123456,32), 123456, "powerfree_count(123456,32) = 123456");
is( powerfree_sum(123456,32), 7620753696, "powerfree_sum(123456,32) = 7620753696");

sub ipf {
  my($n,$k) = @_;
  $k = 2 unless defined $k;
  $n = -$n if $n < 0;
  return 0 if $n == 0;
  return 1 if $n == 1;

  (vecmax(map { $_->[1] } factor_exp($n)) < $k) ? 1 : 0;
}
