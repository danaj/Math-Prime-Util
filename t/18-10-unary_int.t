#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/absint negint signint/;

my @pos = (qw/32767 32768
               1073741823 1073741824 1073741825
               2147483647 2147483648 2147483649
               4294967295 4294967296
               8589934591 8589934592
               9223372036854775807 9223372036854775808
               18446744073709551614 18446744073709551615
               18446744073709551616 18446744073709551617
               36893488147419103231 36893488147419103232
               827639478068904540012
               829342632809347434459
               1726145541361106236340
               2555488174170453670799
               4281633715531559907139
               1178630961471601951655862
               230948092384903284908329048239084023984092384
               982349082340982348502392937523840234029384908325098234/);
my @neg = map { '-' . $_ } @pos;

plan tests => 3 + 2 + 2 + 2 + 2;

{
  my(@absgot,@neggot,@sgngot, @absexp,@negexp,@sgnexp);
  for my $n (-100 .. 100) {
    push @absgot, absint($n);
    push @neggot, negint($n);
    push @sgngot, signint($n);
    push @absexp, abs($n);
    push @negexp, ($n == 0) ? 0 : -$n;
    push @sgnexp, 0 + ($n > 0) - ($n < 0);
  }
  is_deeply( \@absgot, \@absexp, "absint(-100..100)" );
  is_deeply( \@neggot, \@negexp, "negint(-100..100)" );
  is_deeply( \@sgngot, \@sgnexp, "signint(-100..100)" );
}

is_deeply( [absint(0), negint(0), signint(0)], [0,0,0], "absint(0), negint(0), signint(0)" );
is_deeply( [absint(-0), negint(-0), signint(-0)], [0,0,0], "absint(-0), negint(-0), signint(-0)" );

is_deeply([map { absint($_) } @pos], \@pos, "absint with positive inputs");
is_deeply([map { absint($_) } @neg], \@pos, "absint with negative inputs");

is_deeply([map { negint($_) } @pos], \@neg, "negint with positive inputs");
is_deeply([map { negint($_) } @neg], \@pos, "negint with negative inputs");

is_deeply([map {signint($_)} @pos], [map { 1} @pos], "signint with positive inputs");
is_deeply([map {signint($_)} @neg], [map {-1} @neg], "signint with negative inputs");
