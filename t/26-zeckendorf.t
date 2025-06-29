#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/tozeckendorf fromzeckendorf fromdigits todigits/;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};

my @z20 = (qw/0 1 10 100 101 1000 1001 1010 10000 10001 10010 10100 10101 100000 100001 100010 100100 100101 101000 101001 101010/);
my @toz = (
  [24, "1000100"],
  [27, "1001001"],
  [568, "1010010100000"],
  ["4294967295", "1010001000010101000101000100000001000100100100"],
  ["18446744073709551615", "10100101000100000101000100010010001001000000001001000100100010101000100000101000101000001010"],
  ["79228162514264337593543950335", "101010010101000101000100010010101010101010010100001010000100000001010101000010000101000000101000010010010001010010001000101000100000101001"],
);

plan tests => 0 + 1 + scalar(@toz) # tozeckendorf
                + 1 + scalar(@toz) # fromzeckendorf
                + 1;

###### tozeckendorf
is_deeply([map { tozeckendorf($_) } 0..20], \@z20, "tozeckendorf for 0..20");
for my $pair (@toz) {
  my($n,$s) = @$pair;
  is(tozeckendorf($n), $s, "tozeckendorf($n)");
}

###### fromzeckendorf
is_deeply([map { fromzeckendorf($_) } @z20], [0..20], "fromzeckendorf(Z(0..20))");
for my $pair (@toz) {
  my($n,$s) = @$pair;
  my $sp = $s;  $sp = substr($sp,0,30)."..." if length($s) > 33;
  is("".fromzeckendorf($s), $n, "fromzeckendorf($sp)");
}

######
is(fromdigits(tozeckendorf(24),2), 68, "fromdigits(tozeckendorf(24),2) = 68");

# TODO: Test for invalid inputs
