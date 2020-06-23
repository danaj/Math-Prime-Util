#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_smooth is_rough factor vecnone/;

plan tests => 13+13 + 4+4;

for my $n (0 .. 12) {
  my @exp = map { fac_is_smooth($n, $_) } 0 .. 12;
  my @got = map {     is_smooth($n, $_) } 0 .. 12;
  is_deeply( \@got, \@exp, "is_smooth($n, 0..12)" );
}
for my $n (0 .. 12) {
  my @exp = map { fac_is_rough($n, $_) } 0 .. 12;
  my @got = map {     is_rough($n, $_) } 0 .. 12;
  is_deeply( \@got, \@exp, "is_rough($n, 0..12)" );
}

{
  my $n = "1377276413364943226363244108454842276965894752197358387200000"; # 97
  is( is_smooth($n, 23), 0, "large 97-smooth number" );
  is( is_smooth($n, 96), 0, "large 97-smooth number" );
  is( is_smooth($n, 97), 1, "large 97-smooth number" );
  is( is_smooth($n, 98), 1, "large 97-smooth number" );
}
{
  my $n = "172864518041328651521584134678230948270774322090771071422829"; # 2081
  is( is_smooth($n, 4073), 1, "large 4073-smooth, 2081-rough number" );
  is( is_rough($n, 2080), 1, "large 4073-smooth, 2081-rough number" );
  is( is_rough($n, 2081), 1, "large 4073-smooth, 2081-rough number" );
  is( is_rough($n, 2082), 0, "large 4073-smooth, 2081-rough number" );
}

sub fac_is_smooth {
  my($n, $k) = @_;
  # True if no prime factors of n are larger than k
  return 0+(vecnone { $_ > $k } factor($n));
}

sub fac_is_rough {
  my($n, $k) = @_;
  # True if no prime factors of n are smaller than k
  return 0+(vecnone { $_ < $k } factor($n));
}
