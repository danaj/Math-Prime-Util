#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_powerful factor_exp vecall vecnone/;

plan tests => 3+1+10+4+2+1;

{
  my @exp = map { fac_is_powerful($_, 2) } 0 .. 1000;
  is_deeply( [map { is_powerful($_,2) } 0..1000], \@exp, "is_powerful(0..1000,2)");
  is_deeply( [map { is_powerful($_) } 0..1000], \@exp, "is_powerful(0..1000)");
  is_deeply( [map { is_powerful($_,0) } 0..1000], \@exp, "is_powerful(0..1000,0)");
}
{
  ok( vecall(sub { is_powerful($_,1) }, 0..1000), "is_powerful(n,1) = 1");
}

for my $k (3 .. 12) {
  my @exp = map { fac_is_powerful($_, $k) } 0 .. 1000;
  my @got = map {     is_powerful($_, $k) } 0 .. 1000;
  is_deeply(\@got, \@exp, "is_powerful(0..1000,$k)");
}

{
  my @pow2 = map { 5*5 * $_*$_ } 1..50;
  my @npow2 = map { 149 * $_*$_ } 1..50;
  my @pow3 = map { 7*7*7 * $_*$_*$_ } 1..50;
  my @npow3 = map { 4489 * $_*$_*$_ } 1..50;

  ok( vecall(sub{ is_powerful($_,2) }, @pow2), "small is_powerful(n,2), n powerful" );
  ok( vecall(sub{ is_powerful($_,3) }, @pow3), "small is_powerful(n,3), n powerful" );
  ok( vecnone(sub{ is_powerful($_,2) }, @npow2), "small is_powerful(n,2), n not powerful" );
  ok( vecnone(sub{ is_powerful($_,3) }, @npow3), "small is_powerful(n,3), n not powerful" );
}

is( is_powerful("1377276413364943226363244108454842276965894752197358387200000"), 0, "large easy non-powerful number" );
is( is_powerful("2346889178458529643625998598305409091755415961600000"), 1, "large easy powerful number" );

is( is_powerful("56648008573112538662596929676588737208124071038924666321487873929306609840197", 30), 0, "256-bit semiprime is not 30-powerful, without factoring" );

sub fac_is_powerful {
  my($n, $k) = @_;
  return 1 if $n <= 1 || $k <= 1;
  return (vecall { $_->[1] >= $k } factor_exp($n)) ? 1 : 0;
}
