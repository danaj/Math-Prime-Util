#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_powerful powerful_count nth_powerful
                         factor_exp vecall vecnone/;

my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};

plan tests => 3+1+10+4+2+1
            + 5 + 2*$extra
            + 3;

{
  my @exp = map { fac_is_powerful($_, 2) } 0 .. 258;
  is_deeply( [map { is_powerful($_,2) } 0..258], \@exp, "is_powerful(0..258,2)");
  is_deeply( [map { is_powerful($_) } 0..258], \@exp, "is_powerful(0..258)");
  is_deeply( [map { is_powerful($_,0) } 0..258], \@exp, "is_powerful(0..258,0)");
}

ok( vecall(sub { is_powerful($_,1) }, 0..32), "is_powerful(n,1) = 1");

for my $k (3 .. 12) {
  my @nums = (227411960,105218838,79368063,58308379,210322300,44982156,67831696,165946352,243118692,128757041,150085583);
  my @exp = map { fac_is_powerful($_, $k) } 0 .. 32, @nums;
  my @got = map {     is_powerful($_, $k) } 0 .. 32, @nums;
  is_deeply(\@got, \@exp, "is_powerful(n,$k) for 0..32 and 11 larger nums");
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


###### powerful_count
is_deeply( [map { powerful_count($_,1) } 1..20], [1..20], "powerful_count(n,1)=n" );

is_deeply( [map { powerful_count($_) } 1..20],
           [1,1,1,2,2,2,2,3,4,4,4,4,4,4,4,5,5,5,5,5],
           "powerful_count(1..20)" );
is_deeply( [map { powerful_count($_,3) } 1..20],
           [1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3],
           "powerful_count(1..20,3)" );
{
  my $maxk = ($extra) ? 30 : 15;
  my @pow14=(0,14,100,432, 2048, 6561, 16384, 59049, 131072, 262144, 531441, 1594323, 4194304, 8388608, 16777216, 43046721, 129140163, 268435456, 536870912, 1162261467, 3486784401, 8589934592, 17179869184, 34359738368, 68719476736, 137438953472, 274877906944, 549755813888, 1099511627776, 2199023255552, 4398046511104);
  is_deeply( [map { powerful_count($pow14[$_],$_) } 1..$maxk],
             [map { 14 } 1..$maxk],
             "powerful_count(x,1..$maxk) = 14" );
  is_deeply( [map { powerful_count($pow14[$_]-1,$_) } 1..$maxk],
             [map { 13 } 1..$maxk],
             "powerful_count(x-1,1..$maxk) = 13" );
}

# Rather slow in PP
if ($extra) {
  my @exp = (4,14,54,185,619,2027,6553,21044,67231,214122,680330,2158391,6840384,21663503);
  my @got = map { powerful_count(10**$_) } 1..14;
  is_deeply(\@got, \@exp, "2-powerful_count 10^1, 10^2, ..., 10^14");
}
if ($extra) {
  my @exp = (1, 1, 4, 10, 16, 26, 46, 77, 129, 204, 318, 495, 761, 1172, 1799, 2740, 4128, 6200, 9224, 13671, 20205, 29764);
  my @got = map { powerful_count("1".("0"x$_),7) } 1..22;
  is_deeply(\@got, \@exp, "7-powerful_count 10^1, 10^2, ..., 10^22");
}

###### nth_powerful

is(nth_powerful(100),3136,"3136 is the 100th powerful number");
SKIP: {
  skip "Skipping nth_powerful for k>2 in PP", 2 unless $usexs;
  is(nth_powerful(100,6),43046721,"43046721 is the 100th 6-powerful number");
  is(nth_powerful(12,15),16777216,"16777216 is the 12th 15-powerful number");
}


sub fac_is_powerful {
  my($n, $k) = @_;
  return 1 if $n <= 1 || $k <= 1;
  return 0 if $n < (1<<$k);
  return 0 if (!($n%2)) && ($n%4);
  return (vecall { $_->[1] >= $k } factor_exp($n)) ? 1 : 0;
}
