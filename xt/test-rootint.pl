#!/usr/bin/env perl
use strict;
use warnings;
use Carp qw/carp croak confess/;
$| = 1;  # fast pipes

use Math::Prime::Util qw/rootint powint is_perfect_power is_power
                         addint subint/;
use Math::Prime::Util::PP;   # Have it available for comparison

my @rootmax = (0,0,4294967295,2642245,65535,7131,1625,565,255,138,84,56,40,30,23,19,15,13,11,10,9,8,7,6,6,5,5,5,4,4,4,4,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2);


for my $k (reverse 3 .. 63) {

  my($j,$nextroot,$r) = (0,1);
  # 1. Check all small inputs
  for my $n (0..200000) {
    my $flag = 0 + ($n == $nextroot);
    if ($flag) {
      $j++;
      $nextroot = powint($j+1,$k);
    } else {
      $flag |= ($n == 0);
    }
    #croak "is_perfect_power($n) != $flag" if is_perfect_power($n) != $flag;
    croak "is_power($n,$k) != $flag" unless
     ($flag && is_power($n,$k,\$r) && $r == $j) || (!$flag && !is_power($n,$k));
    croak "bad rootint($n,$k)" unless rootint($n,$k) == $j;
  }
  # 2. Check all j^k and a small window around them
  while (++$j <= $rootmax[$k]) {
    my $n = powint($j,$k);
    croak "is_power($n,$k) != 1" unless is_power($n,$k,\$r);
    croak "is_power($n,$k) root $r != $j" if $r != $j;
    for my $i (-10 .. 10) {
      my $rexpect = $j - ($i < 0);
      $r = rootint($n+$i,$k);
      croak "rootint($n + $i,$k) = $r != $rexpect" unless $r == $rexpect;
    }
  }
  # 3. Check values near 2^64
  my $rexpect = $j-1;
  for my $i (0 .. 10000) {
    my $n = subint(~0,$i);
    $r = rootint($n+$i,$k);
    croak "rootint($n,$k) = $r != $rexpect" unless $r == $rexpect;
  }
  print "$k ";

}
print "\n";
