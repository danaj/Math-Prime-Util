#!/usr/bin/env perl
use warnings;
use strict;
use ntheory ":all";
use Math::GMP;

for my $n (8539783 .. 8539788) {
  for (1..2000) {
    my($u,$v) = lucas_sequence($n, 1, -1, $_);
    my $u1 = Math::GMP->new(lucasu(1,-1,$_)) % $n;
    my $v1 = Math::GMP->new(lucasv(1,-1,$_)) % $n;
    die "U $_ $n $u $u1" unless $u==$u1;
    die "V $_ $n $v $v1" unless $v==$v1;
  }
}
