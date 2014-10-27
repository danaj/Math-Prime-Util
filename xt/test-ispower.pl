#!/usr/bin/env perl
use strict;
use warnings;
$| = 1;  # fast pipes

use Math::Prime::Util qw/is_prime is_power/;
use Math::GMPz;

foreach my $e (5 .. 37) {
  next unless is_prime($e);
  print "$e ";
  for (3..1000000) {
    my $n = Math::GMPz->new($_) ** $e;
    last if $n > ~0;
    die "$n\n" unless is_power($n);
    foreach my $o (-10..10) {
      my $m = $n+$o;
      next if $m==$n;
      die "$m\n" if is_power($m) && int(sqrt($m))**2 != $m && $m!=2197;
    }
  }
}
print "\n";

my $int = 100000;
foreach my $i (1 .. 80*$int) {
  print "." unless $i % $int;
  my @iroots = (0,0,map { int($i ** (1.0/$_) + 0.00001) ** $_ } 2 .. 12);
  my $r;
  foreach my $e (2 .. 12) {
    if (is_power($i,$e,\$r)) { die "$i $e" unless $iroots[$e] == $i && $r ** $e == $i; }
    else                     { die "$i $e" unless $iroots[$e] != $i; }
  }
}
print "\n";
