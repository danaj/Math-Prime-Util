#!/usr/bin/env perl
use warnings;
use strict;
use Math::Prime::Util qw/is_prime primes/;

sub is_perm {
  my($a,$b) = @_;
  return length($a) == length($b) &&
         join("",sort split(//,$a)) eq join("",sort split(//,$b));
}

foreach my $inc2 (1 .. 1700) {
  my $inc = $inc2 * 2;
  foreach my $p (@{primes(1000,9999)}) {
    my($p2, $p3) = ($p+$inc, $p+$inc+$inc);
    last if $p3 > 9999;
    next unless is_prime($p2) && is_prime($p3);
    next unless is_perm($p, $p2) && is_perm($p, $p3);
    print "$p/$inc: $p $p2 $p3\n";
  }
}
