#!/usr/bin/env perl
use warnings;
use strict;
use Math::Prime::Util qw/:all/;
use List::Util qw/sum/;

my @tp;
my $p = 7;
while (1) {
  $p = next_prime($p);
  next unless $p =~ /^[2357]/ && $p =~ /[2357]$/;   # p ends are prime
  my $len = 1;
  while (++$len < length($p)) {
    last unless is_prime(substr($p, 0, $len)) && is_prime(substr($p, -$len));
  }
  next unless $len == length($p);
  push @tp, $p;
  last if scalar @tp >= 11;
}
print sum(@tp), "\n";
