#!/usr/bin/env perl
use warnings;
use strict;
use Math::Prime::Util qw/:all/;

sub is_perm {
  my($a,$b) = @_;
  return length($a) == length($b) &&
         join("",sort split(//,$a)) eq join("",sort split(//,$b));
}

my ($maxn, $minratio) = (0, 1000000);
foreach my $n (2 .. 10_000_000) {
  my $totient = euler_phi($n);
  my $ratio = $n / $totient;
  ($maxn, $minratio) = ($n, $ratio) if $ratio < $minratio && is_perm($totient, $n);
}
print "$maxn  $minratio\n";
