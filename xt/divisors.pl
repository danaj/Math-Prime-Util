#!/usr/bin/env perl
use strict;
use warnings;
$| = 1;  # fast pipes

use Math::Prime::Util qw/divisors vecequal/;
use Math::Prime::Util::PP;
use Math::Prime::Util::GMP;

my($lim,$klim) = (200000,800);

print "scalar divisors from 0 to $lim ... ";
for my $n (0..$lim) {
  my $t1 = 0+Math::Prime::Util::divisors($n);
  my $t2 = 0+Math::Prime::Util::GMP::divisors($n);
  my $t3 = 0+Math::Prime::Util::PP::divisors($n);
  die "FAIL  scalar divisors $n:  $t1 $t2 $t3\n"
    unless $t1 == $t2 && $t2 == $t3;
}
print "PASS\n";

print "list divisors from 0 to $lim ... ";
for my $n (0..$lim) {
  my @t1 = Math::Prime::Util::divisors($n);
  my @t2 = Math::Prime::Util::GMP::divisors($n);
  my @t3 = Math::Prime::Util::PP::divisors($n);
  die "FAIL  list divisors $n:  [@t1] [@t2] [@t3]\n"
    unless vecequal(\@t1,\@t2) && vecequal(\@t1,\@t3);
}
print "PASS\n";


print "scalar restricted divisors from 0 to $klim ... ";
for my $n (0..$klim) {
  for my $k (0 .. $n+6) {
    my $t1 = 0+Math::Prime::Util::divisors($n,$k);
    my $t2 = 0+Math::Prime::Util::GMP::divisors($n,$k);
    my $t3 = 0+Math::Prime::Util::PP::divisors($n,$k);
    die "FAIL  scalar divisors $n,$k:  $t1 $t2 $t3\n"
      unless $t1 == $t2 && $t2 == $t3;
  }
}
print "PASS\n";

print "list restricted divisors from 0 to $klim .. ";
for my $n (0..$klim) {
  for my $k (0 .. $n+6) {
    my @t1 = Math::Prime::Util::divisors($n,$k);
    my @t2 = Math::Prime::Util::GMP::divisors($n,$k);
    my @t3 = Math::Prime::Util::PP::divisors($n,$k);
    die "FAIL  list divisors $n,$k:  [@t1] [@t2] [@t3]\n"
      unless vecequal(\@t1,\@t2) && vecequal(\@t1,\@t3);
  }
}
print "PASS\n";

print "PASS all comparisons for divisors() code\n";
