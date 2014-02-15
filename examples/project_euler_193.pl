#!/usr/bin/env perl
use warnings;
use strict;
use Math::Prime::Util qw/moebius mertens/;

my $lim = shift || 2**50;
my $method = shift || 'mertens';

# See http://arxiv.org/pdf/1107.4890v1.pdf

#  2.9s  mertens
#  9.9s  block     
# 10.0s  monolithic
# 36.0s  simple
# lots   brute

my $sum = 0;

if ($method eq 'brute') {
  # Far too slow
  for (1 .. $lim) { $sum++ if moebius($_) }
} elsif ($method eq 'simple') {

  # Basic application of theorem 1.
  for (1 .. int(sqrt($lim)+0.001)) {
    $sum += moebius($_) * int($lim/($_*$_));
  }

} elsif ($method eq 'monolithic') {

  # Efficient theorem 1, but lots of memory.
  my @mob = moebius(0, int(sqrt($lim)+0.001));
  for (1 .. $#mob) { $sum += $mob[$_] * int($lim/($_*$_)) if $mob[$_]; }

} elsif ($method eq 'block') {

  # Break up into chunks to constrain memory.
  my($beg,$end,$mlim) = (1, 1, int(sqrt($lim)+0.001));
  while ($beg < $mlim) {
    $end = $beg + 100_000 - 1;
    $end = $mlim if $end > $mlim;
    my @mob = moebius($beg,$end);
    for ($beg .. $end) {
      $sum += $mob[$_-$beg] * int($lim/($_*$_)) if $mob[$_-$beg];
    }
    $beg = $end+1;
  }
} elsif ($method eq 'mertens') {

  # Pawlewicz's method, using chunked S1, and no optimization for Mertens.

  my $I = 50;   # Tune as desired.
  my $D = int(sqrt($lim/$I)+0.00001);
  my ($S1, $S2) = (0,0);

  # S1
  my $chunk = 100_000;
  for (my $beg = 1; $beg < $D; $beg += $chunk) {
    my $end = $beg + $chunk - 1;
    $end = $D if $end > $D;
    my @mob = moebius($beg,$end);
    for ($beg .. $end) {
      $S1 += $mob[$_-$beg] * int($lim/($_*$_)) if $mob[$_-$beg];
    }
  }
  # S2
  for (1 .. $I-1) {
    my $xi = int(sqrt($lim/$_)+0.00001);
    $S2 += mertens($xi);
  }
  $S2 -= ($I-1) * mertens($D);

  $sum = $S1 + $S2;
}

print "$sum\n";
