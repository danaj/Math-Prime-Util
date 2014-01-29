#!/usr/bin/env perl
use warnings;
use strict;
use Math::Prime::Util qw/:all/;

# Fill in the chains
my @achain = ( [0] );
foreach my $n (0 .. 50_000) {
  next if defined $achain[$n];
  my @seq = aliquot_sequence($n, 1_000_000);
  #print "chain for $n = ", join(",", @seq), "\n";
  while (@seq) {
    my $s = shift @seq;
    $achain[$s] = [$s, @seq] if !defined $achain[$s];
  }
}

# Find max chain length
my ($maxlen, $maxn) = (0, 0);
foreach my $n (0 .. 1_000_000) {
  next unless defined $achain[$n];
  next unless $achain[$n]->[0] == $achain[$n]->[-1];
  my $len = scalar @{$achain[$n]} - 1;
  ($maxlen, $maxn) = ($len, $n) if $len > $maxlen;
}
print "Max length: $maxlen.  n = $maxn\n";
print "Chain for $maxn: ", join(",", @{$achain[$maxn]}), "\n";

sub aliquot_sequence {
  my ($n, $max) = @_;
  my %hash;
  undef $hash{$n};
  my @seq = ($n);
  foreach my $len (1 .. 1000) {
    $n = divisor_sum($n)-$n;
    # Stop if we have exceeded the threshold
    last if $n > $max;
    # If we know how this chain ends, return it now
    return @seq, @{$achain[$n]} if defined $achain[$n];
    push @seq, $n;
    return @seq if exists $hash{$n} || $n == 0;
    undef $hash{$n};
  }
  return ();
}
