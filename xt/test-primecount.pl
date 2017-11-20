#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util qw/prime_count primes prime_precalc urandomm/;
prime_precalc(1e8);
$|=1;

if (1) {
print "Testing prime counts for (a,b) for a in {3..15}, b in units of 1000\n";
for my $bm (0 .. 1000) {
  my $b = 1000 * $bm + urandomm(1000);
  print "$b " unless $bm % 10;
  for my $a (3 .. 15) {
    my($pc1,$pc2) = (prime_count($a,$b), scalar(@{primes($a,$b)}));
    die "($a,$b) =>  $pc1  !=  $pc2\n" unless $pc1 == $pc2;
  }
}
print "\n";
}

if (1) {
print "Testing prime counts for (a,b) for random a,b in 1..1e6\n";
for my $c (1 .. 10000) {
  print "$c " unless $c % 1000;
  my $b = urandomm(1e6);
  my $a = urandomm($b);
  my($pc1,$pc2) = (prime_count($a,$b), scalar(@{primes($a,$b)}));
  die "($a,$b) =>  $pc1  !=  $pc2\n" unless $pc1 == $pc2;
}
print "\n";
}

if (1) {
print "Testing prime counts for (a,b) for all b <= 1000, a <= b\n";
for my $b (0 .. 1000) {
  print "$b " unless $b % 100;
  for my $a (0 .. $b) {
    my($pc1,$pc2) = (prime_count($a,$b), scalar(@{primes($a,$b)}));
    die "($a,$b) =>  $pc1  !=  $pc2\n" unless $pc1 == $pc2;
  }
}
print "\n";
}
