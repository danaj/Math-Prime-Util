#!/usr/bin/env perl
use warnings;
use strict;
use v5.16;
use ntheory ":all";

my $small_nth = 1e7;
my $small_rnth = 1e6;
my $rp_inc = 1e9;

print "Verifying nth prime bounds up to $small_nth\n";
{
  my $n = 1;
  forprimes {
    my $p = $_;
    my($l,$u) = (nth_prime_lower($n),nth_prime_upper($n));
    die "$n: $l $p $u" unless $l <= $p && $u >= $p;
    $n++;
  } $small_nth;
}

print "Verifying nth Ramanujan prime bounds to $small_rnth\n";
{
  my $r = ramanujan_primes($small_rnth);
  for (0 .. $#$r) {
    my $n = $_+1;
    my $rn = $r->[$_];
    my($l,$u) = (nth_ramanujan_prime_lower($n),nth_ramanujan_prime_upper($n));
    die "$n: $l $rn $u" unless $l <= $rn && $u >= $rn;
  }
}

print "Verifying nth Ramanujan prime bounds:\n";
{
  my $s = 0;
  my $n = 1;
  while ($s < 1e12) {
    my $r = ramanujan_primes($s, $s + $rp_inc - 1);
    for (0 .. $#$r) {
      my $rn = $r->[$_];
      my($l,$u) = (nth_ramanujan_prime_lower($n),nth_ramanujan_prime_upper($n));
      die "$n: $l $rn $u" unless $l <= $rn && $u >= $rn;
      #die "$n: $rn" unless $rn == nth_ramanujan_prime($n);
      $n++;
    }
    print "   $s + $rp_inc\n";
    $s += $rp_inc;
  }
}
