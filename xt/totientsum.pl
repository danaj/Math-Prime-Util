#!/usr/bin/env perl
use strict;
use warnings;
$| = 1;  # fast pipes

use 5.020;
use experimental qw(signatures);
use Math::Prime::Util qw/euler_phi moebius divint sqrtint rootint vecsum sumtotient/;

my $n = shift || 1_000_000;
my $s;


# 1  sieve them all and sum
#print "sum totient $n:  ", vecsum(euler_phi(0,$n)), "\n";


# 2  sieve them all and sum
if (0) {
  $s = 0;
  $s += $_ for euler_phi(0,$n);
  print "sum2 $n:         $s\n";
}

# 3  sieve mobius and sum
if (0) {
  my @m = moebius(0,$n);
  $s = 1;
  for (1..$n) {
    if ($m[$_] != 0) {
      my $d = divint($n,$_);
      $s += $m[$_] * $d * $d;
    }
  }
  $s >>= 1;
  print "sum3 $n:         $s\n";
}

# 4  from sidef
my @tsum_cache = (0);

sub R {
  my $n = shift;

  my $s = sqrtint($n);

  if ($#tsum_cache < $s) {
    $tsum_cache[$_] = $tsum_cache[$_-1] + euler_phi($_)  for 1 .. $s;
  }
  return $tsum_cache[$n] if $n <= $#tsum_cache;

  my $L = ($n*($n-1)) >> 1;
  for my $k (2 .. divint($n, $s+1)) {
    $L -= R(divint($n, $k));
  }
  for my $k (1..$s) {
    $L -= $tsum_cache[$k] * (divint($n,$k) - divint($n,$k+1));
  }
  return $L;
}

#print "sum4 $n:         ", R($n),"\n";

# 5  from Andy on stackoverflow
my @x2;
sub R2 {
  my $n = shift;
  return 0 if $n <= 1;
  return $x2[$n] if defined $x2[$n];
  my $sum = ($n*($n-1)) >> 1;
  my $m = 2;
  while (1) {
    my $x = divint($n,$m);
    my $nxt = divint($n,$x);
    if ($nxt >= $n) {
      $x2[$n] = $sum - ($n-$m+1) * R2($x);
      return $x2[$n];
    }
    $sum -= ($nxt-$m+1) * R2($x);
    $m = $nxt + 1;
  }
}
#print "sum5 $n:         ", 1+R2($n),"\n";


# 6  Basic
sub R3 {
  my $n = shift;
  return $n if $n <= 2;
  my $sum = ($n*($n-1)) >> 1;
  for my $m (2 .. $n) {
    $sum -= R3(divint($n,$m));
  }
  $sum;
}
#print "sum6 $n:         ", R3($n),"\n";

# 7  from other definition
my @r4c;
sub R4 {
  my $n = shift;
  return $n if $n <= 2;
  return $r4c[$n] if defined $r4c[$n];
  my $sum = ($n*($n-1)) >> 1;
  my $sqrtn = sqrtint($n);
  for my $i (2 .. $sqrtn) {
    $sum -= R4(divint($n,$i));
  }
  for my $j (1 .. $sqrtn) {
    $sum -= R4($j) * (divint($n,$j) - divint($n,$j+1));
  }
  $r4c[$n] = $sum;
  $sum;
}
#print "sum7 $n:         ", R3($n),"\n";

# 8  from oeis 002088

my @acache;
my $_l = 0;
sub a002088 {
  my $n = shift;
  return $n if $n <= 2;
  return $acache[$n] if defined $acache[$n];
  my $c = 0;
  my $j = 2;
  my $k1 = divint($n, $j);
  while ($k1 > 1) {
    my $j2 = divint($n, $k1) + 1;
    $c += ($j2 - $j) * (2*a002088($k1)-1);
    $j = $j2;
    $k1 = divint($n, $j2);
  }
  my $sum = (($n * ($n-1) - $c + $j) >> 1);
  $acache[$n] = $sum;
  return $sum;
}
#print "sum8 $n:         ", a002088($n),"\n";

# https://raw.githubusercontent.com/trizen/perl-scripts/master/Math/partial_sums_of_euler_totient_function_fast.pl
sub partial_sums_of_euler_totient($n) {
    my $s = sqrtint($n);
    my @euler_sum_lookup = (0);
    my $lookup_size = 2 * rootint($n, 3)**2;
    my @euler_phi   = euler_phi(0, $lookup_size);
    foreach my $i (1 .. $lookup_size) {
        $euler_sum_lookup[$i] = $euler_sum_lookup[$i - 1] + $euler_phi[$i];
    }
    my %seen;
    sub ($n) {
        if ($n <= $lookup_size) {
            return $euler_sum_lookup[$n];
        }
        if (exists $seen{$n}) {
            return $seen{$n};
        }
        my $s = sqrtint($n);
        my $T = ($n * ($n + 1)) >> 1;
        foreach my $k (2 .. int($n / ($s + 1))) {
            $T -= __SUB__->(int($n / $k));
        }
        foreach my $k (1 .. $s) {
            $T -= (int($n / $k) - int($n / ($k + 1))) * __SUB__->($k);
        }
        $seen{$n} = $T;
    }->($n);
}
print "sum9 $n:         ", partial_sums_of_euler_totient($n),"\n";

#  5,8 are fast.  9 is 10x faster than those.

# MPU uses Trizen's algorithm from 9.
#
# MPU PP is a little faster than 9 and doesn't overflow with n > 32-bit.
# MPU XS is about 100x faster, and has 128-bit support.

print "MPU: $n:         ", sumtotient($n),"\n";
