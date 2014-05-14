#!/usr/bin/env perl
use warnings;
use strict;
use 5.14.0;
use Math::Prime::Util qw/:all/;
use Benchmark qw/:all/;

my $lim = shift || 1000;

# Michael B Porter proposed this OEIS sequence:
#
#   a(n) = m such that sigma(m) + sigma(m+1) + ... + sigma(m+n-1) is prime
#
# http://oeis.org/wiki/User:Michael_B._Porter
#
# Charles R Greathouse IV suggested this as an efficient computation:
#   a(n)=my(t=sum(i=1,n,sigma(i)),k=1);while(!isprime(t),t-=sigma(k)-sigma(n+k);k++);k
# which can be turned into a vector as:
#   vector(1000,i,a(i))
#
# Pari does this for 10k elements in ~15 seconds.
# Version opt2 does it in Perl in 3.0s.
# For 20k it's 63s in Pari, 12s in Perl.
# Of course Pari could be optimized as well.

sub simple {
  my $lim = shift;
  my @list;
  foreach my $n (1 .. $lim) {
    my($m, $sum) = (1, 0);
    while (!is_prime($sum)) {
      $sum = 0;
      $sum += divisor_sum($m+$_) for 0..$n-1;
      $m++;
    }
    push @list, $m-1;
  }
  return @list;
}
# perl -MMath::Prime::Util=:all -E 'my @list; foreach my $n (1 .. 1000) { my ($m,$sum) = (1,0); while (!is_prime($sum)) { $sum = 0; $sum += divisor_sum($m+$_) for 0..$n-1; $m++; } push @list, $m-1; } say join ",", @list;'

sub crg4 {
  my $lim = shift;
  my @list;
  foreach my $n (1 .. $lim) {
    my($k, $t) = (1,0);
    $t += divisor_sum($_) for 1..$n;
    while (!is_prime($t)) {
      $t -= divisor_sum($k)-divisor_sum($n+$k);
      $k++;
    }
    push @list,$k;
  }
  return @list;
}
# perl -MMath::Prime::Util=:all -E 'my @list; foreach my $n (1 .. 10000) { my($k,$t)=(1,0); $t += divisor_sum($_) for 1..$n; while (!is_prime($t)) { $t -= divisor_sum($k)-divisor_sum($n+$k); $k++; } push @list, $k; } say join ",", @list;'
# 9.8s for 10k

sub opt1 {
  my $lim = shift;
  my @list = map {
    my($n,$t,$k) = ($_,0,1);
    $t += divisor_sum($_) for 1..$n;
    while (!is_prime($t)) {
      $t -= divisor_sum($k) - divisor_sum($n+$k);
      $k++;
    }
    $k;
  } 1 .. $lim;
  return @list;
}
# perl -MMath::Prime::Util=:all -E 'say join ",", map { my($n,$t,$k) = ($_,0,1); $t += divisor_sum($_) for 1..$n; while (!is_prime($t)) { $t -= divisor_sum($k) - divisor_sum($n+$k); $k++; } $k; } 1 .. 10000'
# 9.5s for 10k

sub opt2 {
  my $lim = shift;
  my @ds;
  my @list = map {
    my($n,$t,$k) = ($_,0,1);
    $ds[$n] //= divisor_sum($n);
    $t += $ds[$_] for 1..$n;
    while (!is_prime($t)) {
      $ds[$n+$k] //= divisor_sum($n+$k);
      $t -= $ds[$k] - $ds[$n+$k];
      $k++;
    }
    $k;
  } 1 .. $lim;
  return @list;
}
# perl -MMath::Prime::Util=:all -E '@ds = (1,1); say join ",", map { my($n,$t,$k) = ($_,0,1); $t += $ds[$_] for 1..$n; while (!is_prime($t)) { $ds[$n+$k] //= divisor_sum($n+$k); $t -= $ds[$k] - $ds[$n+$k]; $k++; } $k; } 1..10000'
# 3.0s for 10k

# Verify
{
  my $vlim = 100;
  my @a1 = simple($vlim);
  my @a2 = crg4($vlim);
  my @a3 = opt1($vlim);
  my @a4 = opt2($vlim);
  foreach my $i (0 .. $vlim-1) {
    die "Mismatch in crg4 at $i" unless $a1[$i] == $a2[$i];
    die "Mismatch in opt1 at $i" unless $a1[$i] == $a3[$i];
    die "Mismatch in opt2 at $i" unless $a1[$i] == $a4[$i];
  }
}

cmpthese(-5, {
  #'simple' => sub { simple($lim) },
  'crg4'   => sub { crg4($lim) },
  'opt1'   => sub { opt1($lim) },
  'opt2'   => sub { opt2($lim) },
});

#say join ", ", opt1($lim);
