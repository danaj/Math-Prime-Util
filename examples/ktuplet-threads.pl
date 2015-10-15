#!/usr/bin/env perl
use strict;
use warnings;
use threads;
use threads::shared;

# TODO: Should have a pretty front end with option arguments etc.
# TODO: Should figure out the number of threads automatically.

use ntheory ":all";
use Math::BigInt;
$|=1;

my $nthreads = 8;
my $low   = Math::BigInt->new(1);
my $high  = Math::BigInt->new(10) ** 16;
my $range = Math::BigInt->new(10) ** 12;    # 10^13 or 10^14 for large clusters

#my @cl = (2,6,8,12,18,20,26,30,32);   # A027569
#my @cl = (6,12,16,18,22,28,30,36,40,42,46,48);   # A214947
#my @cl = (2,6,8,12,18,20,26,30,32,36,42,48,50);   # A257167
#my @cl = (2,6,8,12,18,20,26,30,32,36,42,48,50,56);   # A257304
#my @cl = (2,6,12,14,20,24,26,30,36,42,44,50,54,56,62,66); # 17 number 4
#my @cl = (6,10,12,16,22,24,30,34,36,40,42);  # A213601
#my @cl = (2,6,8,12,18,20,26,30,32,36,42);  # A213645
#my @cl = (6,12,16,18,22,28,30,36,40,42,46,48);   # A214947
#my @cl = (4,6,10,16,18,24,28,30,34,40,46,48);   # A257137
#my @cl = (4,6,10,16,18,24,28,30,34,36,46,48);   # A257138
my @cl = (2,6,8,18,20,30,32,36,38);   # Federighi

my $totresults = int( ($high+$range-1) / $range ) - 1;
#print "totresults: $totresults\n";

my @done :shared;
my @results :shared;
my $n :shared;

my @threads;
push @threads, threads->create('findtuple', $_) for 0 .. $nthreads-1;
$_->join() for (@threads);

while ($n <= $totresults && $done[$n]) {
  print "$results[$n]\n" if length($results[$n]) > 0;
  undef $results[$n];
  $n++;
}

sub findtuple {
  my $tnum = shift;
  my($res, $m, $tlow, $thigh);

  $m = $tnum;
  $tlow = $low + $m*$range;
  $n = 0 if $tnum == 0;

  while ($tlow <= $high) {
    $thigh = $tlow + $range - 1;
    $thigh = $high if $thigh > $high;
    if (scalar(@cl) > 9) {
      $res = join("\n", Math::Prime::Util::GMP::sieve_prime_cluster($tlow, $thigh, @cl));
    } else {
      $res = join("\n", sieve_prime_cluster($tlow, $thigh, @cl));
    }
    {
      lock(@done);
      $done[$m] = 1;
      $results[$m] = $res;
      if (1 && $tnum == 0) {
        while ($n <= $totresults && $done[$n]) {
          print "$results[$n]\n" if length($results[$n]) > 0;
          undef $results[$n];
          $n++;
        }
      }
    }
    $m += $nthreads;
    #$tlow = $low + $m*$range;
    $tlow += $nthreads * $range;
  }

  return 1;
}
