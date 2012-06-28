#!/usr/bin/env perl
use strict;
use warnings;

use Benchmark qw/:all/;
use Devel::Size qw/total_size/;
use Math::Prime::Util;
use Math::Prime::FastSieve;
*mpu_erat = \&Math::Prime::Util::erat_primes;
*fs_erat = \&Math::Prime::FastSieve::primes;

my $upper = shift || 8192;
my $count = shift || -1;
my $countarg;
my $sum;

# This is like counting, but we want an array returned.
# The subs will compute a sum on the results.

# Times for 100k.
# Vs. MPU sieve, as we move from 8k to 10M:
#   Atkin MPTA, Rosetta 3 & 1, Shootout, Scriptol, DO Array, DJ Array, and
#   InMany all slow down.  Atkin 2 speeds up (from 65x slower to 54x slower).
# The DJ string methods have almost no relative slowdown, so stretch out their
# advantage over the others (In Many, DJ Array, DJ Vec, and DO Array).
my $pc_subs = {
  "Rosetta 4" => sub {$sum=0; $sum+=$_ for rosetta4($countarg);$sum;},  #    9/s
  "Atkin MPTA"=> sub {$sum=0; $sum+=$_ for atkin($countarg);$sum;},     #   11/s
  "Merlyn"    => sub {$sum=0; $sum+=$_ for merlyn($countarg);$sum;},    #   15/s
  "Rosetta 2" => sub {$sum=0; $sum+=$_ for rosetta2($countarg);$sum; }, #   16/s
  "DO Vec"    => sub {$sum=0; $sum+=$_ for daos_vec($countarg);$sum;},  #   16/s
  "Atkin 2"   => sub {$sum=0; $sum+=$_ for atkin2($countarg);$sum; },   #   17/s
  "Rosetta 3" => sub {$sum=0; $sum+=$_ for rosetta3($countarg);$sum; }, #   23/s
  "Rosetta 1" => sub {$sum=0; $sum+=$_ for rosetta1($countarg);$sum; }, #   26/s
  "Shootout"  => sub {$sum=0; $sum+=$_ for shootout($countarg);$sum; }, #   30/s
  "Scriptol"  => sub {$sum=0; $sum+=$_ for scriptol($countarg);$sum; }, #   33/s
  "DJ Vec"    => sub {$sum=0; $sum+=$_ for dj1($countarg);$sum; },      #   34/s
  "DO Array"  => sub {$sum=0; $sum+=$_ for daos_array($countarg);$sum;},#   41/s
  "DJ Array"  => sub {$sum=0; $sum+=$_ for dj2($countarg);$sum; },      #   63/s
  "In Many"   => sub {$sum=0; $sum+=$_ for inmany($countarg);$sum; },   #   86/s
  "DJ String1"=> sub {$sum=0; $sum+=$_ for dj3($countarg);$sum; },      #   99/s
  "DJ String2"=> sub {$sum=0; $sum+=$_ for dj4($countarg);$sum; },      #  134/s
  "MPFS Sieve"=> sub {                                                  # 1216/s
               $sum=0; $sum+=$_ for @{fs_erat($countarg)};;$sum; },
  "MPU Sieve" => sub {                                                  # 1290/s
               $sum=0; $sum+=$_ for @{mpu_erat(2,$countarg)};;$sum; },
};

my %verify = (
      10 => 17,
      11 => 28,
     100 => 1060,
     112 => 1480,
     113 => 1593,
     114 => 1593,
    1000 => 76127,
   10000 => 5736396,
  100000 => 454396537,
);

# Verify
while (my($name, $sub) = each (%$pc_subs)) {
  while (my($n, $v_pi_sum) = each (%verify)) {
    $countarg = $n;
    my $pi_sum = $sub->();
    die "$name ($n) = $pi_sum, should be $v_pi_sum" unless $pi_sum == $v_pi_sum;
  }
}
print "Done with verification, starting benchmark\n";

$countarg = $upper;
cmpthese($count, $pc_subs);



# www.scriptol.com/programming/sieve.php
sub scriptol {
  my($max) = @_;
  return 0 if $max < 2;
  return 1 if $max < 3;

  my @flags = (0 .. $max);
  for my $i (2 .. int(sqrt($max)) + 1)
  {
       next unless defined $flags[$i];
       for (my $k=$i+$i; $k <= $max; $k+=$i)
       {
           undef $flags[$k];
       }
  }
  return grep { defined $flags[$_] } 2 .. $max;
}

# http://dada.perl.it/shootout/sieve.perl.html
sub shootout {
  my($max) = @_;
  return 0 if $max < 2;
  return 1 if $max < 3;

  my @primes;
  my @flags = (0 .. $max);
  for my $i (2 .. $max) {
       next unless defined $flags[$i];
       for (my $k=$i+$i; $k <= $max; $k+=$i) {
           undef $flags[$k];
       }
       push @primes, $i;
  }
  @primes;
}

# http://c2.com/cgi/wiki?SieveOfEratosthenesInManyProgrammingLanguages
sub inmany {
  my($max) = @_;
  return 0 if $max < 2;
  return 1 if $max < 3;

  my @c;
  for(my $t=3; $t*$t<=$max; $t+=2) {
     if (!$c[$t]) {
         for(my $s=$t*$t; $s<=$max; $s+=$t*2) { $c[$s]++ }
     }
 }
 my @primes = (2);
 for(my $t=3; $t<=$max; $t+=2) {
     $c[$t] || push @primes, $t;
 }
 @primes;
 # grep { $c[$_] } 3 .. $max;
}

# http://rosettacode.org/wiki/Sieve_of_Eratosthenes#Perl
sub rosetta1 {
  my($max) = @_;
  return 0 if $max < 2;
  return 1 if $max < 3;

  my @primes;
  my @tested = (1);
  my $j = 1;
  while ($j < $max) {
    next if $tested[$j++];
    push @primes, $j;
    for (my $k= $j; $k <= $max; $k+=$j) {
      $tested[$k-1]= 1;
    }
  }
  @primes;
}

# http://rosettacode.org/wiki/Sieve_of_Eratosthenes#Perl
sub rosetta2 {
  my($max) = @_;
  return 0 if $max < 2;
  return 1 if $max < 3;

  my @primes;
  my $nonPrimes = '';
  foreach my $p (2 .. $max) {
    unless (vec($nonPrimes, $p, 1)) {
      for (my $i = $p * $p; $i <= $max; $i += $p) {
        vec($nonPrimes, $i, 1) = 1;
      }
      push @primes, $p;
    }
  }
  @primes;
}

# http://rosettacode.org/wiki/Sieve_of_Eratosthenes#Perl
sub rosetta3 {
  my($max) = @_;
  return 0 if $max < 2;
  return 1 if $max < 3;

  my(@s, $i);
  grep { not $s[ $i  = $_ ] and do
           { $s[ $i += $_ ]++ while $i <= $max; 1 }
       } 2 .. $max;
}

# http://rosettacode.org/wiki/Sieve_of_Eratosthenes#Perl
sub rosetta4 {
  my($max) = @_;
  return 0 if $max < 2;
  return 1 if $max < 3;

  my $i;
  my $s = '';
  grep { not vec $s, $i  = $_, 1 and do
          { (vec $s, $i += $_, 1) = 1 while $i <= $max; 1 }
       } 2 .. $max;
}

# From Math::Primes::TiedArray
sub atkin {
  my($max) = @_;
  return 0 if $max < 2;
  return 1 if $max < 3;
  return 2 if $max < 5;

    my $sqrt     = sqrt($max);
    my %sieve;
    foreach my $x ( 1 .. $sqrt ) {

        foreach my $y ( 1 .. $sqrt ) {

            my $n = 3 * $x**2 - $y**2;
            if (    $x > $y
                and $n <= $max
                and $n % 12 == 11 )
            {
                $sieve{$n} = not $sieve{$n};
            }

            $n = 3 * $x**2 + $y**2;
            if ( $n <= $max and $n % 12 == 7 ) {
                $sieve{$n} = not $sieve{$n};
            }

            $n = 4 * $x**2 + $y**2;
            if (    $n <= $max
                and ( $n % 12 == 1 or $n % 12 == 5 ) )
            {
                $sieve{$n} = not $sieve{$n};
            }
        }
    }
    # eliminate composites by sieving
    foreach my $n ( 5 .. $sqrt ) {

        next unless $sieve{$n};

        my $k = int(1/$n**2) * $n**2;
        while ( $k <= $max ) {
            $sieve{$k} = 0;
            $k += $n**2;
        }
    }
    my @primes = (2, 3);
    push @primes, grep { $sieve{$_} } 5 .. $max;
    @primes;
}

# Naive Sieve of Atkin, basically straight from Wikipedia.
#
# <rant>
#
# First thing to note about SoA, is that people love to quote things like
# "memory use is O(N^(1/2+o(1)))" then proceed to _clearly_ use N bytes in
# their implementation.  If your data structures between SoA and SoE are the
# same, then all talk about comparative O(blah..blah) memory use is stupid.
#
# Secondly, assuming you're not Dan Bernstein, if your Sieve of Atkin is
# faster than your Sieve of Eratosthenes, then I strongly suggest you verify
# your code actually _works_, and secondly I would bet you made stupid mistakes
# in your SoE implementation.  If your SoA code even remotely resembles the
# Wikipedia code and it comes out faster than SoE, then I _guarantee_ your
# SoE is borked.
#
# SoA does have a slightly better asymptotic operation count O(N/loglogN) vs.
# O(N) for SoE.  The Wikipedia-like code that most people use is O(N) so it
# isn't even theoretically better unless you pull lots of stunts like primegen
# does.  Even if you do, loglogN is essentially a small constant for most uses
# (it's under 4 for all 64-bit values), so you need to make sure all the rest
# of your overhead is controlled.
#
# Sumarizing, in practice the SoE is faster, and often a LOT faster.
#
# </rant>
#
sub atkin2 {
  my($max) = @_;
  return 0 if $max < 2;
  return 1 if $max < 3;

  my @sieve;

  my $sqrt = int(sqrt($max));
  for my $x (1 .. $sqrt) {
    for my $y (1 .. $sqrt) {
      my $n;

      $n = 4*$x*$x + $y*$y;
      if ( ($n <= $max) && ( (($n%12) == 1) || (($n%12) == 5) ) ) {
        $sieve[$n] ^= 1;
      }
      $n = 3*$x*$x + $y*$y;
      if ( ($n <= $max) && (($n%12) == 7) ) {
        $sieve[$n] ^= 1;
      }
      $n = 3*$x*$x - $y*$y;
      if ( ($x > $y) && ($n <= $max) && (($n%12) == 11) ) {
        $sieve[$n] ^= 1;
      }
    }
  }

  for my $n (5 .. $sqrt) {
    if ($sieve[$n]) {
      my $k = $n*$n;
      my $z = $k;
      while ($z <= $max) {
        $sieve[$z] = 0;
        $z += $k;
      }
    }
  }

  $sieve[2] = 1;
  $sieve[3] = 1;
  grep { $sieve[$_] } 2 .. $max;
}

# https://github.com/daoswald/Inline-C-Perl-Mongers-Talk/blob/master/primesbench.pl
sub daos_array {
  my($top) = @_;
  return 0 if $top < 2;
  return 1 if $top < 3;
  $top++;

    my @primes = (1) x $top;
    my $i_times_j;
    for my $i ( 2 .. sqrt $top ) {
        if ( $primes[$i] ) {
            for ( my $j = $i; ( $i_times_j = $i * $j ) < $top; $j++ ) {
                undef $primes[$i_times_j];
            }
        }
    }
  return grep { $primes[$_] } 2 .. $#primes;
}

sub daos_vec {
  my($top) = @_;
  return 0 if $top < 2;
  return 1 if $top < 3;

  my $primes = '';
    vec( $primes, $top, 1 ) = 0;
    my $i_times_j;
    for my $i ( 2 .. sqrt $top ) {
        if ( !vec( $primes, $i, 1 ) ) {
            for ( my $j = $i; ( $i_times_j = $i * $j ) <= $top; $j++ ) {
                vec( $primes, $i_times_j, 1 ) = 1;
            }
        }
    }
  return grep { !vec( $primes, $_, 1 ) } 2 .. $top;
}

# Merlyn's Unix Review Column 26, June 1999
# http://www.stonehenge.com/merlyn/UnixReview/col26.html
sub merlyn {
  my($UPPER) = @_;
  return 0 if $UPPER < 2;
  return 1 if $UPPER < 3;

  my @primes;
  my $sieve = "";
  GUESS: for (my $guess = 2; $guess <= $UPPER; $guess++) {
    next GUESS if vec($sieve,$guess,1);
    push @primes, $guess;
    for (my $mults = $guess * $guess; $mults <= $UPPER; $mults += $guess) {
      vec($sieve,$mults,1) = 1;
    }
  }
  @primes;
}


sub dj1 {
  my($end) = @_;
  return 0 if $end < 2;
  return 1 if $end < 3;

  # vector
  my $sieve = '';
  my $n = 3;
  while ( ($n*$n) <= $end ) {
    my $s = $n*$n;
    while ($s <= $end) {
      vec($sieve, $s >> 1, 1) = 1;
      $s += 2*$n;
    }
    do { $n += 2 } while vec($sieve, $n >> 1, 1) != 0;
  }

  my @primes = (2);
  $n = 3;
  while ($n <= $end) {
    push @primes, $n if !vec($sieve, $n >> 1, 1);
    $n += 2;
  }
  @primes;
}

sub dj2 {
  my($end) = @_;
  return 0 if $end < 2;
  return 1 if $end < 3;

  # array
  my @sieve;
  my $n = 3;
  while ( ($n*$n) <= $end ) {
    my $s = $n*$n;
    while ($s <= $end) {
      $sieve[$s>>1] = 1;
      $s += 2*$n;
    }
    do { $n += 2 } while $sieve[$n>>1];
  }
  my @primes = (2);
  $n = 3;
  while ($n <= $end) {
    push @primes, $n if !$sieve[$n>>1];
    $n += 2;
  }
  @primes;
}

sub dj3 {
  my($end) = @_;
  return 0 if $end < 2;
  return 1 if $end < 3;
  $end-- if ($end & 1) == 0;

  # string
  my $sieve = '1' . '0' x ($end>>1);
  my $n = 3;
  while ( ($n*$n) <= $end ) {
    my $s = $n*$n;
    my $filter_s   = $s >> 1;
    my $filter_end = $end >> 1;
    while ($filter_s <= $filter_end) {
      substr($sieve, $filter_s, 1) = '1';
      $filter_s += $n;
    }
    do { $n += 2 } while substr($sieve, $n>>1, 1);
  }
  my @primes = (2);
  $n = 3-2;
  foreach my $s (split("0", substr($sieve, 1), -1)) {
    $n += 2 + 2 * length($s);
    push @primes, $n if $n <= $end;
  }
  @primes;
}

sub dj4 {
  my($end) = @_;
  return 0 if $end < 2;
  return 1 if $end < 3;
  $end-- if ($end & 1) == 0;

  # string with prefill
  my $whole = int( ($end>>1) / 15);
  my $sieve = '100010010010110' . '011010010010110' x $whole;
  substr($sieve, ($end>>1)+1) = '';
  my $n = 7;
  while ( ($n*$n) <= $end ) {
    my $s = $n*$n;
    my $filter_s   = $s >> 1;
    my $filter_end = $end >> 1;
    while ($filter_s <= $filter_end) {
      substr($sieve, $filter_s, 1) = '1';
      $filter_s += $n;
    }
    do { $n += 2 } while substr($sieve, $n>>1, 1);
  }
  my @primes = (2, 3, 5);
  $n = 7-2;
  foreach my $s (split("0", substr($sieve, 3), -1)) {
    $n += 2 + 2 * length($s);
    push @primes, $n if $n <= $end;
  }
  @primes;
}
