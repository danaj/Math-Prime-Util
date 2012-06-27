#!/usr/bin/env perl
use strict;
use warnings;

use Benchmark qw/:all/;
#use Devel::Size qw/total_size/;
#use Math::Prime::Util;
#use Math::Prime::FastSieve;
#*mpu_erat = \&Math::Prime::Util::erat_primes;
#*fs_erat = \&Math::Prime::FastSieve::primes;

my $upper = shift || 8192;
my $count = shift || -1;
my $countarg;

# Shows sizes for sieving to 100k, and rate/second for sieving to 16k
my $pc_subs = {
  "Scriptol"    => sub { scriptol($countarg) },    # 3200k    290/s
  "Shootout"    => sub { shootout($countarg) },    # 3200k    231/s
  "In Many"     => sub { inmany($countarg) },      # 2018k    666/s
  "Rosetta 1"   => sub { rosetta1($countarg) },    # 3449k    187/s
  "Rosetta 2"   => sub { rosetta2($countarg) },    #   13k    109/s
  "Rosetta 3"   => sub { rosetta3($countarg) },    # 4496k    165/s
  "Rosetta 4"   => sub { rosetta4($countarg) },    #   25k     60/s
  "Atkin MPTA"  => sub { atkin($countarg) },       # 3430k     90/s
  "DO Array"    => sub {daoswald_array($countarg)},# 3200k    306/s
  "DO Vec"      => sub {daoswald_vec($countarg)},  #   13k    112/s
  "Merlyn"      => sub { merlyn($countarg)},       #   13k     96/s
  "DJ Vec"      => sub { dj1($countarg) },         #    7k    245/s
  "DJ Array"    => sub { dj2($countarg) },         # 1494k    475/s
  "DJ String1"  => sub { dj3($countarg) },         #   50k    981/s
  "DJ String2"  => sub { dj4($countarg) },         #   50k   1682/s
#  "MPU Sieve"   => sub {
#               scalar @{mpu_erat(2,$countarg)}; }, #    3k  14325/s
#  "MPFS Sieve"  => sub {
#               scalar @{fs_erat($countarg)};    }, #    7k  14325/s
};

my %verify = (
      10 => 4,
      11 => 5,
     100 => 25,
     112 => 29,
     113 => 30,
     114 => 30,
    1000 => 168,
   10000 => 1229,
  100000 => 9592,
);

# Verify
while (my($name, $sub) = each (%$pc_subs)) {
  while (my($n, $pin) = each (%verify)) {
    $countarg = $n;
    my $picount = $sub->();
    die "$name ($n) = $picount, should be $pin" unless $picount == $pin;
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
  #print "scriptol size: ", total_size(\@flags), "\n" if $max > 90000;
  my $count = 0;
  for my $j (2 .. $max) {
    $count++ if defined $flags[$j];
  }
  $count;
}

# http://dada.perl.it/shootout/sieve.perl.html
sub shootout {
  my($max) = @_;
  return 0 if $max < 2;
  return 1 if $max < 3;

  my $count = 0;
  my @flags = (0 .. $max);
  for my $i (2 .. $max) {
       next unless defined $flags[$i];
       for (my $k=$i+$i; $k <= $max; $k+=$i) {
           undef $flags[$k];
       }
       $count++;
  }
  #print "shootout size: ", total_size(\@flags), "\n" if $max > 90000;
  $count;
}

# http://c2.com/cgi/wiki?SieveOfEratosthenesInManyProgrammingLanguages
sub inmany {
  my($max) = @_;
  return 0 if $max < 2;
  return 1 if $max < 3;
  $max++;

  my @c;
  for(my $t=3; $t*$t<$max; $t+=2) {
     if (!$c[$t]) {
         for(my $s=$t*$t; $s<$max; $s+=$t*2) { $c[$s]++ }
     }
 }
 #print "inmany size: ", total_size(\@c), "\n" if $max > 90000;
 my $count = 1;
 for(my $t=3; $t<$max; $t+=2) {
     $c[$t] || $count++;
 }
 $count;
}

# http://rosettacode.org/wiki/Sieve_of_Eratosthenes#Perl
sub rosetta1 {
  my($max) = @_;
  return 0 if $max < 2;
  return 1 if $max < 3;

  my $count = 0;  #my @primes;
  my @tested = (1);
  my $j = 1;
  while ($j < $max) {
    next if $tested[$j++];
    $count++;  #push @primes, $j;
    for (my $k= $j; $k <= $max; $k+=$j) {
      $tested[$k-1]= 1;
    }
  }
  #print "R1 size: ", total_size(\@tested), "\n" if $max > 90000;
  $count;  #scalar @primes;
}

# http://rosettacode.org/wiki/Sieve_of_Eratosthenes#Perl
sub rosetta2 {
  my($max) = @_;
  return 0 if $max < 2;
  return 1 if $max < 3;

  my $count = 0; #my @primes;
  my $nonPrimes = '';
  foreach my $p (2 .. $max) {
    unless (vec($nonPrimes, $p, 1)) {
      for (my $i = $p * $p; $i <= $max; $i += $p) {
        vec($nonPrimes, $i, 1) = 1;
      }
      $count++; #push @primes, $p;
    }
  }
  #print "R2 size: ", total_size(\$nonPrimes), "\n" if $max > 90000;
  $count; #scalar @primes;
}

# http://rosettacode.org/wiki/Sieve_of_Eratosthenes#Perl
sub rosetta3 {
  my($max) = @_;
  return 0 if $max < 2;
  return 1 if $max < 3;

  my $i;
  my @s;
  my $count = scalar
              grep { not $s[ $i  = $_ ] and do
		 { $s[ $i += $_ ]++ while $i <= $max; 1 }
	} 2 .. $max;
  #print "R3 size: ", total_size(\@s), "\n" if $max > 90000;
  $count; #scalar @primes;
}

# http://rosettacode.org/wiki/Sieve_of_Eratosthenes#Perl
sub rosetta4 {
  my($max) = @_;
  return 0 if $max < 2;
  return 1 if $max < 3;

  my $i;
  my $s = '';
  my $count = scalar
              grep { not vec $s, $i  = $_, 1 and do
		 { (vec $s, $i += $_, 1) = 1 while $i <= $max; 1 }
	} 2 .. $max;
  #print "R4 size: ", total_size(\$s), "\n" if $max > 90000;
  $count; #scalar @primes;
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
    $sieve{2} = 1;
    $sieve{3} = 1;
    #print "Atkin size: ", total_size(\%sieve), "\n" if $max > 90000;

    # save the found primes in our cache
    my $count = 0;
    foreach my $n ( 2 .. $max ) {
        next unless $sieve{$n};
        $count++;
    }
    $count;
}

# https://github.com/daoswald/Inline-C-Perl-Mongers-Talk/blob/master/primesbench.pl
sub daoswald_array {
  my($top) = @_;
  return 0 if $top < 2;
  return 1 if $top < 3;

    my @primes = (1) x $top;
    my $i_times_j;
    for my $i ( 2 .. sqrt $top ) {
        if ( $primes[$i] ) {
            for ( my $j = $i; ( $i_times_j = $i * $j ) < $top; $j++ ) {
                undef $primes[$i_times_j];
            }
        }
    }
  #print "do_array size: ", total_size(\@primes), "\n" if $top > 90000;
  my $count = scalar grep { $primes[$_] } 2 .. $#primes;
  $count;
}

sub daoswald_vec {
  my($top) = @_;
  return 0 if $top < 2;
  return 1 if $top < 3;

  my $primes = '';
    vec( $primes, $top, 1 ) = 0;
    my $i_times_j;
    for my $i ( 2 .. sqrt $top ) {
        if ( !vec( $primes, $i, 1 ) ) {
            for ( my $j = $i; ( $i_times_j = $i * $j ) < $top; $j++ ) {
                vec( $primes, $i_times_j, 1 ) = 1;
            }
        }
    }
  #print "do_vec size: ", total_size(\$primes), "\n" if $top > 90000;
  my $count = scalar grep { !vec( $primes, $_, 1 ) } 2 .. $top-1 ;
  $count;
}

# Merlyn's Unix Review Column 26, June 1999
# http://www.stonehenge.com/merlyn/UnixReview/col26.html
sub merlyn {
  my($UPPER) = @_;
  return 0 if $UPPER < 2;
  return 1 if $UPPER < 3;

  my $count = 0;
  my $sieve = "";
  GUESS: for (my $guess = 2; $guess <= $UPPER; $guess++) {
    next GUESS if vec($sieve,$guess,1);
    $count++;
    for (my $mults = $guess * $guess; $mults <= $UPPER; $mults += $guess) {
      vec($sieve,$mults,1) = 1;
    }
  }
  #print "Merlyn size: ", total_size(\$sieve), "\n" if $UPPER > 90000;
  $count;
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
  #print "DJ1 size: ", total_size(\$sieve), "\n" if $end > 90000;
  my $count = 1;
  $n = 3;
  while ($n <= $end) {
    $count++ if !vec($sieve, $n >> 1, 1);
    $n += 2;
  }
  $count;
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
  #print "DJ2 size: ", total_size(\@sieve), "\n" if $end > 90000;
  my $count = 1;
  $n = 3;
  while ($n <= $end) {
    $count++ if !$sieve[$n>>1];
    $n += 2;
  }
  $count;
}

# ~2x faster than inmany, lots faster than the others.  Only loses to dj4,
# which is just this code with a presieve added.
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
  #print "DJ3 size: ", total_size(\$sieve), "\n" if $end > 90000;
  my $count = 1 + $sieve =~ tr/0//;
  $count;
}

# 2-3x faster than inmany, 6-7x faster than any of the other non-DJ methods.
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
  #print "DJ4 size: ", total_size(\$sieve), "\n" if $end > 90000;
  my $count = 1 + $sieve =~ tr/0//;
  $count;
}
