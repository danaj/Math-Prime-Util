#!/usr/bin/env perl
use strict;
use warnings;

use Benchmark qw/:all/;

my $upper = shift || 8192;
my $count = shift || -1;
my $countarg;

my $pc_subs = {
  "Scriptol"    => sub { scriptol($countarg) },
  "Shootout"    => sub { shootout($countarg) },
  "In Many"     => sub { inmany($countarg) },
  "Rosetta 1"   => sub { rosetta1($countarg) },
  "Rosetta 2"   => sub { rosetta2($countarg) },
  "Rosetta 3"   => sub { rosetta3($countarg) },
  "Rosetta 4"   => sub { rosetta4($countarg) },
  "DJ Vec"      => sub { dj1($countarg) },
  "DJ Array"    => sub { dj2($countarg) },
  "DJ String1"  => sub { dj3($countarg) },
  "DJ String2"  => sub { dj4($countarg) },
};

my %verify = (
      10 => 4,
     100 => 25,
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
  return 0 if $upper < 2;
  return 1 if $upper < 3;

  my @flags = (0 .. $max);
  for my $i (2 .. int(sqrt($max)) + 1)
  {
       next unless defined $flags[$i];
       for (my $k=$i+$i; $k <= $max; $k+=$i)
       {
           undef $flags[$k];
       }
  }
  my $count = 0;
  for my $j (2 .. $max) {
    $count++ if defined $flags[$j];
  }
  $count;
}

# http://dada.perl.it/shootout/sieve.perl.html
sub shootout {
  my($max) = @_;
  return 0 if $upper < 2;
  return 1 if $upper < 3;

  my $count = 0;
  my @flags = (0 .. $max);
  for my $i (2 .. $max) {
       next unless defined $flags[$i];
       for (my $k=$i+$i; $k <= $max; $k+=$i) {
           undef $flags[$k];
       }
       $count++;
  }
  $count;
}

# http://c2.com/cgi/wiki?SieveOfEratosthenesInManyProgrammingLanguages
sub inmany {
  my($max) = @_;
  return 0 if $upper < 2;
  return 1 if $upper < 3;

  my @c;
  for(my $t=3; $t*$t<$max; $t+=2) {
     if (!$c[$t]) {
         for(my $s=$t*$t; $s<$max; $s+=$t*2) { $c[$s]++ }
     }
 }
 my $count = 1;
 for(my $t=3; $t<$max; $t+=2) {
     $c[$t] || $count++;
 }
 $count;
}

# http://rosettacode.org/wiki/Sieve_of_Eratosthenes#Perl
sub rosetta1 {
  my($max) = @_;
  return 0 if $upper < 2;
  return 1 if $upper < 3;

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
  $count;  #scalar @primes;
}

# http://rosettacode.org/wiki/Sieve_of_Eratosthenes#Perl
sub rosetta2 {
  my($max) = @_;
  return 0 if $upper < 2;
  return 1 if $upper < 3;

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
  $count; #scalar @primes;
}

# http://rosettacode.org/wiki/Sieve_of_Eratosthenes#Perl
sub rosetta3 {
  my($max) = @_;
  return 0 if $upper < 2;
  return 1 if $upper < 3;

  my $i;
  my @s;
  my $count = scalar
              grep { not $s[ $i  = $_ ] and do
		 { $s[ $i += $_ ]++ while $i <= $max; 1 }
	} 2 .. $max;
  $count; #scalar @primes;
}

# http://rosettacode.org/wiki/Sieve_of_Eratosthenes#Perl
sub rosetta4 {
  my($max) = @_;
  return 0 if $upper < 2;
  return 1 if $upper < 3;

  my $i;
  my $s = '';
  my $count = scalar
              grep { not vec $s, $i  = $_, 1 and do
		 { (vec $s, $i += $_, 1) = 1 while $i <= $max; 1 }
	} 2 .. $max;
  $count; #scalar @primes;
}

sub dj1 {
  my($end) = @_;
  return 0 if $upper < 2;
  return 1 if $upper < 3;

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
  return 0 if $upper < 2;
  return 1 if $upper < 3;

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
  return 0 if $upper < 2;
  return 1 if $upper < 3;
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
  my $count = 1 + $sieve =~ tr/0//;
  $count;
}

# 2-3x faster than inmany, 6-7x faster than any of the other non-DJ methods.
sub dj4 {
  my($end) = @_;
  return 0 if $upper < 2;
  return 1 if $upper < 3;
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
  my $count = 1 + $sieve =~ tr/0//;
  $count;
}
