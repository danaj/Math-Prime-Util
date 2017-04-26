#!/usr/bin/env perl
use strict;
use warnings;

use Benchmark qw/:all/;
use Devel::Size qw/total_size/;
use Math::Prime::Util;
*mpu_isprime = \&Math::Prime::Util::is_prime;

my $count = shift || -1;

my @numlist;
#my @testnums = (0..1000, 5_000_000 .. 5_001_000, 30037, 20359*41117, 92987*65171, 27361*31249, 70790191, 3211717*9673231);
my @testnums = (0..1000, 5_000_000 .. 5_001_000, 50_000_000 .. 50_050_000);

my $ip_subs = {
  #"Abigail" => sub { my$r;$r=abigail($_)     for @numlist; $r;},
  "Monks1"  => sub { my$r;$r=monks1($_)      for @numlist; $r;},
  "Rosetta" => sub { my$r;$r=rosetta($_)     for @numlist; $r;},
  "Rosetta2"=> sub { my$r;$r=rosetta2($_)    for @numlist; $r;},
  "DJ"      => sub { my$r;$r=dj($_)          for @numlist; $r;},
  "DJ1"     => sub { my$r;$r=dj1($_)         for @numlist; $r;},
  "DJ2"     => sub { my$r;$r=dj2($_)         for @numlist; $r;},
  "DJ3"     => sub { my$r;$r=dj3($_)         for @numlist; $r;},
  "DJ4"     => sub { my$r;$r=dj4($_)         for @numlist; $r;},
  "MPU"     => sub { my$r;$r=mpu_isprime($_) for @numlist; $r;},
};

my %verify = (
       0 => 0,
       1 => 0,
       2 => 1,
       3 => 1,
       4 => 0,
       5 => 1,
       6 => 0,
       7 => 1,
      13 => 1,
      20 => 0,
     377 => 0,
70790191 => 1,
);

# Verify
while (my($name, $sub) = each (%$ip_subs)) {
  while (my($n, $v_ip) = each (%verify)) {
    @numlist = ($n);
#print "$name($n):  ", $sub->(), "\n";
    my $isprime = ($sub->() ? 1 : 0);
    die "$name($n) = $isprime, should be $v_ip\n" unless $isprime == $v_ip;
  }
}
for my $n (0 .. 50000) {
  die "dj($n) != mpu($n)" unless dj($n) == mpu_isprime($n);
  die "dj1($n) != mpu($n)" unless dj1($n) == mpu_isprime($n);
  die "dj2($n) != mpu($n)" unless dj2($n) == mpu_isprime($n);
  die "dj3($n) != mpu($n)" unless dj3($n) == mpu_isprime($n);
  die "dj4($n) != mpu($n)" unless dj4($n) == mpu_isprime($n);
  die "rosetta($n) != mpu($n)" unless rosetta($n) == mpu_isprime($n)/2;
  die "rosetta2($n) != mpu($n)" unless rosetta2($n) == mpu_isprime($n)/2;
  die "monks1($n) != mpu($n)" unless monks1($n) == mpu_isprime($n)/2;
}
print "Done with verification, starting benchmark\n";

@numlist = @testnums;
cmpthese($count, $ip_subs);


sub rosetta {
  my $n = shift;
  $n % $_ or return 0 for 2 .. sqrt $n;
  $n > 1;
}

sub rosetta2 {
    my $p = shift;
    if ($p == 2) {
        return 1;
    } elsif ($p <= 1 || $p % 2 == 0) {
        return 0;
    } else {
        my $limit = sqrt($p);
        for (my $i = 3; $i <= $limit; $i += 2) {
            return 0 if $p % $i == 0;
        }
        return 1;
    }
}
sub monks1 {
  my $i = shift;
  use POSIX;
  my ($j,$h,$sentinel) = (0,0,0,0);
  return ($i == 2) if $i <= 2;
  # if $i is an even number, it can't be a prime
  if($i%2==0){}
  else {
    $h=POSIX::floor(sqrt($i));
    $sentinel=0;
    # since $i can't be even -> only divide by odd numbers
    for($j=3; $j<=$h; $j+=2){
      if($i%$j==0){
        $sentinel++;
        # $i is not a prime, we can get out of the loop
        $j=$h;
      }
    }
    return 1 if $sentinel == 0;
  }
  0;
}

# Terrifically clever, but useless for large numbers
sub abigail {
  ('1' x shift) !~ /^1?$|^(11+?)\1+$/
}

sub dj {
  my $n = shift;
  return 2 if $n == 2;
  return 0 if $n <= 1 || $n % 2 == 0;
  my $limit = int(sqrt($n));
  for (my $i = 3; $i <= $limit; $i += 2) {
    return 0 if $n % $i == 0;
  }
  2;
}

sub dj1 {
  my($n) = @_;
  return 0 if $n < 2;        # 0 and 1 are composite
  return 2 if ($n == 2) || ($n == 3) || ($n == 5);  # 2, 3, 5 are prime
  # multiples of 2,3,5 are composite
  return 0 if (($n % 2) == 0) || (($n % 3) == 0) || (($n % 5) == 0);

  my $q;
  foreach my $i (qw/7 11 13 17 19 23 29 31 37 41 43 47 53 59/) {
    $q = int($n/$i); return 2 if $q < $i; return 0 if $n == ($q*$i);
  }

  my $i = 61;  # mod-30 loop
  while (1) {
    $q = int($n/$i); last if $q < $i; return 0 if $n == ($q*$i);  $i += 6;
    $q = int($n/$i); last if $q < $i; return 0 if $n == ($q*$i);  $i += 4;
    $q = int($n/$i); last if $q < $i; return 0 if $n == ($q*$i);  $i += 2;
    $q = int($n/$i); last if $q < $i; return 0 if $n == ($q*$i);  $i += 4;
    $q = int($n/$i); last if $q < $i; return 0 if $n == ($q*$i);  $i += 2;
    $q = int($n/$i); last if $q < $i; return 0 if $n == ($q*$i);  $i += 4;
    $q = int($n/$i); last if $q < $i; return 0 if $n == ($q*$i);  $i += 6;
    $q = int($n/$i); last if $q < $i; return 0 if $n == ($q*$i);  $i += 2;
  }
  2;
}

sub dj2 {
  my($n) = @_;
  return 2 if ($n == 2) || ($n == 3) || ($n == 5);  # 2, 3, 5 are prime
  return 0 if $n < 7;  # everything else below 7 is composite
  # multiples of 2,3,5 are composite
  return 0 if (($n % 2) == 0) || (($n % 3) == 0) || (($n % 5) == 0);

  foreach my $i (qw/7 11 13 17 19 23 29 31 37 41 43 47 53 59/) {
    return 2 if $i*$i > $n;
    return 0 if ($n % $i) == 0;
  }
  my $limit = int(sqrt($n));

  my $i = 61;  # mod-30 loop
  while (1) {
    return 0 if ($n % $i) == 0;  $i += 6;  last if $i > $limit;
    return 0 if ($n % $i) == 0;  $i += 4;  last if $i > $limit;
    return 0 if ($n % $i) == 0;  $i += 2;  last if $i > $limit;
    return 0 if ($n % $i) == 0;  $i += 4;  last if $i > $limit;
    return 0 if ($n % $i) == 0;  $i += 2;  last if $i > $limit;
    return 0 if ($n % $i) == 0;  $i += 4;  last if $i > $limit;
    return 0 if ($n % $i) == 0;  $i += 6;  last if $i > $limit;
    return 0 if ($n % $i) == 0;  $i += 2;  last if $i > $limit;
  }
  2;
}

sub dj3 {
  my($n) = @_;
  return 2 if ($n == 2) || ($n == 3) || ($n == 5);  # 2, 3, 5 are prime
  return 0 if $n < 7;  # everything else below 7 is composite
  # multiples of 2,3,5 are composite
  return 0 if (($n % 2) == 0) || (($n % 3) == 0) || (($n % 5) == 0);

  foreach my $i (qw/7 11 13 17 19 23 29 31 37 41 43 47 53 59/) {
    return 2 if $i*$i > $n;
    return 0 if ($n % $i) == 0;
  }
  my $limit = int(sqrt($n));

  my $i = 61;  # mod-30 loop
  while (($i+30) <= $limit) {
    return 0 if ($n % $i) == 0;  $i += 6;
    return 0 if ($n % $i) == 0;  $i += 4;
    return 0 if ($n % $i) == 0;  $i += 2;
    return 0 if ($n % $i) == 0;  $i += 4;
    return 0 if ($n % $i) == 0;  $i += 2;
    return 0 if ($n % $i) == 0;  $i += 4;
    return 0 if ($n % $i) == 0;  $i += 6;
    return 0 if ($n % $i) == 0;  $i += 2;
  }
  while (1) {
    last if $i > $limit;  return 0 if ($n % $i) == 0;  $i += 6;
    last if $i > $limit;  return 0 if ($n % $i) == 0;  $i += 4;
    last if $i > $limit;  return 0 if ($n % $i) == 0;  $i += 2;
    last if $i > $limit;  return 0 if ($n % $i) == 0;  $i += 4;
    last if $i > $limit;  return 0 if ($n % $i) == 0;  $i += 2;
    last if $i > $limit;  return 0 if ($n % $i) == 0;  $i += 4;
    last if $i > $limit;  return 0 if ($n % $i) == 0;  $i += 6;
    last if $i > $limit;  return 0 if ($n % $i) == 0;  $i += 2;
  }
  2;
}

sub dj4 {
  my($n) = @_;
  return 2 if ($n == 2) || ($n == 3) || ($n == 5);  # 2, 3, 5 are prime
  return 0 if $n < 7;  # everything else below 7 is composite
  # multiples of 2,3,5 are composite
  return 0 if (($n % 2) == 0) || (($n % 3) == 0) || (($n % 5) == 0);

  foreach my $i (qw/7 11 13 17 19 23 29/) {
    return 2 if $i*$i > $n;
    return 0 if ($n % $i) == 0;
  }
  my $limit = int(sqrt($n));

  my $i = 31;
  while (($i+30) <= $limit) {
    return 0 if ($n % $i) == 0;  $i += 6;
    return 0 if ($n % $i) == 0;  $i += 4;
    return 0 if ($n % $i) == 0;  $i += 2;
    return 0 if ($n % $i) == 0;  $i += 4;
    return 0 if ($n % $i) == 0;  $i += 2;
    return 0 if ($n % $i) == 0;  $i += 4;
    return 0 if ($n % $i) == 0;  $i += 6;
    return 0 if ($n % $i) == 0;  $i += 2;
  }
  while (1) {
    last if $i > $limit;  return 0 if ($n % $i) == 0;  $i += 6;
    last if $i > $limit;  return 0 if ($n % $i) == 0;  $i += 4;
    last if $i > $limit;  return 0 if ($n % $i) == 0;  $i += 2;
    last if $i > $limit;  return 0 if ($n % $i) == 0;  $i += 4;
    last if $i > $limit;  return 0 if ($n % $i) == 0;  $i += 2;
    last if $i > $limit;  return 0 if ($n % $i) == 0;  $i += 4;
    last if $i > $limit;  return 0 if ($n % $i) == 0;  $i += 6;
    last if $i > $limit;  return 0 if ($n % $i) == 0;  $i += 2;
  }
  2;
}
