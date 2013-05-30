#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util qw/:all/;
use Time::HiRes qw(gettimeofday tv_interval);
$| = 1;  # fast pipes

my $mpu_limit = shift || 50_000_000;
my $mp_limit = shift || 20_000;

# 1. forprimes does a segmented sieve and calls us for each prime.  This is
#    independent of is_prime and the main sieve.  So for each entry let's
#    compare next_prime and prev_prime.
{
  print "Using MPU forprimes to $mpu_limit\n";
  my $start_time = [gettimeofday];
  my $nextprint = 5000000;
  my $n = 0;
  forprimes {
    die "next $n not $_" unless next_prime($n) == $_;
    die "prev $n" unless prev_prime($_) == $n;
    $n = $_;
    if ($n > $nextprint) { print "$n..";  $nextprint += 5000000; }
  } $mpu_limit;
  my $seconds = tv_interval($start_time);
  my $micro_per_call = ($seconds * 1000000) / (2*prime_count($mpu_limit));
  printf "\nSuccess using forprimes to $mpu_limit.  %6.2f uSec/call\n", $micro_per_call;
}

print "\n";

# 2. Just like before, but now we'll call prime_precalc first.  This makes the
#    prev_prime and next_prime functions really fast since they just look in
#    the cached sieve.
{
  print "Using MPU forprimes to $mpu_limit with prime_precalc\n";
  my $start_time = [gettimeofday];
  prime_precalc($mpu_limit);
  my $nextprint = 5000000;
  my $n = 0;
  forprimes {
    die "next $n not $_" unless next_prime($n) == $_;
    die "prev $n" unless prev_prime($_) == $n;
    $n = $_;
    if ($n > $nextprint) { print "$n..";  $nextprint += 5000000; }
  } $mpu_limit;
  my $seconds = tv_interval($start_time);
  my $micro_per_call = ($seconds * 1000000) / (2*prime_count($mpu_limit));
  printf "\nSuccess using forprimes/precalc to $mpu_limit.  %6.2f uSec/call\n", $micro_per_call;
}

print "\n";

# 3. Now we'll use Math::Primality to compare next_prime, prev_prime, and
#    is_prime.
{
  require Math::Primality;
  print "Using Math::Primality to $mp_limit\n";
  my $start_time = [gettimeofday];
  foreach my $n (0 .. $mp_limit) {
    die "next $n" unless next_prime($n) == Math::Primality::next_prime($n);
    if ($n <= 2) {
      die "prev $n" unless prev_prime($n) == 0;
    } else {
      die "prev $n" unless prev_prime($n) == Math::Primality::prev_prime($n);
    }
    die "is $n" unless is_prime($n) == Math::Primality::is_prime($n);
    print "$n.." unless $n % 10000;
  }
  my $seconds = tv_interval($start_time);
  my $micro_per_call = ($seconds * 1000000) / (6*prime_count($mp_limit));
  printf "\nSuccess using Math::Primality to $mp_limit.  %6.2f uSec/call\n", $micro_per_call;
}
