#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util qw/:all/;
use Time::HiRes qw(gettimeofday tv_interval);
$| = 1;  # fast pipes

my $mpu_limit = shift || 50_000_000;

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
  printf "Success using forprimes to $mpu_limit.  %6.2f uSec/call\n", $micro_per_call;
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
  printf "Success using forprimes/precalc to $mpu_limit.  %6.2f uSec/call\n", $micro_per_call;
}

print "\n\n";

# Now do some more comparative timing.
my @pr = @{primes(10_000_000)};
my $numpr = scalar @pr;
prime_memfree();

{
  print "MPU             prev/next...";
  my $start_time = [gettimeofday];
  my $n = 0;
  foreach my $p (@pr) {
    my $next = next_prime($n);
    die "MPU next($n) is not $p\n" unless $next == $p;
    die "MPU prev($p) is not $n\n" unless $n == prev_prime($p);
    $n = $next;
  }
  my $seconds = tv_interval($start_time);
  my $micro_per_call = ($seconds * 1000000) / (2*prime_count($numpr));
  printf "%8.2f uSec/call\n", $micro_per_call;
}
{
  print "MPU precalc     prev/next...";
  my $start_time = [gettimeofday];
  prime_precalc($pr[-1]+1000);
  my $n = 0;
  foreach my $p (@pr) {
    my $next = next_prime($n);
    die "MPU next($n) is not $p\n" unless $next == $p;
    die "MPU prev($p) is not $n\n" unless $n == prev_prime($p);
    $n = $next;
  }
  my $seconds = tv_interval($start_time);
  my $micro_per_call = ($seconds * 1000000) / (2*prime_count($numpr));
  printf "%8.2f uSec/call\n", $micro_per_call;
  prime_memfree();
}

# 3. Now Math::Pari.
if (eval { require Math::Pari; 1; }) {
  print "Math::Pari      prec/next...";
  my $start_time = [gettimeofday];
  my $n = 0;
  foreach my $p (@pr) {
    my $next = Math::Pari::nextprime($n+1);
    die "MPU next($n) is not $p\n" unless $next == $p;
    die "MPU prev($p) is not $n\n" unless $n == Math::Pari::precprime($p-1);
    $n = $next;
  }
  my $seconds = tv_interval($start_time);
  my $micro_per_call = ($seconds * 1000000) / (2*prime_count($numpr));
  printf "%8.2f uSec/call\n", $micro_per_call;
} else {
  print "Math::Pari not installed.  Skipping\n";
}

# 4. Math::Primality
if (eval { require Math::Primality; 1; }) {
  print "Math::Primality prev/next...";
  my @mppr = @pr[0..50000];
  my $nummppr = scalar @mppr;
  my $start_time = [gettimeofday];
  my $n = 0;
  foreach my $p (@mppr) {
    my $next = Math::Primality::next_prime($n);
    my $prev = ($p == 2) ? 0 : Math::Primality::prev_prime($p);
    die "MPU next($n) is not $p\n" unless $next == $p;
    die "MPU prev($p) is not $n\n" unless $n == $prev;
    $n = $next;
  }
  my $seconds = tv_interval($start_time);
  my $micro_per_call = ($seconds * 1000000) / (2*prime_count($nummppr));
  printf "%8.2f uSec/call\n", $micro_per_call;
} else {
  print "Math::Primality not installed.  Skipping\n";
}

