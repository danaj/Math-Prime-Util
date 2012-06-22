#!/usr/bin/env perl
use strict;
use warnings;

use Config;
BEGIN {
  if (! $Config{useithreads} || $] < 5.008) {
    print("1..0 # Skip Threads not supported\n");
    exit(0);
  }
  # Should be be looking for newer than 5.008?
  if (! eval { require threads }) {
    print "1..0 # Skip threads.pm not installed\n";
    exit 0;
  }
}

use Test::More 'tests' => 9;
use Math::Prime::Util ":all";
my $numthreads = 8;

srand(50);
my @randn;
push @randn, rand(100000) for (0..377);

thread_test(
  sub { my $sum = 0;  $sum += prime_count($_) for (@randn); return $sum;},
  $numthreads, "sum prime_count");

thread_test(
  sub { my $sum = 0;  for (@randn) {$sum += prime_count($_); prime_memfree; } return $sum;},
  $numthreads, "sum prime_count with overlapping memfree calls");

thread_test(
  sub { my $sum = 0; for my $d (@randn) { for my $f (factor($d)) { $sum += $f; } } return $sum; },
  $numthreads, "factor");

thread_test(
  sub { my $sum = 0;  $sum += nth_prime($_) for (@randn); return $sum;},
  $numthreads, "nth_prime");

thread_test(
  sub { my $sum = 0;  $sum += next_prime($_) for (@randn); return $sum;},
  $numthreads, "next_prime");

thread_test(
  sub { my $sum = 0;  $sum += prev_prime($_) for (@randn); return $sum;},
  $numthreads, "prev_prime");

thread_test(
  sub { my $sum = 0;  $sum += is_prime($_) for (@randn); return $sum;},
  $numthreads, "is_prime");

thread_test(
  sub { my $sum = 0;  for (@randn) { srand($_); $sum += random_ndigit_prime(6); } return $sum;},
  $numthreads, "random 7-digit prime");

thread_test(
  sub { my $sum = 0;  $sum += int(RiemannR($_)) for (@randn); return $sum;},
  $numthreads, "RiemannR");

sub thread_test {
  my $tsub = shift;
  my $nthreads = shift;
  my $text = shift;

  my @threads;
  # Fire off all our threads
  push @threads, threads->create($tsub) for (1..$nthreads);
  # Get results
  my $par_sum = 0;
  $par_sum += $_->join() for (@threads);
  prime_memfree;

  # Now do the same operation sequentially
  my $seq_sum = 0;
  $seq_sum += $tsub->() for (1..$nthreads);
  prime_memfree;

  is($par_sum, $seq_sum, "$nthreads threads $text");
}
