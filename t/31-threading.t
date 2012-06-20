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

use Test::More 'tests' => 1;
use Math::Prime::Util ":all";
my $nthreads = 64;

srand(50);
my @digits;
push @digits, random_ndigit_prime(6) for (0..9);

my $tsub = sub { my $sum = 0;  $sum += prime_count($_) for (@digits); return $sum;};
my @threads;
# Fire off all our threads
push @threads, threads->create($tsub) for (1..$nthreads);
# Retrieve results
my $par_sum = 0;
$par_sum += $_->join() for (@threads);

# Now try it on main
my $seq_sum = 0;
$seq_sum += $tsub->() for (1..$nthreads);

is($par_sum, $seq_sum, "$nthreads threads summed prime count");
