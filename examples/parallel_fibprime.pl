#!/usr/bin/env perl
use strict;
use warnings;
use threads;
use threads::shared;
use Math::BigInt lib => 'GMP';
use Math::Prime::Util ':all';
use Time::HiRes qw(gettimeofday tv_interval);
$| = 1;

# Find Fibonacci primes in parallel, using Math::Prime::Util and Perl threads.
#
# Dana Jacobsen, 2012.
#
# This will fully utilize however many cores you choose (using the $nthreads
# variable).  It spreads the numbers across threads, where each one runs a
# BPSW test.  A separate thread handles the in-order display.  I have tested
# it on machines with 2, 4, 8, 12, 24, and 64 cores.
#
# On my 12-core computer:
#    24    5387              0.65488
#    25    9311              4.39227
#    26    9677              4.54363
#    27   14431             18.82531
#    28   25561            121.34584
#    29   30757            212.99409
#    30   35999            376.59567
#    31   37511            432.10713
#    32   50833           1151.85562
#
# Though not as pretty as the Haskell solution on haskell.org, it is a
# different way of solving the problem that is faster and more scalable.

my $time_start = [gettimeofday];
my $nthreads = 12;
prime_precalc(1_000_000);

my @found :shared;     # push the primes found here
my @karray : shared;   # array of min k for each thread

my @threads;
push @threads, threads->create('fibprime', $_) for (1..$nthreads);

# Let the threads work for a little before starting the display loop
sleep 2;
my $n = 0;
lock(@karray);
while (1) {
  cond_wait(@karray);
  {
    lock(@found);
    next if @found == 0;
    # Someone has found a result.  Discover min k processed so far.
    my $mink = $karray[1] || 0;
    for my $t (2..$nthreads) {
      my $progress = $karray[$t] || 0;
      $mink = $progress if $progress < $mink;
    }
    next unless $mink > 0;  # someone hasn't even started
    @found = sort { (split(/ /, $a))[0] <=> (split(/ /, $b))[0] } @found;
    while ( @found > 0 && (split(/ /, $found[0]))[0] <= $mink ) {
      my($k, $time_int) = split(/ /, shift @found);
      printf "%3d %7d %20.5f\n", ++$n, $k, $time_int;
    }
  }
}
$_->join() for (@threads);

sub fib_n {
  my ($n, $fibstate) = @_;
  @$fibstate = (1, Math::BigInt->new(0), Math::BigInt->new(1))
     unless defined $fibstate->[0];
  my ($curn, $a, $b) = @$fibstate;
  die "fib_n only increases" if $n < $curn;
  do { ($a, $b) = ($b, $a+$b); } for (1 .. $n-$curn);
  @$fibstate = ($n, $a, $b);
  $b;
}

sub fibprime {
  my $tnum = shift;
  my @fibstate;
  my $nth = $tnum;
  while (1) {
    # Exploit knowledge that excepting k=4, all prime F_k have a prime k.
    my $k = ($nth <= 2) ?  2 + $nth  :  nth_prime($nth);
    $nth += $nthreads;
    my $Fk = fib_n($k, \@fibstate);
    if (is_prob_prime($Fk)) {
      lock(@found);
      push @found, $k . " " . tv_interval($time_start);
    }
    {
      lock(@karray);
      $karray[$tnum] = $k;
      cond_signal(@karray);
    }
  }
}
