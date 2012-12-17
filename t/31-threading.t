#!/usr/bin/env perl
use strict;
use warnings;

use Config;
BEGIN {
  if (! $Config{useithreads} || $] < 5.008) {
    print("1..0 # Skip Threads not supported\n");
    exit(0);
  }
  # All these tests used to run on Cygwin, but now they're all giving me
  # random panics in mutexes.  Same with NetBSD.
  if ($Config{osname} eq 'cygwin') {
    print "1..0 # Skip Cygwin threads are unstable\n";
    exit 0;
  }
  if ($Config{osname} eq 'netbsd') {
    print "1..0 # Skip NetBSD threads have panic issues\n";
    exit 0;
  }
  # Should be be looking for newer than 5.008?
  if (! eval { require threads }) {
    print "1..0 # Skip threads.pm not installed\n";
    exit 0;
  }
}

use Test::More 'tests' => 9;
use Math::Prime::Util ":all";

my $extra = defined $ENV{RELEASE_TESTING} && $ENV{RELEASE_TESTING};
my $is_win32 = ($Config{osname} eq 'MSWin32')  ?  1  :  0;

# threads are memory hogs, so we want few of them.  But for testing purposes,
# we want a lot of them.  4-8 perhaps.
my $numthreads = 4;

# Random numbers, pregenerated
my @randn = (
 qw/71094  1864 14650 58418 46196 45215 70355 80402 70420 33457 73424 45744
    22229 61529 82574 61578 26369 76750 15724 61272 52128 77252  2207  3579
    69734 14488 20846 46906  6992 43938 34945 51978 11336 58462 11973 75248
    39165  8147 62890 63565 39279 47830 43617 40299 65372 37479   884 27007
    24978 55716 38115 71502 30134 40903 71231 40095  9054 54133 13876 55660
    44544  1880 39217 36609 38711 49576 55029 21831 75022 69128  2311 16321
     1400  9659  6010  8206 78113 76539 17430 69393 26519 50162 49317 20231/);

if ($extra) {
  $numthreads *= 2;
  push @randn, (qw/
    11019 28515 73527 50147 33512 28347 19122 66580 14286 81842 38344 10886
    52253 57834 37446 49360 24401 45815 54811  1703 38180 22473 17946 58833
    29700 55366 35155 31902 28299 34139 51961 75210  9126 30331 54658 50208
    13936 57086 27118 75817 31571 76715 53441 31118 22091 47356 67284 37756
    67826   819 78844 64174 53566 28410 40455 76690 69141  2001  1251 39140
     2328 49159 14379 73801 30195 19745 72355 51038 76557 63516 54486 45951
    65586 61730  6310 73490 71132 25970 51034 27856 11490 25817 24283 52759
    68248  9246 52896 72365 31983 74001 16616 63960 70718 43518 27054  6397
     1247 64241 27517  2927  3557 76192 36376 21334  1395 20926 36088 65519
     2650  9739 23758 74720 34458 41096 51926 45932 14850 38181 60833 53481
     8086 43319 11891 22674 22916 72918  3650 35246 39543 25544 35578 67023
    50752 29653 76351 64909  9425 27547 10108 13399 69540  3833 12748  6386
    76511 28041 31586 50034  8828 17845 44376 74301 39762 40216  5092 16261
     7434 29908 18671  7189 18373 31608 67155 19129 20586  6713 73424 20568
    64299 71411 53762 20070 56014  3560  9129 50993 44983 15434  5035 77815
    22836  9786 24808 50756 15298 48358 36466  4308   195 69058 55813 18840
    23284 41448 37349 59268 36894 79674 31694 73975 71738 18344 26328  5264
    79976 26714 27187 65237 18881 74975 28505 16271 51390 22598 65689 65512
    20357 68743 72422 69481 26714  6824 30012/);
}


thread_test(
  sub { my $sum = 0;  $sum += prime_count($_) for (@randn); return $sum;},
  $numthreads, "sum prime_count");

SKIP: {
  skip "Win32 needs precalc, skipping alloc/free stress test", 1 if $is_win32;

  thread_test(
    sub { my $sum = 0;  for (@randn) {$sum += prime_count($_); prime_memfree; } return $sum;},
    $numthreads, "sum prime_count with overlapping memfree calls");
}

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

# Override rand because some systems don't use thread-safe rand, which will
# make this test go beserk (it doesn't make *sense* without TS rand).
# This should provide a semi-random result that depends on the seed given.
# rand should return a number in the range [0,1) --- 0 is ok, 1 is not.  Since
# this is all in Perl, it all gets duped with the interpreter so is TS.
{
  my $seed = 1;
  sub mysrand { $seed = $_[0]; }
  sub rand { $seed = (1103515245*$seed + 12345) & 0xFFFF_FFFF; $seed/(0xFFFF_FFFF+1); }
}

thread_test(
  sub { my $sum = 0;  for (@randn) { mysrand($_); $sum += random_ndigit_prime(6); } return $sum;},
  $numthreads, "random 7-digit prime");

thread_test(
  sub { my $sum = 0;  $sum += int(RiemannR($_)) for (@randn); return $sum;},
  $numthreads, "RiemannR");

sub thread_test {
  my $tsub = shift;
  my $nthreads = shift;
  my $text = shift;

  if ($is_win32) {
    prime_precalc(1_200_000);  # enough for all these tests
  }

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
