#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util qw/:all/;
use Time::HiRes qw(gettimeofday tv_interval);
$| = 1;  # fast pipes

my $nprimes = shift || 50_000_000;

# 1. forprimes does a segmented sieve and calls us for each prime.  This is
#    independent of is_prime and the main sieve.  So for each entry let's
#    compare next_prime and prev_prime.
{
  print "Using MPU forprimes to $nprimes\n";
  my $start_time = [gettimeofday];
  my $nextprint = 5000000;
  my $n = 0;
  forprimes {
    die "next $n not $_" unless next_prime($n) == $_;
    die "prev $n" unless prev_prime($_) == $n;
    $n = $_;
    if ($n > $nextprint) { print "$n..";  $nextprint += 5000000; }
  } $nprimes;
  my $seconds = tv_interval($start_time);
  my $micro_per_call = ($seconds * 1000000) / (2*prime_count($nprimes));
  printf "Success using forprimes to $nprimes.  %6.2f uSec/call\n", $micro_per_call;
}

print "\n";

# 2. Just like before, but now we'll call prime_precalc first.  This makes the
#    prev_prime and next_prime functions really fast since they just look in
#    the cached sieve.
{
  print "Using MPU forprimes to $nprimes with prime_precalc\n";
  my $start_time = [gettimeofday];
  prime_precalc($nprimes);
  my $nextprint = 5000000;
  my $n = 0;
  forprimes {
    die "next $n not $_" unless next_prime($n) == $_;
    die "prev $n" unless prev_prime($_) == $n;
    $n = $_;
    if ($n > $nextprint) { print "$n..";  $nextprint += 5000000; }
  } $nprimes;
  my $seconds = tv_interval($start_time);
  my $micro_per_call = ($seconds * 1000000) / (2*prime_count($nprimes));
  printf "Success using forprimes/precalc to $nprimes.  %6.2f uSec/call\n", $micro_per_call;
}

print "\n\n";

# Now do some more comparative timing.
my @pr = @{primes($nprimes)};
my $numpr = scalar @pr;
prime_memfree();

{
  print "MPU             forprimes...";
  my $start_time = [gettimeofday];
  my $i = 0;
  forprimes {
    die "next $_ not ", $pr[$i-1] unless $pr[$i++] == $_;
  } $nprimes;
  my $seconds = tv_interval($start_time);
  my $micro_per_call = ($seconds * 1000000) / (1*prime_count($nprimes));
  printf "%8.2f uSec/call\n", $micro_per_call;
  prime_memfree();
}
{
  print "MPU             prev/next...";
  my $start_time = [gettimeofday];
  my $n = 0;
  foreach my $p (@pr) {
    my $next = next_prime($n);
    my $prev = prev_prime($p);
    die "MPU next($n) is not $p\n" unless $next == $p;
    die "MPU prev($p) is not $n\n" unless $prev == $n;
    $n = $next;
  }
  my $seconds = tv_interval($start_time);
  my $micro_per_call = ($seconds * 1000000) / (2*$numpr);
  printf "%8.2f uSec/call\n", $micro_per_call;
}
{
  print "MPU precalc     prev/next...";
  my $start_time = [gettimeofday];
  prime_precalc($pr[-1]+1000);
  my $n = 0;
  foreach my $p (@pr) {
    my $next = next_prime($n);
    my $prev = prev_prime($p);
    die "MPU next($n) is not $p\n" unless $next == $p;
    die "MPU prev($p) is not $n\n" unless $prev == $n;
    $n = $next;
  }
  my $seconds = tv_interval($start_time);
  my $micro_per_call = ($seconds * 1000000) / (2*$numpr);
  printf "%8.2f uSec/call\n", $micro_per_call;
  prime_memfree();
}

# Math::Prime::FastSieve
if (eval { require Math::Prime::FastSieve; Math::Prime::FastSieve->import(); Inline->init(); 1; }) {
  print "Math::Prime::FastSieve......";
  my $start_time = [gettimeofday];
  my $sieve = Math::Prime::FastSieve::Sieve->new( $pr[-1]+1000 );
  my $n = 0;
  foreach my $p (@pr) {
    my $next = $sieve->nearest_ge($n+1);
    my $prev = $sieve->nearest_le($p-1);
    die "MPFS next($n) is not $p\n" unless $next == $p;
    die "MPFS prev($p) is not $n\n" unless $prev == $n;
    $n = $next;
  }
  my $seconds = tv_interval($start_time);
  my $micro_per_call = ($seconds * 1000000) / (2*$numpr);
  printf "%8.2f uSec/call\n", $micro_per_call;
} else {
  print "Math::Prime::FastSieve not installed.  Skipping\n";
}

# Math::Pari.
if (eval { require Math::Pari; 1; }) {
  print "Math::Pari      prec/next...";
  my @pari_pr = grep { $_ < 5_000_000 } @pr;
  my $pari_numpr = scalar @pari_pr;
  my $start_time = [gettimeofday];
  my $n = 0;
  foreach my $p (@pari_pr) {
    my $next = Math::Pari::nextprime($n+1);
    my $prev = Math::Pari::precprime($p-1);
    die "Pari next($n) is not $p\n" unless $next == $p;
    die "Pari prec($p) is not $n\n" unless $prev == $n;
    $n = $next;
  }
  my $seconds = tv_interval($start_time);
  my $micro_per_call = ($seconds * 1000000) / (2*$pari_numpr);
  printf "%8.2f uSec/call\n", $micro_per_call;
} else {
  print "Math::Pari not installed.  Skipping\n";
}

# Math::NumSeq::Primes
if (eval { require Math::NumSeq::Primes; 1; }) {
  print "Math::NumSeq::Primes next...";
  my $start_time = [gettimeofday];
  my $seq = Math::NumSeq::Primes->new();
  my $n = 0;
  foreach my $p (@pr) {
    my $next = ($seq->next)[1];
    die "MNP next($n) is not $p\n" unless $next == $p;
    $n = $next;
  }
  my $seconds = tv_interval($start_time);
  my $micro_per_call = ($seconds * 1000000) / (1*$numpr);
  printf "%8.2f uSec/call\n", $micro_per_call;
} else {
  print "Math::NumSeq::Primes not installed.  Skipping\n";
}

# Math::Primality
if (eval { require Math::Primality; 1; }) {
  print "Math::Primality prev/next...";
  my @mp_pr = grep { $_ < 100_000 } @pr;
  my $mp_numpr = scalar @mp_pr;
  my $start_time = [gettimeofday];
  my $n = 0;
  foreach my $p (@mp_pr) {
    my $next = Math::Primality::next_prime($n);
    my $prev = ($p == 2) ? 0 : Math::Primality::prev_prime($p);
    die "MP next($n) is not $p\n" unless $next == $p;
    die "MP prev($p) is not $n\n" unless $prev == $n;
    $n = $next;
  }
  my $seconds = tv_interval($start_time);
  my $micro_per_call = ($seconds * 1000000) / (2*$mp_numpr);
  printf "%8.2f uSec/call\n", $micro_per_call;
} else {
  print "Math::Primality not installed.  Skipping\n";
}
