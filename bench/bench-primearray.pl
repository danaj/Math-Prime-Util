#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util qw/:all/;
use Math::Prime::Util::PrimeArray;
use Math::NumSeq::Primes;
use Math::Prime::TiedArray;
use Benchmark qw/:all/;
use List::Util qw/min max/;
my $count = shift || -2;

my ($s, $nlimit, $ilimit, $expect);

if (1) {
print '-' x 79, "\n";
print "summation to 100k, looking for best methods (typically slice)\n";
$nlimit = 100000;
$ilimit = prime_count($nlimit)-1;
$expect = 0; forprimes { $expect += $_ } $nlimit;

cmpthese($count,{
  'pa fetch'  => sub { $s=0; my $o = tie my @p, "Math::Prime::Util::PrimeArray";
                       $s += $o->FETCH($_) for 0..$ilimit;
                       die unless $s == $expect; },
  'pa index'  => sub { $s=0; tie my @primes, "Math::Prime::Util::PrimeArray";
                       $s += $primes[$_] for 0..$ilimit;
                       die unless $s == $expect; },
  'pa loop'   => sub { $s=0; tie my @primes, "Math::Prime::Util::PrimeArray";
                       for (@primes) { last if $_ > $nlimit; $s += $_; }
                       die $s unless $s == $expect; },
  'pa slice'  => sub { $s=0; tie my @primes, "Math::Prime::Util::PrimeArray";
                       $s += $_ for @primes[0..$ilimit];
                       die unless $s == $expect; },
  'pa each'   => sub { $s=0; tie my @primes, "Math::Prime::Util::PrimeArray";
                       # Note: using last inside each is Very Bad Stuff.
                       while(my(undef,$v) = each @primes) { last if $v > $nlimit; $s += $v; }
                       die $s unless $s == $expect; },
  'pa shift'  => sub { $s=0; tie my @primes, "Math::Prime::Util::PrimeArray";
                       while ((my $p = shift @primes) <= $nlimit) { $s += $p; }
                       die unless $s == $expect; },
});
}

if (1) {
print '-' x 79, "\n";
print "summation to 100k, looking for best MPTA extension (typically ~1000)\n";
$nlimit = 100000;
$ilimit = prime_count($nlimit)-1;
$expect = 0; forprimes { $expect += $_ } $nlimit;

cmpthese($count,{
  'MPTA'       => sub { $s=0; tie my @primes, "Math::Prime::TiedArray";
                       $s += $primes[$_] for 0..$ilimit;
                       die unless $s == $expect; },
  'MPTA 400'   => sub { $s=0; tie my @primes, "Math::Prime::TiedArray", extend_step => 400;
                       $s += $primes[$_] for 0..$ilimit;
                       die unless $s == $expect; },
  'MPTA 1000'  => sub { $s=0; tie my @primes, "Math::Prime::TiedArray", extend_step => 1000;
                       $s += $primes[$_] for 0..$ilimit;
                       die unless $s == $expect; },
  'MPTA 4000'  => sub { $s=0; tie my @primes, "Math::Prime::TiedArray", extend_step => 4000;
                       $s += $primes[$_] for 0..$ilimit;
                       die unless $s == $expect; },
});
}

if (1) {
print '-' x 79, "\n";
print "summation to 100k\n";
print "Note: MPU::PrimeArray is about 30x faster than MPTA here.\n";
print "      Math::NumSeq::Primes is reasonable fast (not random access)\n";
print "      MPU's forprimes smashes everything else (not random access)\n";
$nlimit = 100000;
$ilimit = prime_count($nlimit)-1;
$expect = 0; forprimes { $expect += $_ } $nlimit;

cmpthese($count,{
  'primes'    => sub { $s=0; $s += $_ for @{primes($nlimit)}; die unless $s == $expect; },
  'forprimes' => sub { $s=0; forprimes { $s += $_ } $nlimit;  die unless $s == $expect; },
  'iterator'  => sub { $s=0; my $it = prime_iterator();
                       $s += $it->() for 0..$ilimit;
                       die unless $s == $expect; },
  'OO iter'   => sub { $s=0; my $it = prime_iterator_object();
                       $s += $it->iterate() for 0..$ilimit;
                       die unless $s == $expect; },
  'pa slice'  => sub { $s=0; tie my @primes, "Math::Prime::Util::PrimeArray";
                       $s += $_ for @primes[0..$ilimit];
                       die unless $s == $expect; },
  'NumSeq'    => sub { $s=0; my $seq = Math::NumSeq::Primes->new;
                       while (1) { my($undev,$v) = $seq->next; last if $v > $nlimit; $s += $v; }
                       die $s unless $s == $expect; },
  # This was slightly faster than slice or shift
  'MPTA'      => sub { $s=0; tie my @primes, "Math::Prime::TiedArray", extend_step => 1000;
                       $s += $primes[$_] for 0..$ilimit;
                       die unless $s == $expect; },
});
}

if (0) {
print '-' x 79, "\n";
print "summation to 10M\n";
print "Note: Math::Prime::TiedArray takes too long\n";
print "      Math::NumSeq::Primes is now ~2x slower than PrimeArray\n";
print "      forprimes is still the fastest solution for sequential access\n";
$nlimit = 10_000_000;
$ilimit = prime_count($nlimit)-1;
$expect = 0; forprimes { $expect += $_ } $nlimit;

cmpthese($count,{
  'primes'    => sub { $s=0; $s += $_ for @{primes($nlimit)}; die unless $s == $expect; },
  'forprimes' => sub { $s=0; forprimes { $s += $_ } $nlimit;  die unless $s == $expect; },
  'pa index'  => sub { $s=0; tie my @primes, "Math::Prime::Util::PrimeArray";
                       $s += $primes[$_] for 0..$ilimit;
                       die unless $s == $expect; },
  'pa loop'   => sub { $s=0; tie my @primes, "Math::Prime::Util::PrimeArray";
                       for (@primes) { last if $_ > $nlimit; $s += $_; }
                       die $s unless $s == $expect; },
  'pa slice'  => sub { $s=0; tie my @primes, "Math::Prime::Util::PrimeArray";
                       $s += $_ for @primes[0..$ilimit];
                       die unless $s == $expect; },
  'pa each'   => sub { $s=0; tie my @primes, "Math::Prime::Util::PrimeArray";
                       while(my(undef,$v) = each @primes) { last if $v > $nlimit; $s += $v; }
                       die $s unless $s == $expect; },
  'pa shift'  => sub { $s=0; tie my @primes, "Math::Prime::Util::PrimeArray";
                       while ((my $p = shift @primes) <= $nlimit) { $s += $p; }
                       die unless $s == $expect; },
  'numseq'    => sub { $s=0; my $seq = Math::NumSeq::Primes->new;
                       while (1) { my($undev,$v) = $seq->next; last if $v > $nlimit; $s += $v; }
                       die $s unless $s == $expect; },
});
}

if (1) {
print '-' x 79, "\n";
print "Walk primes backwards from 1M\n";
print "Note: MPTA takes 4x longer than just calling MPU's nth_prime!\n";
$nlimit = 1_000_000;
$ilimit = prime_count($nlimit)-1;
$expect = 0; forprimes { $expect += $_ } $nlimit;

cmpthese($count,{
  'rev primes'=> sub { $s=0; $s += $_ for reverse @{primes($nlimit)}; die unless $s == $expect; },
  'nthprime'  => sub { $s=0; $s += nth_prime($_) for reverse 1..$ilimit+1; die unless $s == $expect; },
  'pa index'  => sub { $s=0; tie my @primes, "Math::Prime::Util::PrimeArray";
                       $s += $primes[$_] for reverse 0..$ilimit;
                       die unless $s == $expect; },
  'OO iter'   => sub { $s=0; my $it = prime_iterator_object($nlimit);
                       $s += $it->prev->value() for 0..$ilimit;
                       die unless $s == $expect; },
  'tiedarray' => sub { $s=0; tie my @primes, "Math::Prime::TiedArray", extend_step => 1000;
                       $s += $primes[$_] for reverse 0..$ilimit;
                       die unless $s == $expect; },
});
}

if (1) {
print '-' x 79, "\n";
print "Random walk in 1M\n";
print "MPTA takes about 2 minutes and lots of RAM per iteration.\n";
srand(29);
my @rindex;
do { push @rindex, int(rand(1000000)) } for 1..10000;
$expect = 0; $expect += nth_prime($_+1) for @rindex;

cmpthese($count,{
  'nthprime'  => sub { $s=0; $s += nth_prime($_+1) for @rindex; },
  'pa index'  => sub { $s=0; tie my @primes, "Math::Prime::Util::PrimeArray";
                       $s += $primes[$_] for @rindex;
                       die unless $s == $expect; },
   # Argh!  Is it possible to write a slower sieve than the one MPTA uses?
  #'tiedarray' => sub { $s=0; tie my @primes, "Math::Prime::TiedArray", extend_step => 10000;
  #                     $s += $primes[$_] for @rindex;
  #                     die unless $s == $expect; },
});
}

print '-' x 79, "\n";
