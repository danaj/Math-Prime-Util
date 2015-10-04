#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/sieve_prime_cluster is_prime primes twin_primes/;
use Math::BigInt try => "GMP,Pari";
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};

my @tests = (
  [ "A001359", [0, 2], [0,200], [3, 5, 11, 17, 29, 41, 59, 71, 101, 107, 137, 149, 179, 191, 197] ],
  [ "A022004", [0,2,6], [317321,319727], [qw/317321 317771 317957 318677 319127 319727/] ],
  [ "A022005", [0,4,6], [557857,560293], [qw/557857 558787 559213 560233 560293/] ],
);

my @patterns = (
  [0,2],
  [0,2,6],
  [0,4,6],
  [0,2,6,8],
  [0,2,6,8,12],
  [0,4,6,10,12],
  [0,4,6,10,12,16],
  [0,2,8,12,14,18,20],
  [0,2,6,8,12,18,20],
);

my @high_check = (
  [ "999999217031", 'A022006', [2,6,8,12] ],
  [ "999998356957", 'A022007', [4,6,10,12] ],
  [ "999961920817", 'A022008', [4,6,10,12,16] ],
  [ "9999956467211", 'A022009', [2,6,8,12,18,20] ],
  [ "9996858589169", 'A022010', [2,8,12,14,18,20] ],
  [ "99996813484481", 'A022010', [2,6,8,12,18,20,26] ],
  [ "99997194198047", 'A022012', [2,6,12,14,20,24,26] ],
  [ "99996215495153", 'A022013', [6,8,14,18,20,24,26] ],
  [ "999897629673401", 'A022545', [2,6,8,12,18,20,26,30] ],
  [ "506946970236647", 'A022546', [2,6,12,14,20,24,26,30] ],
  [ "291985604331973", 'A022547', [4,6,10,16,18,24,28,30] ],
  [ "999823346788939", 'A022548', [4,10,12,18,22,24,28,30] ],
  [ "29997979809623711", 'A027569', [2,6,8,12,18,20,26,30,32] ],
  [ "29998892234668517", 'A027570', [2,6,12,14,20,24,26,30,32] ],
  [ "9996248338169127877", 'A213601', [6,10,12,16,22,24,30,34,36,40,42] ],
  [ "2830868185774228331", 'A213645', [2,6,8,12,18,20,26,30,32,36,42] ],
  [ "999955337060684083", 'A213646', [4,6,10,16,18,24,28,30,34,36] ],
  [ "999930334493085881", 'A213647', [2,6,8,12,18,20,26,30,32,36] ],
);

#[2,6,8,18,20,30,32,36,38);   # Federighi
#[2,6,8,12,18,20,26,30,32,36,42,48,50,56);   # A257304
#[4,6,10,16,18,24,28,30,34,40,46,48,54,58,60,66);   # A257375
#[6,12,16,18,22,28,30,36,40,42,46,48);   # A214947

plan tests => scalar(@tests) + 2 + 2 * scalar(@patterns) + scalar(@high_check);

for my $t (@tests) {
  my($what, $tuple, $range, $expect) = @$t;
  shift @$tuple if $tuple->[0] == 0;
  my @res = sieve_prime_cluster($range->[0],$range->[1], @$tuple );
  is_deeply( \@res, $expect, "$what @$range" );
}

is_deeply( [sieve_prime_cluster(1,1e10,2,4)], [3], "Inadmissible pattern (0,2,4) finds (3,5,7)");
is_deeply( [sieve_prime_cluster(1,1e10,2,8,14,26)], [3,5], "Inadmissible pattern (0,2,8,14,26) finds (3,5,11,17,29) and (5,7,13,19,31)");

my($sbeg,$send) = (0, 100000);
my $mbeg = Math::BigInt->new(10)**21;
my $mend = $mbeg + 10000 + int(rand(100000));

$send += 1000000 if $extra;
$mend += 100000 if $extra;

my $small = [ primes($sbeg,$send), twin_primes($sbeg,$send) ];
my $large = [ primes($mbeg,$mend), twin_primes($mbeg,$mend) ];

for my $pat (@patterns) {
  my @pat = @$pat;
  shift @pat if $pat[0] == 0;
  my @sieve = sieve_prime_cluster($sbeg,$send,@pat);
  my @tuple = ktuple($sbeg,$send,$small,@pat);
  my $num = scalar(@tuple);

  is_deeply( \@sieve, \@tuple, "Pattern [@pat] $num in range $sbeg .. $send");
}
for my $pat (@patterns) {
  my @pat = @$pat;
  shift @pat if $pat[0] == 0;
  my @sieve = sieve_prime_cluster($mbeg,$mend,@pat);
  my @tuple = ktuple($mbeg,$mend,$large,@pat);
  my $num = scalar(@tuple);

  is_deeply( \@sieve, \@tuple, "Pattern [@pat] $num in range $mbeg .. $mend");
}

for my $test (@high_check) {
  my($n,$name,$cl) = @$test;
  my @res = sieve_prime_cluster($n-1e6, $n+1e6, @$cl);
  is_deeply(\@res, [$n], "Window around $name high cluster finds the cluster");
}

sub ktuple {
  my($beg, $end, $prset, @pat) = @_;
  my @p;
  if (@pat && $pat[0] == 2) {
    @p = @{$prset->[1]};
    shift @pat;
  } else {
    @p = @{$prset->[0]};
  }
  for my $c (@pat) {
    @p = grep { is_prime($_+$c) } @p;
  }
  shift @p while @p && $p[0] < $beg;
  pop @p while @p && $p[-1] > $end;
  @p;
}
