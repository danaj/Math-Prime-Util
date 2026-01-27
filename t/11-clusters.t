#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/sieve_prime_cluster is_prime primes twin_primes addint subint powint/;
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $usegmp = Math::Prime::Util::prime_get_config->{'gmp'};
my $usexs = Math::Prime::Util::prime_get_config->{'xs'};

my @tests = (
  [ "A001359", [0, 2], [0,200], [3, 5, 11, 17, 29, 41, 59, 71, 101, 107, 137, 149, 179, 191, 197] ],
  [ "A022004", [0,2,6], [317321,319727], [qw/317321 317771 317957 318677 319127 319727/] ],
  [ "A022005", [0,4,6], [557857,560293], [qw/557857 558787 559213 560233 560293/] ],
);

my @patterns = (
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
);
my @high_check2 = (
  [ "9999956467211", 'A022009', [2,6,8,12,18,20] ],
  [ "9996858589169", 'A022010', [2,8,12,14,18,20] ],
  [ "99996813484481", 'A022011', [2,6,8,12,18,20,26] ],
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
push @high_check, @high_check2 if $extra;

#[2,6,8,18,20,30,32,36,38);   # Federighi
#[2,6,8,12,18,20,26,30,32,36,42,48,50,56);   # A257304
#[4,6,10,16,18,24,28,30,34,40,46,48,54,58,60,66);   # A257375
#[6,12,16,18,22,28,30,36,40,42,46,48);   # A214947

plan tests => scalar(@tests) + 2 + 2 + 2 * scalar(@patterns) + scalar(@high_check);

for my $t (@tests) {
  my($what, $tuple, $range, $expect) = @$t;
  shift @$tuple if $tuple->[0] == 0;
  my @res = sieve_prime_cluster($range->[0],$range->[1], @$tuple );
  is_deeply( \@res, $expect, "$what @$range" );
}

is_deeply( [sieve_prime_cluster(1,1e9,2,4)], [3], "Inadmissible pattern (0,2,4) finds (3,5,7)");
is_deeply( [sieve_prime_cluster(1,1e9,2,8,14,26)], [3,5], "Inadmissible pattern (0,2,8,14,26) finds (3,5,11,17,29) and (5,7,13,19,31)");

my @pcache;  # holds primes in two ranges

my($sbeg,$send) = (0, 100000);
$send += 1000000 if $extra;
$pcache[0] = primes($sbeg,$send+256);

my $mbeg = powint(10,21);
my $mend = $mbeg + 10000 + int(rand(100000));
$mend += 100000 if $extra;
if ($usegmp) {
  $pcache[1] = primes($mbeg,$mend+256);
} else {
  # Without GMP and using the Calc backend, this is just painful slow
  $mend = $mbeg + 5000;
  $pcache[1] = [map { $mbeg+$_ } (qw/117 193 213 217 289 327 367 373 399 409 411 427 433 447 471 553 609 723 733 951 1063 1081 1213 1237 1311 1383 1411 1417 1459 1521 1573 1581 1687 1731 1749 1867 1897 2001 2011 2041 2049 2203 2209 2257 2259 2307 2317 2343 2349 2583 2611 2673 2701 2713 2719 2761 2803 2823 2961 3007 3021 3271 3289 3327 3331 3369 3399 3423 3483 3657 3759 3777 3861 3897 3973 3999 4011 4017 4039 4063 4081 4119 4123 4197 4231 4297 4353 4359 4381 4437 4521 4581 4591 4671 4743 4749 4791 4813 4851 4891 4897 4977 5203 5277 5317 5371 5427 5437 5499 5577 5683 5719 5751 5763 5913 5959 6003 6009 6103 6247 6297 6309 6493 6531 6727 6747 6759 6781 6783 6853 6871 6883 6993 7039 7059 7069 7147 7231 7269 7413 7467 7471 7509 7527 7639 7681 7689 7711 7741 7761 7887 8011 8071 8143 8173 8187 8221 8223 8283 8299 8343 8407 8467 8497 8587 8623 8761 8799 8973 9069 9111 9121 9159 9183 9187 9211 9217 9271 9333 9349 9369 9477 9501 9723 9847 9861 9961 9999/)];
}

###### twin primes native
{
  my $beg = 0;
  my $end = $extra ? $send : 20000;
  my @sieve = sieve_prime_cluster($beg,$end,2);
  my @tuple = ktuple($beg,$end,$pcache[0],2);
  my $num = scalar(@tuple);
  is_deeply( \@sieve, \@tuple, "Pattern [2] $num in range $beg .. $end");
}

###### twin primes bigint
{
  my $beg = $mbeg;
  my $end = $extra ? $mend : $mbeg + 1000;
  my @sieve = sieve_prime_cluster($beg,$end,2);
  my @tuple = ktuple($beg,$end,$pcache[1],2);
  my $num = scalar(@tuple);
  is_deeply( \@sieve, \@tuple, "Pattern [2] $num in range $beg .. $end");
}

###### extended patterns native
for my $pat (@patterns) {
  my @pat = @$pat;
  shift @pat if $pat[0] == 0;
  my @sieve = sieve_prime_cluster($sbeg,$send,@pat);
  my @tuple = ktuple($sbeg,$send,$pcache[0],@pat);
  my $num = scalar(@tuple);

  is_deeply( \@sieve, \@tuple, "Pattern [@pat] $num in range $sbeg .. $send");
}

###### extended patterns bigint
for my $pat (@patterns) {
  my @pat = @$pat;
  shift @pat if $pat[0] == 0;
  my @sieve = sieve_prime_cluster($mbeg,$mend,@pat);
  my @tuple = ktuple($mbeg,$mend,$pcache[1],@pat);
  my $num = scalar(@tuple);

  is_deeply( \@sieve, \@tuple, "Pattern [@pat] $num in range $mbeg .. $mend");
}

####### target a small window around known large clusters
for my $test (@high_check) {
  my($n,$name,$cl) = @$test;
  my $window = ($usexs && $usegmp) ? 1000000 : 1000;
  my @res = sieve_prime_cluster(subint($n,$window), addint($n,$window), @$cl);
  is_deeply(\@res, [$n], "Window around $name high cluster finds the cluster");
}

sub ktuple {
  my($beg, $end, $pcache, @pat) = @_;

  my @p = grep { $_ >= $beg && $_ <= $end } @$pcache;
  my %prhash = map { $_ => 1; } @$pcache;

  foreach my $c (@pat) {
    @p = grep { $prhash{$_+$c} } @p;
  }
  @p;
}
