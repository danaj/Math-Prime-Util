#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/abundance aliquot_sum divisor_sum/;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
#my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
#my $usegmp= Math::Prime::Util::prime_get_config->{'gmp'};
#my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
#$use64 = 0 if $use64 && 18446744073709550592 == ~0;

my %sigmak = (
  # A0000005
  0 => [1,2,2,3,2,4,2,4,3,4,2,6,2,4,4,5,2,6,2,6,4,4,2,8,3,4,4,6,2,8,2,6,4,4,4,9,2,4,4,8,2,8,2,6,6,4,2,10,3,6,4,6,2,8,4,8,4,4,2,12,2,4,6,7,4,8,2,6,4,8,2,12,2,4,6,6,4,8,2,10,5,4,2,12,4,4,4,8,2,12,4,6,4,4,4,12,2,6,6,9,2,8,2,8],
  # A000203
  1 => [1, 3, 4, 7, 6, 12, 8, 15, 13, 18, 12, 28, 14, 24, 24, 31, 18, 39, 20, 42, 32, 36, 24, 60, 31, 42, 40, 56, 30, 72, 32, 63, 48, 54, 48, 91, 38, 60, 56, 90, 42, 96, 44, 84, 78, 72, 48, 124, 57, 93, 72, 98, 54, 120, 72, 120, 80, 90, 60, 168, 62, 96, 104, 127, 84, 144, 68, 126, 96, 144],
  # A001157
  2 => [1, 5, 10, 21, 26, 50, 50, 85, 91, 130, 122, 210, 170, 250, 260, 341, 290, 455, 362, 546, 500, 610, 530, 850, 651, 850, 820, 1050, 842, 1300, 962, 1365, 1220, 1450, 1300, 1911, 1370, 1810, 1700, 2210, 1682, 2500, 1850, 2562, 2366, 2650, 2210, 3410, 2451, 3255],
  # A001158
  3 => [1, 9, 28, 73, 126, 252, 344, 585, 757, 1134, 1332, 2044, 2198, 3096, 3528, 4681, 4914, 6813, 6860, 9198, 9632, 11988, 12168, 16380, 15751, 19782, 20440, 25112, 24390, 31752, 29792, 37449, 37296, 44226, 43344, 55261, 50654, 61740, 61544],
);

plan tests => 1     # aliquot_sum(n)
            + 1     # abundance(n)
            + 1     # divisor_sum(0,k)
            + 1     # divisor_sum(n)
            + 1     # divisor_sum(n,k)
            + 1     # divisor_sum(n,sub{})
            + 1     # tau
            + 1;    # overflow

subtest 'aliquot_sum' => sub {
  # A001065 starting at n=1
  is_deeply( [map { aliquot_sum($_) } 0..20],
             [0, 0,1,1,3,1,6,1,7,4,8,1,16,1,10,9,15,1,21,1,22],
             "aliquot_sum(0..20)" );
  # Perfect numbers: aliquot_sum(n) == n
  is_deeply( [map { aliquot_sum($_) } (6, 28, 496, 8128)],
             [6, 28, 496, 8128], "perfect numbers" );
  # Amicable pairs: each is the aliquot_sum of the other
  is( aliquot_sum(220), 284, "aliquot_sum(220) = 284" );
  is( aliquot_sum(284), 220, "aliquot_sum(284) = 220" );
  # Prime: aliquot_sum(p) = 1
  is( aliquot_sum(999999937), 1, "aliquot_sum of prime = 1" );
  # Prime power: aliquot_sum(2^k) = 2^k - 1
  is( "".aliquot_sum("4611686018427387904"), "4611686018427387903",
      "aliquot_sum(2^62) = 2^62-1" );
  # Large n where sigma overflows UV, dispatches to PP
  is( "".aliquot_sum("9600000000000000688"),
      "9000000000000000676",
      "aliquot_sum where sigma overflows UV" );
  # simple bigint
  is( "".aliquot_sum("927208363107752634625925"),"317454036909046584758395",
      "aliquot_sum(927208363107752634625925) = 317454036909046584758395" );
};


subtest 'abundance' => sub {
  # abundance(n) = sigma(n) - 2n: positive=abundant, zero=perfect, negative=deficient
  is_deeply( [map { abundance($_) } 1..20],
             [-1,-1,-2,-1,-4,0,-6,-1,-5,-2,-10,4,-12,-4,-6,-1,-16,3,-18,2],
             "abundance(1..20)" );
  # Perfect numbers: abundance = 0  (A000396)
  is_deeply( [map { abundance($_) } (6, 28, 496, 8128)],
             [0, 0, 0, 0], "perfect numbers have abundance 0" );
  # First odd abundant number (A005231): 945 = 3^3 * 5 * 7
  is( abundance(945), 30, "abundance(945) = 30 (first odd abundant)" );
  # For a prime p, sigma(p) = p+1, so abundance(p) = 1-p
  is( abundance(999999937), -999999936, "abundance of large prime = 1-p" );
  # bigint
  is("".abundance("927208363107752634625925"),"-609754326198706049867530", "abundance(927208363107752634625925)");
};

###### Divisor sum
is_deeply( [map { divisor_sum(0,$_) } 0..5], [0,0,0,0,0,0], "divisor_sum(0,k) = 0" );

subtest 'divisor_sum(n)', sub {
  my @slist = map { divisor_sum($_) } 1 .. scalar @{$sigmak{1}};
  is_deeply(\@slist, $sigmak{1}, "divisor_sum(n)");
};

subtest 'divisor_sum(n,k)', sub {
  while (my($k, $sigmaref) = each (%sigmak)) {
    my @slist = map { divisor_sum($_,$k) } 1 .. scalar @$sigmaref;
    is_deeply( \@slist, $sigmaref, "Sigma_$k using divisor_sum(n,k)" );
  }
};
subtest 'divisor_sum(n,sub{})', sub {
  while (my($k, $sigmaref) = each (%sigmak)) {
    my @slist = map {divisor_sum($_,sub{int($_[0])**$k})} 1..scalar @$sigmaref;
    is_deeply( \@slist, $sigmaref, "Sum of divisors to the ${k}th power: Sigma_$k" );
  }
};

subtest 'tau using divisor_sum', sub {
  my $len = scalar @{$sigmak{0}};
  my @slist1 = map { divisor_sum($_, sub {1}) } 1 .. $len;
  my @slist2 = map { divisor_sum($_, 0      ) } 1 .. $len;
  is_deeply( \@slist1, $sigmak{0}, "tau as divisor_sum(n, sub {1})" );
  is_deeply( \@slist2, $sigmak{0}, "tau as divisor_sum(n, 0)" );

  # tau_4 A007426
  my @tau4 = (1,4,4,10,4,16,4,20,10,16,4,40,4,16,16,35,4,40,4,40,16,16,4,80,10,16,20,40,4,64,4,56,16,16,16,100);
  push @tau4, (4,16,16,80,4,64,4,40,40,16,4,140,10,40,16,40,4,80,16,80,16,16,4,160,4,16,40,84,16,64,4,40,16,64,4,200,4,16,40,40,16) if $extra;
  my @t;
  foreach my $n (1 .. scalar @tau4) {
    push @t, divisor_sum($n, sub { divisor_sum($_[0],sub { divisor_sum($_[0],0) }) });
  }
  is_deeply( \@t, \@tau4, "Tau4 (A007426), nested divisor sums" );
};

subtest 'integer overflow check', sub {
  # Check some cases for integer overflow
  is( "".divisor_sum( 1<<27,2),   "24019198012642645", "divisor_sum(2^27,2)" );
  is( "".divisor_sum(262144,3),   "20587884010836553", "divisor_sum(2^18,3)" );
  is( "".divisor_sum( 16384,4),   "76861433640456465", "divisor_sum(2^16,4)" );
  is( "".divisor_sum(  2048,5),   "37191016277640225", "divisor_sum(2^11,5)" );
  is( "".divisor_sum(  5003,5), "3134386256752025244", "divisor_sum(5003,5)" );
};
