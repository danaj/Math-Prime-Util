#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/euler_phi jordan_totient carmichael_lambda
                         divisor_sum moebius inverse_totient
                         sumtotient dedekind_psi/;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
#my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
#my $usegmp= Math::Prime::Util::prime_get_config->{'gmp'};
my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
$use64 = 0 if $use64 && 18446744073709550592 == ~0;


my @A002322 = (0,1,1,2,2,4,2,6,2,6,4,10,2,12,6,4,4,16,6,18,4,6,10,22,2,20,12,18,6,28,4,30,8,10,16,12,6,36,18,12,4,40,6,42,10,12,22,46,4,42,20,16,12,52,18,20,6,18,28,58,4,60,30,6,16,12,10,66,16,22,12,70,6,72,36,20,18,30,12,78,4,54,40,82,6,16,42,28,10,88,12,12,22,30,46,36,8,96,42,30,20,100,16,102,12,12,52,106,18,108,20,36,12,112,18,44,28,12,58,48,4,110,60,40,30,100,6,126,32,42,12,130,10,18,66,36,16,136,22,138,12,46,70,60,12,28,72,42,36,148,20,150,18,48,30,60,12,156,78,52,8,66,54,162,40,20,82,166,6,156,16,18,42,172,28,60,20,58,88,178,12,180,12,60,22,36,30,80,46,18,36,190,16,192,96,12,42,196,30,198,20);

plan tests => 1     # euler_phi
            + 1     # carmichael_lambda
            + 1     # inverse_totient
            + 1     # sumtotient
            + 1;    # dedekind_psi

subtest 'Totient (euler_phi)', sub {
  my @A000010 = (0,1,1,2,2,4,2,6,4,6,4,10,4,12,6,8,8,16,6,18,8,12,10,22,8,20,12,18,12,28,8,30,16,20,16,24,12,36,18,24,16,40,12,42,20,24,22,46,16,42,20,32,24,52,18,40,24,36,28,58,16,60,30,36,32,48,20,66,32,44);
  my @phi1 = map { euler_phi($_) } (0 .. $#A000010);
  is_deeply( \@phi1, \@A000010, "euler_phi 0 .. $#A000010" );
  my @phir = euler_phi(0, $#A000010);
  is_deeply( \@phir, \@A000010, "euler_phi with range: 0, $#A000010" );
  {
    my $s = 0;
    $s += $_ for euler_phi(1, 240);
    is($s, 17544, "sum of totients to 240");
  }
  is(euler_phi(-123456), 0, "euler_phi(-n) = 0");
  is(euler_phi(123456), 41088, "euler_phi(123456) = 41088");
  is(euler_phi(123457), 123456, "euler_phi(123457) = 123456");
  is(euler_phi(123456789), 82260072, "euler_phi(123456789) = 82260072");

  is_deeply([euler_phi(0,0)], [0],     "euler_phi(0,0)");
  is_deeply([euler_phi(1,0)], [],      "euler_phi with end < start");
  is_deeply([euler_phi(0,1)], [0,1],   "euler_phi 0-1");
  is_deeply([euler_phi(1,2)], [1,1],   "euler_phi 1-2");
  is_deeply([euler_phi(1,3)], [1,1,2], "euler_phi 1-3");
  is_deeply([euler_phi(2,3)], [1,2],   "euler_phi 2-3");
  is_deeply([euler_phi(10,20)], [4,10,4,12,6,8,8,16,6,18,8], "euler_phi 10-20");
  is_deeply( [euler_phi(1513,1537)],
    [qw/1408 756 800 756 1440 440 1260 576 936 760 1522 504 1200 648
        1016 760 1380 384 1530 764 864 696 1224 512 1456/],
    "euler_phi(1513,1537)" );
  # negative euler_phi returns zero
  is_deeply( [euler_phi(-5,5)], [0,0,0,0,0,0,1,1,2,2,4], "euler_phi -5 to 5" );

  is_deeply([[euler_phi(4294967293,4294967295)],
             [euler_phi(4294967293,4294967296)],
             [euler_phi(4294967295,4294967297)],
             [euler_phi(4294967296,4294967298)]],
            [[4294493280,2147483646,2147483648],[4294493280,2147483646,2147483648,2147483648],[2147483648,2147483648,4288266240],[2147483648,4288266240,1431655764]],
            "euler_phi ranges around 2^32");

  SKIP: {
    skip "ranges around 2^64 only on 64-bit",1 unless $use64;
    is_deeply([euler_phi("18446744073709551613","18446744073709551615"),
               euler_phi("18446744073709551613","18446744073709551616"),
               euler_phi("18446744073709551615","18446744073709551616")],
              [qw/17023385317621506048 7713001620195508224 9208981628670443520 17023385317621506048 7713001620195508224 9208981628670443520 9223372036854775808 9208981628670443520 9223372036854775808/],
              "euler_phi ranges around 2^64");
  }
};

{ ###### Carmichael Lambda
  my @lambda = map { carmichael_lambda($_) } (0 .. $#A002322);
  is_deeply( \@lambda, \@A002322, "carmichael_lambda 0 .. $#A002322" );
}

subtest 'Inverse Totient', sub {
  my $tot = 0;
  $tot += 0+inverse_totient($_) for 0..100;
  is($tot, 198, "Totient count 0-100 = 198");
  is(0+inverse_totient(1728), 62, "inverse_totient(1728) = 62");
  is(0+inverse_totient(362880), 1138, "inverse_totient(9!) = 1138");
  SKIP: {
    skip "Larger inverse totient with EXTENDED_TESTING",1 unless $extra;
    is(0+inverse_totient(3978374400), 63600, "inverse_totient(3978374400) = 63600");
  }

  is_deeply( [inverse_totient(0)], [], "inverse_totient(0)" );
  is_deeply( [inverse_totient(1)], [1,2], "inverse_totient(1)" );
  is_deeply( [inverse_totient(2)], [3,4,6], "inverse_totient(2)" );
  is_deeply( [inverse_totient(3)], [], "inverse_totient(3)" );
  is_deeply( [inverse_totient(4)], [5,8,10,12], "inverse_totient(4)" );
  is_deeply( [inverse_totient(2*12135413)], [], "inverse_totient(2*12135413)" );
  is_deeply( [inverse_totient(2*10754819)], [21509639,43019278], "inverse_totient(2*10754819)" );

  is_deeply( [inverse_totient(10000008)], [10555583,15000039,21111166,30000078], "inverse_totient(10000008)" );
  is_deeply( [inverse_totient(10000)], [10291,12625,13805,18825,20582,20625,22088,25000,25100,25250,27500,27610,33132,37500,37650,41250], "inverse_totient(10000)" );
  ok( scalar(grep { $_ == 123456789} inverse_totient(82260072)) == 1, "inverse_totient(82260072) includes 123456789" );
};

subtest 'sumtotient', sub {
  is_deeply([map { sumtotient($_) } 0..10], [0,1,2,4,6,10,12,18,22,28,32], "sumtotient(0..10)");
  is(sumtotient(12345),46326398,"sumtotient(12345)");
  is("".sumtotient(654321),130137945644,"sumtotient(654321)");
};

subtest 'dedekind_psi', sub {
  is( dedekind_psi(0),  0, "dedekind_psi(0) == 0" );
  is( dedekind_psi(-7), 0, "dedekind_psi(-7) == 0" );
  # A001615: psi(n) for n=1..20
  is_deeply( [map { dedekind_psi($_) } 1..20],
             [qw/1 3 4 6 6 12 8 12 12 18 12 24 14 24 24 24 18 36 20 36/],
             "dedekind_psi(1..20)" );
  is( dedekind_psi(30), 72, "dedekind_psi(30) = 72" );   # from pod example
  is( dedekind_psi(1000), 1800, "dedekind_psi(1000) = 1800" );
  is( dedekind_psi(8589934592), 12884901888, "dedekind_psi(2^33) = 3*2^32" );
  if ($use64) {
    is("".dedekind_psi("927208363107752634625925"), "1204512000016257309081600", "dedekind_psi(927208363107752634625925) = 1204512000016257309081600" );
  }
};
