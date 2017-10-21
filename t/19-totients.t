#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/euler_phi jordan_totient carmichael_lambda
                         divisor_sum moebius/;

#my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
#my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
#my $usegmp= Math::Prime::Util::prime_get_config->{'gmp'};
my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
$use64 = 0 if $use64 && 18446744073709550592 == ~0;

my %totients = (
     123456 => 41088,
     123457 => 123456,
  123456789 => 82260072,
);
my @A000010 = (0,1,1,2,2,4,2,6,4,6,4,10,4,12,6,8,8,16,6,18,8,12,10,22,8,20,12,18,12,28,8,30,16,20,16,24,12,36,18,24,16,40,12,42,20,24,22,46,16,42,20,32,24,52,18,40,24,36,28,58,16,60,30,36,32,48,20,66,32,44);
#@totients{0..$#A000010} = @A000010;

my @A002322 = (0,1,1,2,2,4,2,6,2,6,4,10,2,12,6,4,4,16,6,18,4,6,10,22,2,20,12,18,6,28,4,30,8,10,16,12,6,36,18,12,4,40,6,42,10,12,22,46,4,42,20,16,12,52,18,20,6,18,28,58,4,60,30,6,16,12,10,66,16,22,12,70,6,72,36,20,18,30,12,78,4,54,40,82,6,16,42,28,10,88,12,12,22,30,46,36,8,96,42,30,20,100,16,102,12,12,52,106,18,108,20,36,12,112,18,44,28,12,58,48,4,110,60,40,30,100,6,126,32,42,12,130,10,18,66,36,16,136,22,138,12,46,70,60,12,28,72,42,36,148,20,150,18,48,30,60,12,156,78,52,8,66,54,162,40,20,82,166,6,156,16,18,42,172,28,60,20,58,88,178,12,180,12,60,22,36,30,80,46,18,36,190,16,192,96,12,42,196,30,198,20);

plan tests => 2 + 9 + scalar(keys %totients)
                + 1 # Small Carmichael Lambda
                ;

###### euler_phi (totient)
{
  my @phi = map { euler_phi($_) } (0 .. $#A000010);
  is_deeply( \@phi, \@A000010, "euler_phi 0 .. $#A000010" );
}
{
  my @phi = euler_phi(0, $#A000010);
  is_deeply( \@phi, \@A000010, "euler_phi with range: 0, $#A000010" );
}
{
  my $s = 0;
  $s += $_ for euler_phi(1, 240);
  is($s, 17544, "sum of totients to 240");
}
while (my($n, $phi) = each (%totients)) {
  is( euler_phi($n), $phi, "euler_phi($n) == $phi" );
}
is_deeply( [euler_phi(0,0)], [0],     "euler_phi(0,0)" );
is_deeply( [euler_phi(1,0)], [],      "euler_phi with end < start" );
is_deeply( [euler_phi(0,1)], [0,1],   "euler_phi 0-1" );
is_deeply( [euler_phi(1,2)], [1,1],   "euler_phi 1-2" );
is_deeply( [euler_phi(1,3)], [1,1,2], "euler_phi 1-3" );
is_deeply( [euler_phi(2,3)], [1,2],   "euler_phi 2-3" );
is_deeply( [euler_phi(10,20)], [4,10,4,12,6,8,8,16,6,18,8], "euler_phi 10-20" );
is_deeply( [euler_phi(1513,1537)],
   [qw/1408 756 800 756 1440 440 1260 576 936 760 1522 504 1200 648
       1016 760 1380 384 1530 764 864 696 1224 512 1456/],
           "euler_phi(1513,1537)" );

###### Carmichael Lambda
{
  my @lambda = map { carmichael_lambda($_) } (0 .. $#A002322);
  is_deeply( \@lambda, \@A002322, "carmichael_lambda with range: 0, $#A000010" );
}
