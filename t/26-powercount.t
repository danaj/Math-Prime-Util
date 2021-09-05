#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/perfect_power_count/;

my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};


my @A069623 = (1, 1, 1, 2, 2, 2, 2, 3, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12);
my @A070428 = (1, 4, 13, 41, 125, 367, 1111, 3395, 10491, 32670, 102231, 320990, 1010196, 3184138, 10046921, 31723592, 100216745, 316694005, 1001003332, 3164437425, 10004650118, 31632790244, 100021566157, 316274216762, 1000100055684);

$#A069623 = 40;
$#A070428 = 10;

plan tests => 0
            + 4    # perfect_power_count  basic tests
            + 1    # perfect_power_count  large value
            + 2    # perfect_power_count  ranges
            + 0;

is(perfect_power_count(0), 0, "perfect_power_count(0) = 0");
is(perfect_power_count(1), 1, "perfect_power_count(1) = 1");
is_deeply( [map { perfect_power_count(1+$_) } 0..$#A069623], \@A069623,  "perfect_power_count(n) for 1..".scalar(@A069623) );
is_deeply( [map { perfect_power_count(10**$_) } 0..$#A070428], \@A070428,  "perfect_power_count(10^n) for 0..$#A070428" );

# mpu 'say 1+vecsum(map{!!is_power($_)}1..12345678)'
is(perfect_power_count(12345678), 3762, "perfect_power_count(12345678) = 3762");


is( perfect_power_count(123456, 133332), 17, "perfect_power_count(123456,133332) = 17" );
is_deeply( [map { perfect_power_count($_,16) } 8,9,10],
           [3,2,1],
           "perfect_power_count(8..10,16) = 3,2,1" );
