#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/perfect_power_count prime_power_count/;

my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};


my @A069623 = (1, 1, 1, 2, 2, 2, 2, 3, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12);
my @A070428 = (1, 4, 13, 41, 125, 367, 1111, 3395, 10491, 32670, 102231, 320990, 1010196, 3184138, 10046921, 31723592, 100216745, 316694005, 1001003332, 3164437425, 10004650118, 31632790244, 100021566157, 316274216762, 1000100055684);
my @A025528 = (0, 1, 2, 3, 4, 4, 5, 6, 7, 7, 8, 8, 9, 9, 9, 10, 11, 11, 12, 12, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 18, 18, 18, 18, 18, 19, 19, 19, 19, 20, 20, 21, 21, 21, 21, 22, 22, 23, 23, 23, 23, 24, 24, 24, 24, 24, 24, 25, 25, 26, 26, 26, 27, 27, 27, 28, 28, 28, 28, 29, 29, 30, 30);
my @A267712 = (7, 35, 193, 1280, 9700, 78734, 665134, 5762859, 50851223, 455062595, 4118082969, 37607992088, 346065767406, 3204942420923, 29844572385358, 279238346816392, 2623557174778438, 24739954338671299, 234057667428388198, 2220819603016308079);

$#A069623 = 40;
$#A070428 = 10;
$#A025528 = 40;
$#A267712 = ($usexs || $extra) ? 7 : 6;

plan tests => 0
            + 4
            + 4
            + 2
            + 0;

is(perfect_power_count(0), 0, "perfect_power_count(0) = 0");
is(perfect_power_count(1), 1, "perfect_power_count(1) = 1");
is_deeply( [map { perfect_power_count(1+$_) } 0..$#A069623], \@A069623,  "perfect_power_count(n) for 1..".scalar(@A069623) );
is_deeply( [map { perfect_power_count(10**$_) } 0..$#A070428], \@A070428,  "perfect_power_count(10^n) for 0..$#A070428" );

is(prime_power_count(0), 0, "prime_power_count(0) = 0");
is(prime_power_count(1), 0, "prime_power_count(1) = 0");
is_deeply( [map { prime_power_count(1+$_) } 0..$#A025528], \@A025528,  "prime_power_count(n) for 1..".scalar(@A025528) );
is_deeply( [map { prime_power_count(10**(1+$_)) } 0..$#A267712], \@A267712,  "prime_power_count(10^n) for 1..".scalar(@A267712) );

# mpu 'say 1+vecsum(map{!!is_power($_)}1..12345678)'
is(perfect_power_count(12345678), 3762, "perfect_power_count(12345678) = 3762");
# mpu 'say vecsum(map{!!is_prime_power($_)}1..12345678)'
is(prime_power_count(12345678), 809830, "prime_power_count(12345678) = 809830");
