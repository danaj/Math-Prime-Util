#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
#  TODO: nth_prime_power
use Math::Prime::Util qw/prime_powers prime_power_count
                         next_prime_power prev_prime_power
                         is_prime_power/;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $usexs = Math::Prime::Util::prime_get_config->{'xs'};

my @small = (2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25, 27, 29, 31, 32, 37, 41, 43, 47, 49, 53, 59, 61, 64, 67, 71, 73, 79, 81, 83, 89, 97);

plan tests =>   3   # prime_powers
              + 3   # next_prime_power
              + 4   # prev_prime_power
              + 0   # prime_power_count
              + 0   # nth_prime_power
                ;



###### prime_powers

is_deeply( prime_powers(998), [grep { is_prime_power($_) } 1..998], "prime_powers(998)" );

is_deeply( [map { prime_powers($_) } 1..100],
           [map { my $n=$_; [grep { $_ <= $n } @small] } 1..100],
           "prime_powers(1..100)" );

my $base1 = 1441897;
is_deeply( prime_powers($base1,$base1+1000),
           [map { $base1+$_ } (34,36,52,66,84,106,112,120,156,160,172,174,190,246,262,276,294,312,330,354,370,382,402,420,424,430,436,444,454,480,496,504,514,532,540,556,562,612,616,630,634,652,682,694,702,714,724,730,736,744,756,760,772,820,826,834,846,886,900,930,952,966,972,976,990,1000)],
           "prime_powers($base1, $base1 + 1000)" );


###### next_prime_power

is_deeply( [map { next_prime_power($_-1) } @small],
           \@small,
           "next_prime_power");

is_deeply( [map { next_prime_power(2**$_-1) } 1..9],
           [2, 4, 8, 16, 32, 64, 128, 256, 512],
           "next_prime_power(2^i+1)" );
is_deeply( [map { next_prime_power(2**$_) } 0..9],
           [2, 3, 5, 9, 17, 37, 67, 131, 257, 521],
           "next_prime_power(2^i)" );

###### prev_prime_power

is_deeply( [map { prev_prime_power($_) } 0..2], [undef,undef,undef], "prev_prime_power(0..2) = undef" );

is_deeply( [map { prev_prime_power($_+1) } @small],
           \@small,
           "prev_prime_power");

is_deeply( [map { prev_prime_power(2**$_+1) } 1..9],
           [2, 4, 8, 16, 32, 64, 128, 256, 512],
           "prev_prime_power(2^i+1)" );
is_deeply( [map { prev_prime_power(2**$_) } 1..9],
           [undef, 3, 7, 13, 31, 61, 127, 251, 509],
           "prev_prime_power(2^i)" );

###### prime_power_count

# See t/26-powercount.t

###### nth_prime_power

# TODO add function
