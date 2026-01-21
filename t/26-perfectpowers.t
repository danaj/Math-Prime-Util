#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_perfect_power
                         next_perfect_power prev_perfect_power
                         perfect_power_count perfect_power_count_approx
                         perfect_power_count_lower perfect_power_count_upper
                         nth_perfect_power nth_perfect_power_approx
                         nth_perfect_power_lower nth_perfect_power_upper/;

my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
$use64 = 0 if $use64 && 18446744073709550592 == ~0;


my @A069623 = (1, 1, 1, 2, 2, 2, 2, 3, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12);
my @A070428 = (1, 4, 13, 41, 125, 367, 1111, 3395, 10491, 32670, 102231, 320990, 1010196, 3184138, 10046921, 31723592, 100216745, 316694005, 1001003332, 3164437425, 10004650118, 31632790244, 100021566157, 316274216762, 1000100055684);
my @A001597 = (1, 4, 8, 9, 16, 25, 27, 32, 36, 49, 64, 81, 100, 121, 125, 128, 144, 169, 196, 216, 225, 243, 256, 289, 324, 343, 361, 400, 441, 484, 512, 529, 576, 625, 676, 729, 784, 841, 900, 961, 1000, 1024, 1089, 1156, 1225, 1296, 1331, 1369, 1444, 1521, 1600, 1681, 1728, 1764);

$#A069623 = 40;
$#A070428 = 10;

my @pp100 = (qw/-64 -32 -27 -8 -1 0 1 4 8 9 16 25 27 32 36 49 64 81 100/);

my @uviv = ( [qw/-18446745128696702936 -18446724184312856125/],
             [qw/-9223372036854775808 -9223358842721533951/],
             [qw/-4298942376 -4291015625/],
             [qw/-2147483648 -2146689000/],
             [qw/2147395600 2147483648/],
             [qw/4294836225 4294967296/],
             [qw/9223372030926249001 9223372036854775808/],
             [qw/18446744065119617025 18446744073709551616/], );

plan tests => 0
            + 4    # is_perfect_power
            + 8    # next / prev
            + 4    # count  basic tests
            + 1    # count  large value
            + 2    # count  ranges
            + 3    # nth
            + 1    # small bounds
            + 6    # upper, lower, approx
            + 0;

######  is_perfect_power

is_deeply( [map { is_perfect_power($_) } 0..10], [1,1,0,0,1,0,0,0,1,1,0], "is_perfect_power(0 .. 10)" );
is_deeply( [grep { is_perfect_power($_) } -100..100], \@pp100, "is_perfect_power(-100 .. 100)" );
is( is_perfect_power("18446744065119617025"), 1, "is_perfect_power(18446744065119617025)" );
is( is_perfect_power("18446744073709551616"), 1, "is_perfect_power(18446744073709551616)" );


######  next / prev

is_deeply( [map { next_perfect_power($_) } 0..20],
           [1,4,4,4,8,8,8,8,9,16,16,16,16,16,16,16,25,25,25,25,25],
           "next perfect power with small inputs" );
is_deeply( [map { prev_perfect_power($_) } 0..20],
           [-1,0,1,1,1,4,4,4,4,8,9,9,9,9,9,9,9,16,16,16,16],
           "prev perfect power with small inputs" );

is_deeply( [map { next_perfect_power($_) } -9 .. 9],
           [-8,-1,-1,-1,-1,-1,-1,-1,0,1,4,4,4,8,8,8,8,9,16],
           "next perfect power with small inputs around zero" );
is_deeply( [map { prev_perfect_power($_) } -9 .. 9],
           [-27,-27,-8,-8,-8,-8,-8,-8,-8,-1,0,1,1,1,4,4,4,4,8],
           "prev perfect power with small inputs around zero" );

is_deeply( [map { next_perfect_power($_) } @pp100],
           [@pp100[1..$#pp100], 121],
           "next_perfect_power on perfect powers -100 to 100" );
is_deeply( [map { prev_perfect_power($_) } @pp100],
           [-125, @pp100[0..$#pp100-1]],
           "prev_perfect_power on perfect powers -100 to 100" );

{
  my(@gotprev, @expprev,  @gotnext, @expnext);
  for my $pair (@uviv) {
    push @expprev, $pair->[0];
    push @gotprev, prev_perfect_power($pair->[1]);
    push @expnext, $pair->[1];
    push @gotnext, next_perfect_power($pair->[0]);
    #is($gotprev[-1],$expprev[-1],"prev_perfect_power($pair->[1]) = $pair->[0]");
    #is($gotnext[-1],$expnext[-1],"next_perfect_power($pair->[0]) = $pair->[1]");
  }
  is_deeply( \@gotprev, \@expprev, "prev_perfect_power on numbers crossing 32-bit/64-bit boundaries" );
  is_deeply( \@gotnext, \@expnext, "next_perfect_power on numbers crossing 32-bit/64-bit boundaries" );
}

######  perfect_power_count
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

######  nth_perfect_power

is_deeply( [map { nth_perfect_power($_) } 1 .. scalar(@A001597)],
           \@A001597,
           "nth perfect_powers creates A001597" );

is_deeply( [map { nth_perfect_power($_) } 67224..67229],
           [qw/4294574089 4294705156 4294836225 4294967296 4295098369 4295229444/],
           "nth perfect powers with results around 2^32" );
SKIP: {
  skip "ranges around 2^64 only on 64-bit",1 unless $use64;
  is_deeply( [map { nth_perfect_power($_) } 4297615579,4297615580,4297615581,4297615582],
             [qw/18446744047939747849 18446744056529682436 18446744065119617025 18446744073709551616/],
             "nth perfect powers with results around 2^64" );
}

######  approx and bounds for count and nth

is_deeply( [map{cmp_rn($_+1,$A001597[$_])} 0..$#A001597], \@A001597, "small perfect power limits" );

is( cmp_rn(1571,2048383), 2048383, "perfect power limits for 1571" );
is( cmp_rn(59643,3373286400), 3373286400, "perfect power limits for 59643" );
is( cmp_rn(15964377,"252826822479841"), "252826822479841", "perfect power limits for 15964377" );

is( approx_in_range(1571,2048383), 2048383, "perfect power approx for 1571" );
is( approx_in_range(59643,3373286400), 3373286400, "perfect power approx for 59643" );
is( approx_in_range(15964377,"252826822479841"), "252826822479841", "perfect power approx for 15964377" );



sub cmp_rn {
  my($n,$rn) = @_;
  return 'nth lower' unless nth_perfect_power_lower($n) <= $rn;
  return 'nth upper' unless nth_perfect_power_upper($n) >= $rn;
  return 'pc lower' unless perfect_power_count_lower($rn) <= $n;
  return 'pc upper' unless perfect_power_count_upper($rn) >= $n;
  $rn;
}

sub approx_in_range {
  my($n,$rn) = @_;
  my $arn = nth_perfect_power_approx($n);
  my $an  = perfect_power_count_approx($rn);
  return 'nth approx too low' if $arn < ($rn-$rn/100);
  return 'nth approx too high' if $arn > ($rn+$rn/100);
  return 'count approx too low' if $an < ($n-$n/100);
  return 'count approx too high' if $an > ($n+$n/100);
  $rn;
}
