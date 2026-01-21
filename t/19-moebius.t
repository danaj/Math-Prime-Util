#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/moebius mertens/;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
my $usegmp= Math::Prime::Util::prime_get_config->{'gmp'};
my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
$use64 = 0 if $use64 && 18446744073709550592 == ~0;

my %mertens = (
        1 =>    1,
        2 =>    0,
        3 =>   -1,
        4 =>   -1,
        5 =>   -2,
       10 =>   -1,
      100 =>    1,
     1000 =>    2,
    10000 =>  -23,
        8 =>   -2,
       16 =>   -1,
       32 =>   -4,
       64 =>   -1,
      128 =>   -2,
      256 =>   -1,
      512 =>   -4,
     1024 =>   -4,
     2048 =>    7,
     4096 =>  -19,
     8192 =>   22,
);
my %big_mertens = (
   100000 =>  -48,
   444444 =>  -37,
  1000000 =>  212,
 10000000 => 1037,
);
delete $big_mertens{10000000} unless $extra || $usexs;

if ($extra && $use64) {
  %big_mertens = ( %big_mertens,
          2 =>  0,      # A087987, mertens at primorials
          6 => -1,
         30 => -3,
        210 => -1,
       2310 => -1,
      30030 => 16,
     510510 => -25,
    9699690 => 278,
  223092870 => 3516,

    6433477 => 900,     # 30^2
  109851909 => -4096,   # A084235, 2^12

      2**14 =>  -32,    # A084236
      2**15 =>   26,
      2**16 =>   14,
      2**17 =>  -20,
      2**18 =>   24,
      2**19 => -125,
      2**20 =>  257,
      2**21 => -362,
      2**22 =>  228,
      2**23 =>  -10,

     10**8  => 1928,
     10**9  => -222,
  1*10**10  => -33722,  # From Deleglise and Rivat
  2*10**10  =>  48723,
  3*10**10  =>  42411,
  4*10**10  => -25295,
    10**11  => -87856,
  );
}
# These are slow with XS, and *really* slow with PP.
if (!$usexs) {
  %big_mertens = map { $_ => $big_mertens{$_} }
                 grep { $_ < 100000000 }
                 keys %big_mertens;
}

plan tests => 1 + 5 + 2 + 2 + 3 + scalar(keys %big_mertens);

ok(!eval { moebius(0); }, "moebius(0)");

is_deeply( [map { moebius($_) } 1 .. 20],
           [1,-1,-1,0,-1,1,-1,0,0,1,-1,0,-1,1,1,0,-1,0,-1,0],
           "moebius 1 .. 20 (single)" );
is_deeply( [moebius(1,20)],
           [1,-1,-1,0,-1,1,-1,0,0,1,-1,0,-1,1,1,0,-1,0,-1,0],
           "moebius 1 .. 20 (range)" );
# moebius uses |n| so negative inputs reflect about zero.
is_deeply( [map { moebius(-$_) } 1 .. 20],
           [1,-1,-1,0,-1,1,-1,0,0,1,-1,0,-1,1,1,0,-1,0,-1,0],
           "moebius -1 .. -20 (single)" );
is_deeply( [moebius(-14,-9)],
           [1,-1,0,-1,1,0],
           "moebius -14 .. -9 (range)" );
is_deeply( [moebius(-7,5)],
           [-1,1,-1,0,-1,-1,1,0,1,-1,-1,0,-1],
           "moebius -7 .. 5 (range)" );

is( moebius(3*5*7*11*13), -1, "moebius(3*5*7*11*13) = -1" );
is( moebius("20364840299624512075310661735"), 1, "moebius(73#/2) = 1" );

# near end points
is_deeply([[moebius(4294967293,4294967295)],
           [moebius(4294967293,4294967296)],
           [moebius(4294967295,4294967297)],
           [moebius(4294967296,4294967298)]],
          [[1,1,-1],[1,1,-1,0],[-1,0,1],[0,1,-1]],
          "moebius ranges around 2^32");

SKIP: {
  skip "ranges around 2^64 only on 64-bit",1 unless $use64;
  is_deeply([[moebius("18446744073709551613","18446744073709551615")],
             [moebius("18446744073709551613","18446744073709551616")],
             [moebius("18446744073709551615","18446744073709551617")]],
            [[-1,0,-1],[-1,0,-1,0],[-1,0,1]],
            "moebius ranges around 2^64");
}

{
  my(@mert_sum1, @mert_sum2, @mertens, @expect, $M);
  while (my($n, $val) = each (%mertens)) {
    $M = 0;  $M += moebius($_) for 1 .. $n;
    push @mert_sum1, $M;
    $M = 0;  $M += $_ for moebius(1,$n);
    push @mert_sum2, $M;
    push @mertens, mertens($n);
    push @expect, $val;
  }
  is_deeply( \@mert_sum1, \@expect, "sum(moebius(k) for k=1..n)   small n" );
  is_deeply( \@mert_sum2, \@expect, "sum(moebius(1,n))   small n" );
  is_deeply( \@mertens, \@expect, "mertens(n)   small n" );
}
while (my($n, $mertens) = each (%big_mertens)) {
  is( mertens($n), $mertens, "mertens($n)" );
}
