#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util
   qw/todigits fromdigits todigitstring sumdigits vecsum factorial/;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
#my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
#my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
#my $usegmp= Math::Prime::Util::prime_get_config->{'gmp'};
#$use64 = 0 if $use64 && 18446744073709550592 == ~0;

plan tests => 0 + 6 + 3 + 2 + 1  # fromdigits
                + 6 + 4 + 1      # todigits
                + 4 + 2*$extra   # sumdigits
                + 12;

###### fromdigits
is(fromdigits([0,1,1,0,1],2), 13, "fromdigits binary with leading 0");
is(fromdigits([1,1,0,1],2), 13, "fromdigits binary");
is(fromdigits([0,1,1,0,1]), 1101, "fromdigits decimal");
is(fromdigits([0,1,1,0,1],3), 37, "fromdigits base 3");
is(fromdigits([0,1,1,0,1],16), 4353, "fromdigits base 16");
is(fromdigits([0,1,1,0,2216],16), 6568, "fromdigits base 16 with overflow");

# Pari/GP seems to have a carry issue with base 2
is(fromdigits([7,999,44],5), 7*5**2 + 999*5 + 44*1, "fromdigits base 5 with carry");
is(fromdigits([7,999,44],3), 7*3**2 + 999*3 + 44*1, "fromdigits base 3 with carry");
is(fromdigits([7,999,44],2), 7*2**2 + 999*2 + 44*1, "fromdigits base 2 with carry");

is(fromdigits("1f",16), 31, "fromdigits hex string");
is(fromdigits("24"), 24, "fromdigits decimal");

is(fromdigits("zzzyzzzyzzzyzzzy",36), "7958656371562241451187966", "fromdigits with Large base 36 number");

###### todigits
is_deeply([todigits(0)], [], "todigits 0");
is_deeply([todigits(1)], [1], "todigits 1");
is_deeply([todigits(77)], [7,7], "todigits 77");
is_deeply([todigits(77,2)], [1,0,0,1,1,0,1], "todigits 77 base 2");
is_deeply([todigits(77,3)], [2,2,1,2], "todigits 77 base 3");
is_deeply([todigits(77,21)], [3,14], "todigits 77 base 21");

is_deeply([todigits(900,2)], [1,1,1,0,0,0,0,1,0,0], "todigits 900 base 2");
is_deeply([todigits(900,2,0)], [], "todigits 900 base 2 len 0");
is_deeply([todigits(900,2,3)], [1,0,0], "todigits 900 base 2 len 3");
is_deeply([todigits(900,2,32)], [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,1,0,0], "todigits 900 base 2 len 32");

is(vecsum(todigits("293852387239761276234029385230912847923872323")), 201, "vecsum of todigits of bigint");

###### sumdigits
is(sumdigits("-45.36"), 4+5+3+6, "sumdigits(-45.36)");
{
  my @sumd   = map { sumdigits($_) } 0 .. 1000;
  my @splitd = map { vecsum(split(//,$_)) } 0 .. 1000;
  is_deeply( \@sumd, \@splitd, "sumdigits 0 to 1000");
}
is(sumdigits("0x3290f8E"), 51, "sumdigits hex");
is(sumdigits("293852387239761276234029385230912847923872323"), 201, "sumdigits bigint");
if ($extra) {
  is(sumdigits(factorial(1000)), 10539, "sumdigits 1000!");
  is(sumdigits(factorial(10000)), 149346, "sumdigits 10000!");
}

###### examples from Wolfram docs
is_deeply([todigits(1234135634,16)], [4,9,8,15,6,10,5,2], "todigits 1234135634 base 16");
is_deeply([todigits(56,2,8)], [0,0,1,1,1,0,0,0], "todigits 56 base 2 len 8");
is(fromdigits([todigits(56,2,8)],2), 56, "fromdigits of previous");
is(todigitstring(56,2), "111000", "56 as binary string");
is(fromdigits(todigitstring(56,2),2), 56, "fromdigits of previous");
is(todigitstring(37,2), "100101", "todigitstring 37");
is(fromdigits([5,1,2,8]), 5128, "fromdigits 5128 base 10");
is(fromdigits([1,0,1,1,0,1,1],2), 91, "fromdigits 91 base 2");
is(fromdigits("1923"), 1923, "fromdigits 1923 base 10");
is(fromdigits("1011011",2), 91, "fromdigits 91 base 2");
is(fromdigits([7,11,0,0,0,122]), 810122, "fromdigits with carry");
is_deeply([todigits(6345354, 10, 4)], [5,3,5,4], "only last 4 digits");
