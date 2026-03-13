#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/todigits fromdigits todigitstring sumdigits
                         is_palindrome digital_root mult_digital_root
                         vecsum factorial/;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};

plan tests => 1    # fromdigits
            + 1    # todigits
            + 1    # sumdigits
            + 1    # todigitstring
            + 1    # is_palindrome
            + 1    # digital_root
            + 1    # mult_digital_root
            + 1;   # examples from Wolfram docs


subtest 'fromdigits', sub {
  is(fromdigits([0,1,1,0,1],2), 13, "fromdigits binary with leading 0");
  is(fromdigits([1,1,0,1],2), 13, "fromdigits binary");
  is(fromdigits([0,1,1,0,1]), 1101, "fromdigits decimal");
  is(fromdigits([0,1,1,0,1],3), 37, "fromdigits base 3");
  is(fromdigits([0,1,1,0,1],16), 4353, "fromdigits base 16");
  is(fromdigits([0,1,1,0,2216],16), 6568, "fromdigits base 16 with overflow");
  is(fromdigits([7,999,44],5), 7*5**2 + 999*5 + 44*1, "fromdigits base 5 with carry");
  is(fromdigits([7,999,44],3), 7*3**2 + 999*3 + 44*1, "fromdigits base 3 with carry");
  is(fromdigits([7,999,44],2), 7*2**2 + 999*2 + 44*1, "fromdigits base 2 with carry");
  is("".fromdigits([1..15,1..15,1..15],16), "108977460683796539709587792812439445667270661579197935", "fromdigits base 16 with many digits");

  is(fromdigits("1f",16), 31, "fromdigits hex string");
  is(fromdigits("24"), 24, "fromdigits decimal");

  is("".fromdigits("zzzyzzzyzzzyzzzy",36), "7958656371562241451187966", "fromdigits with Large base 36 number");
};

subtest 'todigits', sub {
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

  is_deeply([todigits(-143)], [1,4,3], "todigits ignores negative sign");
};

subtest 'sumdigits', sub {
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

  is(sumdigits(-143), 8, "sumdigits ignores negative sign");
};

subtest 'todigitstring', sub {
  is(todigitstring("3" x 21, 3), "10001020211011120202011020201202220201012100", "todigitstring base 3");
  is(todigitstring("7" x 26, 9), "1303055203367717374834745502", "todigitstring base 9");
  is(todigitstring("9" x 27, 11), "92586630a001888a8112250349", "todigitstring base 11");
  is(todigitstring(-143,16), "8f", "todigitstring ignores negative sign");
  is(todigitstring(12345,8,10), "0000030071", "todigitstring will 0 pad");
};

subtest 'is_palindrome', sub {
  is_deeply([map {is_palindrome($_)} 0..9,11,22,33,121,959,12321],
            [map { 1 }               0..9,11,22,33,121,959,12321],
            "palindromes");
  is_deeply([map {is_palindrome($_)} 10,12,19,34,98,12322],
            [map { 0 }               10,12,19,34,98,12322],
            "not palindromes");

  is(is_palindrome(9,2),1,"9 is a base-2 palindrome");
  is(is_palindrome(10,2),0,"10 is not a base-2 palindrome");
  is(is_palindrome(130,3),1,"130 is a base-3 palindrome");
  is(is_palindrome(130,4),1,"130 is a base-4 palindrome");

  my $m = 719848917;
  ok(is_palindrome($m) && is_palindrome($m,2) && is_palindrome($m,8),"$m is a palindrome in bases 2, 8, and 10");
  my $k = 532900;
  ok(is_palindrome($k,3) && is_palindrome($k,7),"$k is a palindrome in bases 3 and 7");
};

subtest 'digital_root' => sub {
  is_deeply( [map { digital_root($_) } 0..19],
             [0,1,2,3,4,5,6,7,8,9,1,2,3,4,5,6,7,8,9,1],
             "digital_root 0..19" );
  is( digital_root(123456789),    9, "digital_root(123456789)" );
  is( digital_root(1000000000),   1, "digital_root(10^9)" );
  is( digital_root(999999998),    8, "digital_root(999999998)" );
  is_deeply( [map { digital_root($_, 16) } 0..16],
             [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,1],
             "digital_root 0..16 base 16" );
  is( digital_root(49,  6),  4, "digital_root(49, 6)" );
  is( digital_root(255, 16), 15, "digital_root(255, 16)" );
  is( digital_root("99999999999999999999"),  9, "digital_root bigint (20 nines)" );
  is( digital_root("10000000000000000000"),  1, "digital_root bigint (10^19)" );
};

subtest 'mult_digital_root' => sub {
  is_deeply( [map { mult_digital_root($_) } 0..9],
             [0,1,2,3,4,5,6,7,8,9],
             "mult_digital_root 0..9" );
  is( mult_digital_root(23),  6, "23 -> 6" );         # 2*3=6
  is( mult_digital_root(39),  4, "39 -> 4" );         # 27->14->4
  is( mult_digital_root(77),  8, "77 -> 8" );         # 49->36->18->8
  is( mult_digital_root(679), 6, "679 -> 6" );        # 378->168->48->32->6
  is( mult_digital_root(999), 2, "999 -> 2" );        # 729->126->12->2
  is( mult_digital_root(13311), 9, "13311 -> 9" );    # 13311->9
  is( mult_digital_root(531537), 5, "531537 -> 5" );  # 531537,1575,175,35,15,5
  is( mult_digital_root(6,  6),   0, "mdr(6,  6)=0" );   # "10"_6 -> 0
  is( mult_digital_root(7,  6),   1, "mdr(7,  6)=1" );   # "11"_6 -> 1
  is( mult_digital_root(14, 6),   4, "mdr(14, 6)=4" );   # "22"_6 -> 4
  is( mult_digital_root(255, 16), 14, "mdr(255,16)=14" ); # FF->225->E1->14
  is( mult_digital_root("1" x 50), 1, "mdr('1'x50)=1" );
};

subtest 'Wolfram examples', sub {
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
};
