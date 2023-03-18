#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_delicate_prime/;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};

my @a050249 = (294001, 505447, 584141, 604171, 971767, 1062599, 1282529, 1524181, 2017963, 2474431, 2690201, 3085553, 3326489, 4393139, 5152507, 5564453, 5575259, 6173731, 6191371, 6236179, 6463267, 6712591, 7204777, 7469789, 7469797);

my @a186995 = (127, 2, 373, 83, 28151, 223, 6211, 2789, 294001, 3347, 20837899, 4751, 6588721, 484439, 862789, 10513, 2078920243, 10909, 169402249, 2823167, 267895961, 68543, 1016960933671, 181141, 121660507, 6139219, 11646280537, 488651);

plan tests => 4 + 5 + 3*$extra;

# If we forget to change digits to zeros:
is(is_delicate_prime(504991), 0, "is_delicate_prime(504991) = 0");
# If we don't change the leading digit to zero (OEIS A158124):
is(is_delicate_prime(929573), 0, "is_delicate_prime(929573) = 0");

is_deeply( [map { is_delicate_prime($_) } @a050249], [map { 1 } @a050249], "is_delicate_prime(n) = 1 for first ".scalar(@a050249)." known.");

is_deeply( [map { is_delicate_prime($a186995[$_-2],$_) } 2..2+$#a186995], [map { 1 } @a186995], "first delicate primes for bases 2 to ".(2+$#a186995).".");

is_deeply( [grep { is_delicate_prime($_,2) } 2..277], [qw/127 173 191 223 233 239 251 257 277/], "First 9 delicate primes base 2" );
is_deeply( [grep { is_delicate_prime($_,3) } 2..283], [qw/2 7 13 19 31 41 149 239 283/], "First 9 delicate primes base 3" );
is_deeply( [grep { is_delicate_prime($_,5) } 2..739], [qw/83 233 277 397 487 509 593 647 739/], "First 9 delicate primes base 5" );
is_deeply( [grep { is_delicate_prime($_,7) } 2..1361], [qw/223 409 491 587 701 1051 1163 1237 1361/], "First 9 delicate primes base 7" );
is_deeply( [grep { is_delicate_prime($_,11) } 2..11411], [qw/3347 3761 5939 6481 8831 9257 9749 10487 11411/], "First 9 delicate primes base 11" );

if ($extra) {
  is(is_delicate_prime("999999999998832431"), 1, "is_delicate_prime(999999999998832431) = 1");
  is(is_delicate_prime("999999999999999543767"), 1, "is_delicate_prime(999999999999999543767) = 1");
  is_deeply( [grep { is_delicate_prime($_) } 0..1062599], [grep { $_ <= 1062599 } @a050249], "is_practical(0 .. 1062599) returns first 6 delicate primes");
}

#  500 digits: 2021 * powint(10, 500-4) + 7543997
# 1000 digits: 2021 * powint(10,1000-4) + 2550219
# 2021 digits: 2021 * powint(10,2021-4) + 4523733
#  500 digits: 2021 * 10^(500-4) + 7543997
# 1000 digits: 2021 * 10^(1000-4) + 2550219
# 2021 digits: 2021 * 10^(2021-4) + 4523733
