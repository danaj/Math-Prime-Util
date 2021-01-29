#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_delicate_prime/;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};

my @a050249 = (294001, 505447, 584141, 604171, 971767, 1062599, 1282529, 1524181, 2017963, 2474431, 2690201, 3085553, 3326489, 4393139, 5152507, 5564453, 5575259, 6173731, 6191371, 6236179, 6463267, 6712591, 7204777, 7469789, 7469797);

plan tests => 3 + 3*$extra;

# If we forget to change digits to zeros:
is(is_delicate_prime(504991), 0, "is_delicate_prime(504991) = 0");
# If we don't change the leading digit to zero (OEIS A158124):
is(is_delicate_prime(929573), 0, "is_delicate_prime(929573) = 0");

is_deeply( [map { is_delicate_prime($_) } @a050249], [map { 1 } @a050249], "is_delicate_prime(n) = 1 for first ".scalar(@a050249)." known.");

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
