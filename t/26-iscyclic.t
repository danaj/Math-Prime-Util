#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/ is_cyclic vecsum /;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};

plan tests => 4;

is_deeply( [grep { is_cyclic($_) } -20 .. 20],
           [1,2,3,5,7,11,13,15,17,19],
           "Cyclic numbers -20 to 20: primes plus 15" );

is_deeply( [grep { is_cyclic($_) } (qw/9 15 21 25 27 33 35 39 45 49 51 55 57 63 65 69 75 77 81 85 87 91 93 95 99 105 111 115 117 119 121 123 125 129 133 135 141 143 145 147 153 155 159 161 165 169 171 175 177 183 185 187 189 195/)],
           [qw/15 33 35 51 65 69 77 85 87 91 95 115 119 123 133 141 143 145 159 161 177 185 187/],
           "Cyclic composites under 200" );

is( is_cyclic(32753), 1, "32753 (the 10,000th cyclic number) is cyclic" );

SKIP: {
  skip "count 10k cyclic numbers only with extended testing", 1 unless $extra;
  is( vecsum(map { is_cyclic($_) } 1..32753), 10000, "count of cyclic numbers to 32753 = 10000" );
}
