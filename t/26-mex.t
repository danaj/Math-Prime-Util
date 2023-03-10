#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/vecmex vecpmex powint forpart/;

plan tests => 0
            + 8    # vecmex
            + 8    # vecpmex
            + 1    # vecpmex
            + 0;

###### vecmex
is(vecmex(), 0, "vecmex() = 0");
is(vecmex(0), 1, "vecmex(0) = 1");
is(vecmex(1), 0, "vecmex(1) = 0");
is(vecmex(1,2,4), 0, "vecmex(1,2,4) = 0");
is(vecmex(0,1,2,4), 3, "vecmex(0,1,2,4) = 3");
is(vecmex(0,1,24,4), 2, "vecmex(0,1,24,4) = 2");
is(vecmex(4,2,1,0), 3, "vecmex(4,2,1,0) = 3");
is(vecmex(3,powint(10,20),0,2), 1, "vecmex(3,10^20,0,2) = 1");

###### vecpmex
is(vecpmex(), 1, "vecpmex() = 1");
is(vecpmex(1), 2, "vecpmex(1) = 2");
is(vecpmex(2), 1, "vecpmex(2) = 1");
is(vecpmex(2,3,5), 1, "vecpmex(2,3,5) = 1");
is(vecpmex(1,2,3,5), 4, "vecpmex(1,2,3,5) = 4");
is(vecpmex(1,2,24,5), 3, "vecpmex(1,2,24,5) = 3");
is(vecpmex(5,3,2,1), 4, "vecpmex(5,3,2,1) = 4");
is(vecpmex(4,powint(10,20),1,3), 2, "vecpmex(4,10^20,1,3) = 2");

###### A022567

# See https://arxiv.org/pdf/2303.05332.pdf

sub sigmamex { my $s=0; forpart { $s += vecpmex(@_); } $_[0]; $s; }
is_deeply([map { sigmamex($_) } 1..10], [2,3,6,9,14,22,32,46,66,93], "sigmaxmex(1..10)");
