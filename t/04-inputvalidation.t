#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/next_prime/;

plan tests => 22;

eval { next_prime(undef); };
like($@, qr/^Parameter must be defined/, "next_prime(undef)");
eval { next_prime(""); };
like($@, qr/^Parameter must be a positive integer/, "next_prime('')");
eval { next_prime(-4); };
like($@, qr/^Parameter '-4' must be a positive integer/, "next_prime(-4)");
eval { next_prime("-"); };
like($@, qr/^Parameter '-' must be a positive integer/, "next_prime('-')");
eval { next_prime("+"); };
like($@, qr/^Parameter '\+' must be a positive integer/, "next_prime('+')");

# +4 is fine
is(next_prime("+4"), 5, "next_prime('+4') works");
# ++4 does not
eval { next_prime("++4"); };
like($@, qr/^Parameter '\+\+4' must be a positive integer/, "next_prime('++4')");
eval { next_prime("+-4"); };
like($@, qr/^Parameter '\+\-4' must be a positive integer/, "next_prime('+-4')");

# Test leading zeros
is(next_prime("0004"), 5, "next_prime('0004') works");
is(next_prime("+0004"), 5, "next_prime('+0004') works");
eval { next_prime("-0004"); };
like($@, qr/^Parameter '\-0004' must be a positive integer/, "next_prime('-0004')");

eval { next_prime("a"); };
like($@, qr/^Parameter 'a' must be a positive integer/, "next_prime('a')");
eval { next_prime(5.6); };
like($@, qr/^Parameter '5.6' must be a positive integer/, "next_prime('5.6')");

# 5.0 should be ok.
is(next_prime(5.0), 7, "next_prime(5.0) works");
eval { next_prime("4e"); };
like($@, qr/^Parameter '4e' must be a positive integer/, "next_prime('4e')");
eval { next_prime("1.1e12"); };
like($@, qr/^Parameter '1.1e12' must be a positive integer/, "next_prime('1.1e12')");

# 1e8 as a string will fail, as a number will work.
eval { next_prime("1e8"); };
like($@, qr/^Parameter '1e8' must be a positive integer/, "next_prime('1e8')");
is(next_prime(1e8), 100000007, "next_prime(1e8) works");

eval { next_prime("NaN"); };
like($@, qr/^Parameter 'NaN' must be a positive integer/, "next_prime('NaN')");

# The actual strings can be implementation specific
eval { next_prime(0+'inf'); };
like($@, qr/must be a positive integer/, "next_prime(0+'inf')");
eval { next_prime(20**20**20); };
like($@, qr/must be a positive integer/, "next_prime(20**20**20)");

eval { next_prime("11111111111111111111111111111111111111111x"); };
like($@, qr/must be a positive integer/, "next_prime('111...111x')");
