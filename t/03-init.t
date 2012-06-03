#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util qw/prime_precalc prime_free/;

use Test::More  tests => 3;


# How do we test these?  Use some private functions and assume things about
# internal behavior.  Not the best idea.

can_ok( 'Math::Prime::Util', '_get_prime_cache_size' );

my $init_size = Math::Prime::Util::_get_prime_cache_size;

prime_precalc(10_000_000);

cmp_ok( Math::Prime::Util::_get_prime_cache_size, '>', $init_size, "Internal space grew after large precalc" );

prime_free;

is( Math::Prime::Util::_get_prime_cache_size, $init_size, "Internal space went back to original size after free" );
