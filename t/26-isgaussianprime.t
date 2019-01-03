#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_gaussian_prime/;

plan tests => 7;

ok( !is_gaussian_prime(29,0), "29 is not a Gaussian Prime" );
ok(  is_gaussian_prime(31,0), "31 is a Gaussian Prime" );
ok( !is_gaussian_prime(0,-29), "0-29i is not a Gaussian Prime" );
ok(  is_gaussian_prime(0,-31), "0-31i is a Gaussian Prime" );
ok(  is_gaussian_prime(58924,132000511), "58924+132000511i is a Gaussian Prime" );
ok(  is_gaussian_prime(519880,-2265929), "519880-2265929i is a Gaussian Prime" );
ok( !is_gaussian_prime(20571,150592260), "20571+150592260i is not a Gaussian Prime" );
