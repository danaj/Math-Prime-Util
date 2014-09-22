#!/usr/bin/env perl
use strict;
use warnings;
use ntheory qw/is_prime/;

use Test::More  tests => 1;

ok(is_prime(7), "ntheory can do is_prime");
