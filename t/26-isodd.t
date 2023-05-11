#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_odd is_even/;

plan tests => 4;

is_deeply( [map { is_odd($_) } -50..50], [map { $_%2 } -50..50], "is_odd(-50..50)" );
is_deeply( [map { is_even($_) } -50..50], [map { 1-($_%2) } -50..50], "is_even(-50..50)" );
is( is_odd("9832494848388408230982349082309487"), 1, "is_odd bigint" );
is( is_even("98324948483884082309823490823094875"), 0, "is_even bigint" );
