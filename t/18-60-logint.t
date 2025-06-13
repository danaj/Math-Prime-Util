#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/logint/;

plan tests => 2 + 4 + 4;

ok( !defined eval { logint(0,10) } && !defined eval { logint(-5,10) },
    "logint(n,base): n must be at least 1" );
ok( !defined eval { logint(10,1) } && !defined eval { logint(10,0) },
    "logint(n,base): base must be at least 2" );

is_deeply( [map { logint($_,2) } 1..200], [map { int(log($_)/log(2)+1e-10) } 1..200], "logint base 2: 0 .. 200" );
is_deeply( [map { logint($_,3) } 1..200], [map { int(log($_)/log(3)+1e-10) } 1..200], "logint base 3: 0 .. 200" );
is_deeply( [map { logint($_,5) } 1..200], [map { int(log($_)/log(5)+1e-10) } 1..200], "logint base 5: 0 .. 200" );
is_deeply( [map { logint($_,10) } 1..200], [map { length($_)-1 } 1..200], "logint base 10: 0 .. 200" );
{
  my $be;
  is( logint(19284098234,16,\$be), 8, "logint(19284098234,16) = 8" );
  is( $be, 16**8, "power is 16^8" );
}
is(logint(58,~0), 0, "logint(58,~0) = 0");
is(logint("1329227995784915872903807060280344576",2), 120, "logint(2^120,2) = 120");
