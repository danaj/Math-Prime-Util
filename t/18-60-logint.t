#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/logint powint/;

plan tests => 2 + 4 + 2 + 3 + 1 + 4;

ok( !defined eval { logint(0,10) } && !defined eval { logint(-5,10) },
    "logint(n,base): n must be at least 1" );
ok( !defined eval { logint(10,1) } && !defined eval { logint(10,0) },
    "logint(n,base): base must be at least 2" );

is_deeply( [map { logint($_,2) } 1..200], [map { int(log($_)/log(2)+1e-10) } 1..200], "logint base 2: 0 .. 200" );
is_deeply( [map { logint($_,3) } 1..200], [map { int(log($_)/log(3)+1e-10) } 1..200], "logint base 3: 0 .. 200" );
is_deeply( [map { logint($_,5) } 1..200], [map { int(log($_)/log(5)+1e-10) } 1..200], "logint base 5: 0 .. 200" );
is_deeply( [map { logint($_,10) } 1..200], [map { length($_)-1 } 1..200], "logint base 10: 0 .. 200" );
{
  eval { logint(16, 2, undef); 1 };
  like($@, qr/scalar reference/i, "logint rejects undef power reference");

  eval { logint(16, 2, []); 1 };
  like($@, qr/scalar reference/i, "logint rejects array power reference");
}
{
  my $be;
  is( logint(19284098234,16,\$be), 8, "logint(19284098234,16) = 8" );
  is( "".$be, 16**8, "power is 16^8" );
}
{
  my $huge_base = "18446744073709551617";
  my $be;
  is( logint(16, $huge_base), 0, "logint with huge base returns 0 without power reference" );
  is( logint(16, $huge_base, \$be), 0, "logint with huge base returns 0 with power reference" );
  is( $be, 1, "logint huge base power reference is 1" );
}
is(logint(powint(3,1000), powint(10,24)), 19, "logint with bigint n and bigint base");
is(logint(58,~0), 0, "logint(58,~0) = 0");
is(logint("1329227995784915872903807060280344576",2), 120, "logint(2^120,2) = 120");
