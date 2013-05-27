#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/primes forprimes prime_iterator/;

my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;

plan tests => 8 + 5
            + 3 + 7
            + 2;

ok(!eval { forprimes { 1 } undef; },   "forprimes undef");
ok(!eval { forprimes { 1 } 2, undef; },   "forprimes 2,undef");
ok(!eval { forprimes { 1 } undef, 2; },   "forprimes 2,undef");
# This is caught at compile type because of the prototype
#ok(!eval { forprimes { 1 } 2, 3, 4; },   "forprimes 2,3,4");
ok(!eval { forprimes { 1 } -2, 3; },   "forprimes -2,3");
ok(!eval { forprimes { 1 } 2, -3; },   "forprimes 2,-3");
ok(!eval { forprimes { 1 } "abc"; },   "forprimes abc");
ok(!eval { forprimes { 1 } 2, "abc"; },   "forprimes 2, abc");
ok(!eval { forprimes { 1 } 5.6; },   "forprimes abc");

{ my @t; forprimes { push @t, $_ } 50;
  is_deeply( [@t], [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47], "forprimes 50" );
}
{ my @t; forprimes { push @t, $_ } 2,20;
  is_deeply( [@t], [2,3,5,7,11,13,17,19], "forprimes 2,20" );
}
{ my @t; forprimes { push @t, $_ } 20,30;
  is_deeply( [@t], [23,29], "forprimes 20,30" );
}
{ my @t; forprimes { push @t, $_ } 199, 223;
  is_deeply( [@t], [199,211,223], "forprimes 199,223" );
}
{ my @t; forprimes { push @t, $_ } 31398, 31468;
  is_deeply( [@t], [], "forprimes 31398,31468 (empty region)" );
}

ok(!eval { prime_iterator(-2); }, "iterator -2");
ok(!eval { prime_iterator("abc"); }, "iterator abc");
ok(!eval { prime_iterator(4.5); }, "iterator 4.5");

{ my $it = prime_iterator();
  is_deeply( [map { $it->() } 1..10], [2,3,5,7,11,13,17,19,23,29], "iterator first 10 primes" );
}
{my $it = prime_iterator(47);
  is_deeply( [map { $it->() } 1..5], [47,53,59,61,67], "iterator 5 primes starting at 47" );
}
{my $it = prime_iterator(199);
  is_deeply( [map { $it->() } 1..3], [199,211,223], "iterator 3 primes starting at 199" );
}
{my $it = prime_iterator(200);
  is_deeply( [map { $it->() } 1..3], [211,223,227], "iterator 3 primes starting at 200" );
}
{my $it = prime_iterator(31397);
  is_deeply( [map { $it->() } 1..3], [31397,31469,31477], "iterator 3 primes starting at 31397" );
}
{my $it = prime_iterator(31396);
  is_deeply( [map { $it->() } 1..3], [31397,31469,31477], "iterator 3 primes starting at 31396" );
}
{my $it = prime_iterator(31398);
  is_deeply( [map { $it->() } 1..3], [31469,31477,31481], "iterator 3 primes starting at 31398" );
}

# For fun, nest them.
{
  my @t;
  forprimes {
    forprimes {
      forprimes { push @t, $_ } $_,$_+10;
    } 10*$_,10*$_+10;
  } 10;
  is_deeply( [@t], [qw/23 29 31 29 31 37 31 37 41 37 41 43 47 53 59 61 59 61 67 71 73 79 73 79 83 79 83 89/], "triple nested forprimes" );
}
{
  my @t;
  my $ita = prime_iterator();
  while ((my $a = $ita->()) <= 10) {
    my $itb = prime_iterator(10*$a);
    while ((my $b = $itb->()) <= 10*$a+10) {
      my $itc = prime_iterator($b);
      while ((my $c = $itc->()) <= $b+10) {
        push @t, $c;
      }
    }
  }
  is_deeply( [@t], [qw/23 29 31 29 31 37 31 37 41 37 41 43 47 53 59 61 59 61 67 71 73 79 73 79 83 79 83 89/], "triple nested iterator" );
}
