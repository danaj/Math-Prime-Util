#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/primes prev_prime next_prime
                         forprimes prime_iterator prime_iterator_object/;
use Math::BigInt try => "GMP,Pari";
use Math::BigFloat;

my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;

plan tests => 8        # forprimes errors
            + 12 + 5   # forprimes simple
            + 3        # iterator errors
            + 7        # iterator simple
            + 2        # forprimes/iterator nesting
            + 2        # forprimes BigInt/BigFloat
            + 3        # oo iterator errors
            + 7        # oo iterator simple
            + 25       # oo iterator methods
            + 0;

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

{my @t; forprimes {push @t,$_} 1; is_deeply( [@t], [], "forprimes 1" ); }
{my @t; forprimes {push @t,$_} 2; is_deeply( [@t], [2], "forprimes 3" ); }
{my @t; forprimes {push @t,$_} 3; is_deeply( [@t], [2,3], "forprimes 3" ); }
{my @t; forprimes {push @t,$_} 4; is_deeply( [@t], [2,3], "forprimes 4" ); }
{my @t; forprimes {push @t,$_} 5; is_deeply( [@t], [2,3,5], "forprimes 5" ); }
{my @t; forprimes {push @t,$_} 3,5; is_deeply( [@t], [3,5], "forprimes 3,5" ); }
{my @t; forprimes {push @t,$_} 3,6; is_deeply( [@t], [3,5], "forprimes 3,6" ); }
{my @t; forprimes {push @t,$_} 3,7; is_deeply( [@t], [3,5,7], "forprimes 3,7" ); }
{my @t; forprimes {push @t,$_} 5,7; is_deeply( [@t], [5,7], "forprimes 5,7" ); }
{my @t; forprimes {push @t,$_} 6,7; is_deeply( [@t], [7], "forprimes 6,7" ); }
{my @t; forprimes {push @t,$_} 5,11; is_deeply( [@t], [5,7,11], "forprimes 5,11" ); }
{my @t; forprimes {push @t,$_} 7,11; is_deeply( [@t], [7,11], "forprimes 7,11" ); }

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

# With BigInt and BigFloat objects
{ my @t;
  forprimes { push @t, $_ } Math::BigInt->new("5"), Math::BigInt->new("11");
  is_deeply( [@t], [5,7,11], "forprimes with BigInt range" );
}
{ my @t;
  forprimes { push @t, $_ } Math::BigFloat->new("5"), Math::BigFloat->new("11");
  is_deeply( [@t], [5,7,11], "forprimes with BigFloat range" );
}

# Test new object iterator
ok(!eval { prime_iterator_object(-2); }, "iterator -2");
ok(!eval { prime_iterator_object("abc"); }, "iterator abc");
ok(!eval { prime_iterator_object(4.5); }, "iterator 4.5");

{ my $it = prime_iterator_object();
  is_deeply( [map { $it->iterate() } 1..10], [2,3,5,7,11,13,17,19,23,29], "iterator first 10 primes" );
}
{my $it = prime_iterator_object(47);
  is_deeply( [map { $it->iterate() } 1..5], [47,53,59,61,67], "iterator 5 primes starting at 47" );
}
{my $it = prime_iterator_object(199);
  is_deeply( [map { $it->iterate() } 1..3], [199,211,223], "iterator 3 primes starting at 199" );
}
{my $it = prime_iterator_object(200);
  is_deeply( [map { $it->iterate() } 1..3], [211,223,227], "iterator 3 primes starting at 200" );
}
{my $it = prime_iterator_object(31397);
  is_deeply( [map { $it->iterate() } 1..3], [31397,31469,31477], "iterator 3 primes starting at 31397" );
}
{my $it = prime_iterator_object(31396);
  is_deeply( [map { $it->iterate() } 1..3], [31397,31469,31477], "iterator 3 primes starting at 31396" );
}
{my $it = prime_iterator_object(31398);
  is_deeply( [map { $it->iterate() } 1..3], [31469,31477,31481], "iterator 3 primes starting at 31398" );
}

{
  my $it = prime_iterator_object;
  do { $it->next } for 1..10;
  is( $it->value(), 31, "iterator object moved forward 10 now returns 31");
  $it->prev;
  is( $it->value(), 29, "iterator object moved back now returns 29");
  is( $it->iterate(), 29, "iterator object iterates to 29");
  is( $it->iterate(), 31, "iterator object iterates to 31");
  $it->rewind->next->next->next->prev;
  is( $it->value(), 5, "iterator object rewind and move returns 5");
  # Validate that it automatically handles bigint range traversal.
  my $top_prime = prev_prime(~0);
  my $big_prime = next_prime(Math::BigInt->new(''.~0));
  ok( $big_prime > ~0, "internal check, next_prime on big int works");
  $it->rewind($top_prime);
  is( $it->value(), $top_prime, "iterator object can rewind to $top_prime");
  $it->next;
  is( $it->value(), $big_prime, "iterator object next is $big_prime");
  $it->rewind(~0);
  is( $it->value(), $big_prime, "iterator object rewound to ~0 is $big_prime");
  $it->prev;
  is( $it->value(), $top_prime, "iterator object prev goes back to $top_prime");

  # Validation for the Math::NumSeq compatiblity stuff
  $it->rewind;
  do { $it->next } for 1..200;
  is( $it->tell_i(), 201, "iterator object tell_i");
  is( $it->i_start, 1, "iterator object i_start = 1");
  like( $it->description, qr/prime numbers/, "iterator object description");
  is( $it->values_min, 2, "iterator object values_min = 2");
  is( $it->values_max, undef, "iterator object values_max = undef");
  # missing: characteristic
  is( $it->oeis_anum, "A000040", "iterator object oeis_anum = A000040");
  # missing: parameter_info_array / parameter_info_list
  is( $it->seek_to_i(156)->value, 911, "iterator object seek_to_i goes to nth prime");
  is( $it->seek_to_value(156)->value, 157, "iterator object seek_to_value goes to value");
  is( $it->ith(589), 4289, "iterator object ith returns nth prime");
  ok( $it->pred(577), "iterator object pred returns true if is_prime");
  is( $it->value_to_i(4289), 589, "iterator object value_to_i works");
  is( $it->value_to_i(4290), undef, "iterator object value_to_i for non-prime returns undef");
  is( $it->value_to_i_floor(4290), 589, "iterator object value_to_i_floor");
  is( $it->value_to_i_ceil(4290), 590, "iterator object value_to_i_ceil");
  my $est = $it->value_to_i_estimate( 4171510507 );
  my $act = 197710788;
  # We will get an estimate that is much, much closer than Math::NumSeq
  ok( ($est > ($act-500)) && ($est < ($act+500)),
      "iterator object value_to_i_estimage is in range");
}
