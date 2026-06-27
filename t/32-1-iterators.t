#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/prime_iterator prime_iterator_object/;
use Math::BigInt try => "GMP,GMPz";

my $use64 = 4294967295 < ~0;
my $broken64 = (18446744073709550592 == ~0);

subtest 'prime_iterator' => sub {
  ok(!eval { prime_iterator(-2); }, "iterator -2");
  ok(!eval { prime_iterator("abc"); }, "iterator abc");
  ok(!eval { prime_iterator(4.5); }, "iterator 4.5");

  my $it1 = prime_iterator();
  is_deeply( [map { $it1->() } 1..10], [2,3,5,7,11,13,17,19,23,29], "iterator first 10 primes" );

  my @T = (
    [   47, 5, [47,53,59,61,67]],
    [  199, 3, [199,211,223]],
    [  200, 3, [211,223,227]],
    [31397, 3, [31397,31469,31477]],
    [31396, 3, [31397,31469,31477]],
    [31398, 3, [31469,31477,31481]],
  );

  for my $i (@T) {
    my($start, $count, $exp) = @$i;
    my $it = prime_iterator($start);
    my $got = [map { $it->() } 1..$count];
    is_deeply($got, $exp, "start at $start, iterate $count times");
  }

  my $it2 = prime_iterator(Math::BigInt->new("68719476736"));
  is_deeply( [map { $it2->() } 1..3],
             [68719476767,68719476851,68719476853],
             "iterator 3 primes with BigInt start" );

  my @bgot;
  my $ita = prime_iterator();
  while ((my $a = $ita->()) <= 10) {
    my $itb = prime_iterator(10*$a);
    while ((my $b = $itb->()) <= 10*$a+10) {
      my $itc = prime_iterator($b);
      while ((my $c = $itc->()) <= $b+10) {
        push @bgot, $c;
      }
    }
  }
  is_deeply( \@bgot, [qw/23 29 31 29 31 37 31 37 41 37 41 43 47 53 59 61 59 61 67 71 73 79 73 79 83 79 83 89/], "triple nested iterator" );
};

subtest 'prime_iterator_object basic iteration' => sub {
  ok(!eval { prime_iterator_object(-2); }, "iterator -2");
  ok(!eval { prime_iterator_object("abc"); }, "iterator abc");
  ok(!eval { prime_iterator_object(4.5); }, "iterator 4.5");

  my $it1 = prime_iterator_object();
  is_deeply( [map { $it1->iterate() } 1..10], [2,3,5,7,11,13,17,19,23,29], "iterator first 10 primes" );

  my @T = (
    [   47, 5, [47,53,59,61,67]],
    [  199, 3, [199,211,223]],
    [  200, 3, [211,223,227]],
    [31397, 3, [31397,31469,31477]],
    [31396, 3, [31397,31469,31477]],
    [31398, 3, [31469,31477,31481]],
  );

  for my $i (@T) {
    my($start, $count, $exp) = @$i;
    my $it = prime_iterator_object($start);
    my $got = [map { $it->iterate() } 1..$count];
    is_deeply($got, $exp, "start at $start, iterate $count times");
  }
};

subtest 'prime_iterator_object methods' => sub {
  my $it = prime_iterator_object;
  do { $it->next } for 1..10;
  is( $it->value(), 31, "iterator object moved forward 10 now returns 31");
  $it->prev;
  is( $it->value(), 29, "iterator object moved back now returns 29");
  is( $it->peek(), 31, "iterator object peek shows 31");
  is( $it->iterate(), 29, "iterator object iterates to 29");
  is( $it->iterate(), 31, "iterator object iterates to 31");
  $it->rewind->next->next->next->prev;
  is( $it->value(), 5, "iterator object rewind and move returns 5");
  $it->rewind(1);
  is( $it->value(), 2, "iterator object rewind(1) goes to 2");
  $it->rewind(0);
  is( $it->value(), 2, "iterator object rewind(0) goes to 2");

  SKIP: {
    skip "This Perl can't add large numbers correctly",4 if $broken64;
    my $top_prime = $use64 ? "18446744073709551557" : "4294967291";
    my $big_prime = $use64 ? "18446744073709551629" : "4294967311";
    $it->rewind($top_prime);
    is( $it->value(), $top_prime, "iterator object can rewind to $top_prime");
    $it->next;
    is( "".$it->value(), $big_prime, "iterator object next is $big_prime");
    $it->rewind(~0);
    is( "".$it->value(), $big_prime, "iterator object rewound to ~0 is $big_prime");
    $it->prev;
    is( $it->value(), $top_prime, "iterator object prev goes back to $top_prime");
  }

  $it->rewind;
  do { $it->next } for 1..100;
  is( $it->tell_i(), 101, "iterator object tell_i");
  is( $it->i_start, 1, "iterator object i_start = 1");
  like( $it->description, qr/prime numbers/, "iterator object description");
  is( $it->values_min, 2, "iterator object values_min = 2");
  is( $it->values_max, undef, "iterator object values_max = undef");
  is( $it->oeis_anum, "A000040", "iterator object oeis_anum = A000040");
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
  ok( ($est > ($act-500)) && ($est < ($act+500)),
      "iterator object value_to_i_estimate is in range");
};

done_testing();
