#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util qw/prime_precalc prime_memfree/;

use Test::More  tests => 3 + 3 + 3 + 6;


my $bigsize = 10_000_000;

# How do we test these?  Use some private functions and assume things about
# internal behavior.  Not the best idea.

can_ok( 'Math::Prime::Util', '_get_prime_cache_size' );

my $init_size = Math::Prime::Util::_get_prime_cache_size;

prime_precalc($bigsize);

cmp_ok( Math::Prime::Util::_get_prime_cache_size, '>', $init_size, "Internal space grew after large precalc" );

prime_memfree;

is( Math::Prime::Util::_get_prime_cache_size, $init_size, "Internal space went back to original size after memfree" );


# Now do the object way.
{
  #my $mf = new_ok( 'Math::Prime::Util::MemFree');  # Better 0.88+ way
  my $mf = Math::Prime::Util::MemFree->new;
  isa_ok $mf, 'Math::Prime::Util::MemFree';
  prime_precalc($bigsize);
  cmp_ok( Math::Prime::Util::_get_prime_cache_size, '>', $init_size, "Internal space grew after large precalc" );
}
is( Math::Prime::Util::_get_prime_cache_size, $init_size, "Memory released after MemFree object goes out of scope");

# Wrap multiple calls, make sure we wait until the last one is done.
{
  my $mf = Math::Prime::Util::MemFree->new;
  prime_precalc($bigsize);
  cmp_ok( Math::Prime::Util::_get_prime_cache_size, '>', $init_size, "Internal space grew after large precalc" );
  {
    my $mf2 = Math::Prime::Util::MemFree->new;
    prime_precalc( 2 * $bigsize );
  }
  cmp_ok( Math::Prime::Util::_get_prime_cache_size, '>', $init_size, "Memory not freed yet because a MemFree object still live." );
}
is( Math::Prime::Util::_get_prime_cache_size, $init_size, "Memory released after last MemFree object goes out of scope");

# Show how an eval death can leak
eval { prime_precalc($bigsize); cmp_ok( Math::Prime::Util::_get_prime_cache_size, '>', $init_size, "Internal space grew after large precalc" ); prime_memfree; };
is( Math::Prime::Util::_get_prime_cache_size, $init_size, "Memory freed after successful eval");

eval { prime_precalc($bigsize); cmp_ok( Math::Prime::Util::_get_prime_cache_size, '>', $init_size, "Internal space grew after large precalc" ); die; prime_memfree; };
isnt( Math::Prime::Util::_get_prime_cache_size, $init_size, "Memory normally not freed after eval die");
prime_memfree;

eval { my $mf = Math::Prime::Util::MemFree->new; prime_precalc($bigsize); cmp_ok( Math::Prime::Util::_get_prime_cache_size, '>', $init_size, "Internal space grew after large precalc" ); die; };
is( Math::Prime::Util::_get_prime_cache_size, $init_size, "Memory is freed after eval die using object scoper");
