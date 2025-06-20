#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util qw/prime_precalc prime_memfree prime_get_config/;
use Math::Prime::Util::MemFree;

use Test::More  tests => 3 + 3 + 3 + 6;


my $bigsize = 10_000_000;

# This is still a slightly dubious assumption, that the precalc size _must_
# go up when we request it.

can_ok( 'Math::Prime::Util', 'prime_get_config' );

my $diag = "Using " .
  ((Math::Prime::Util::prime_get_config->{xs})
    ? "XS"
    : "PP") .
  ((Math::Prime::Util::prime_get_config->{gmp})
    ? " with MPU::GMP version $Math::Prime::Util::GMP::VERSION."
    : ".") .
  "\n";
diag $diag;

my $init_size = prime_get_config->{'precalc_to'};

prime_precalc($bigsize);

cmp_ok( prime_get_config->{'precalc_to'}, '>', $init_size, "Internal space grew after large precalc" );

prime_memfree;

is( prime_get_config->{'precalc_to'}, $init_size, "Internal space went back to original size after memfree" );


# Now do the object way.
{
  #my $mf = new_ok( 'Math::Prime::Util::MemFree');  # Better 0.88+ way
  my $mf = Math::Prime::Util::MemFree->new;
  isa_ok $mf, 'Math::Prime::Util::MemFree';
  prime_precalc($bigsize);
  cmp_ok( prime_get_config->{'precalc_to'}, '>', $init_size, "Internal space grew after large precalc" );
}
is( prime_get_config->{'precalc_to'}, $init_size, "Memory released after MemFree object goes out of scope");

# Wrap multiple calls, make sure we wait until the last one is done.
{
  my $mf = Math::Prime::Util::MemFree->new;
  prime_precalc($bigsize);
  cmp_ok( prime_get_config->{'precalc_to'}, '>', $init_size, "Internal space grew after large precalc" );
  {
    my $mf2 = Math::Prime::Util::MemFree->new;
    prime_precalc( 2 * $bigsize );
  }
  cmp_ok( prime_get_config->{'precalc_to'}, '>', $init_size, "Memory not freed yet because a MemFree object still live." );
}
is( prime_get_config->{'precalc_to'}, $init_size, "Memory released after last MemFree object goes out of scope");

# Show how an eval death can leak
eval { prime_precalc($bigsize); cmp_ok( prime_get_config->{'precalc_to'}, '>', $init_size, "Internal space grew after large precalc" ); prime_memfree; };
is( prime_get_config->{'precalc_to'}, $init_size, "Memory freed after successful eval");

eval { prime_precalc($bigsize); cmp_ok( prime_get_config->{'precalc_to'}, '>', $init_size, "Internal space grew after large precalc" ); die; prime_memfree; };
isnt( prime_get_config->{'precalc_to'}, $init_size, "Memory normally not freed after eval die");
prime_memfree;

eval { my $mf = Math::Prime::Util::MemFree->new; prime_precalc($bigsize); cmp_ok( prime_get_config->{'precalc_to'}, '>', $init_size, "Internal space grew after large precalc" ); die; };
is( prime_get_config->{'precalc_to'}, $init_size, "Memory is freed after eval die using object scoper");
