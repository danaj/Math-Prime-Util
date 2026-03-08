#!/usr/bin/env perl
use strict;
use warnings;

BEGIN {
  unless ($ENV{RELEASE_TESTING}) {
    require Test::More;
    Test::More::plan(skip_all => 'these tests are for release candidate testing');
  }
}

#---------------------------------------------------------------------


use Test::More;
eval "use Test::Weaken";
plan skip_all => "Test::Weaken required for testing leaks" if $@;

use Math::Prime::Util qw/prime_iterator primes factor moebius
                         forprimes divisors divrem/;
use Math::Prime::Util::PrimeIterator;

my $leaks;

$leaks = Test::Weaken::leaks(
  { constructor => sub {
      my $it = prime_iterator(100);
      return \$it;
    }
  });
ok(!$leaks, "Prime iterator doesn't leak");

$leaks = Test::Weaken::leaks(
  { constructor => sub {
      my $it = Math::Prime::Util::PrimeIterator->new(100);
      return \$it;
    }
  });
ok(!$leaks, "Prime iterator object doesn't leak");

$leaks = Test::Weaken::leaks(
  { constructor => sub {
      my $it = prime_iterator("1000000000000000000000000");
      my $p = $it->();
      return \$it;
    }
  });
ok(!$leaks, "Prime iterator object doesn't leak");

$leaks = Test::Weaken::leaks(
  sub {
    my $primes = primes(66000);
    return $primes
  });
ok(!$leaks, "primes array doesn't leak");

$leaks = Test::Weaken::leaks(
  sub {
    return [factor(2938424874)];
  });
ok(!$leaks, "factor array doesn't leak");

$leaks = Test::Weaken::leaks(
  sub {
    return [moebius(500,1000)];
  });
ok(!$leaks, "moebius range doesn't leak");

$leaks = Test::Weaken::leaks(
  sub {
    my @p;
    forprimes { push @p, $_ } 1000;
    return \@p;
  });
ok(!$leaks, "forprimes block doesn't leak");

$leaks = Test::Weaken::leaks(
  sub {
    return [divrem(1000000007, 12345)];
  });
ok(!$leaks, "divrem doesn't leak");

$leaks = Test::Weaken::leaks(
  sub {
    return [divisors(720720)];
  });
ok(!$leaks, "divisors doesn't leak");

done_testing();
