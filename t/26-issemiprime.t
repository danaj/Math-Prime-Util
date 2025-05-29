#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_semiprime/;
my $usegmp = Math::Prime::Util::prime_get_config->{'gmp'};

plan tests => 6;

{
  is_deeply( [map { is_semiprime($_) } 121,341,2047,7009,28153],
             [qw/1 1 1 1 1/],
             "Semiprimes that were incorrectly calculated in v0.70" );
  is_deeply( [grep { is_semiprime($_) } 10000..10100],
             [qw/10001 10003 10006 10015 10018 10019 10021 10022 10027 10029 10031 10033 10041 10042 10046 10049 10055 10057 10063 10073 10077 10078 10081 10083 10085 10097/],
             "Identify semiprimes from 10000 to 10100" );
  is(is_semiprime("669386384129397653"), 1, "is_semiprime(669386384129397653)");
  is(is_semiprime("10631816576169524657"), 1, "is_semiprime(10631816576169524657)");
  is(is_semiprime("1814186289136250293214268090047441301"), 0, "is_semiprime(1814186289136250293214268090047441301)");
  SKIP: {
    skip "Skipping difficult is_semiprime", 1 unless $usegmp;
    is(is_semiprime("42535430147496493121551759"), 0, "is_semiprime(42535430147496493121551759)");
  }
}
