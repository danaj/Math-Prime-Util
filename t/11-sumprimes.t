#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/sum_primes vecsum primes/;

my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};

my %sums = (
  "189695660 to 189695892" => 0,
  "0 to 300000" => 3709507114,
  "12345 to 54321" => 132980191,
  "10000000 to 10001000" => 610034659,
);

my @large = (
  ["100", "1060"],
  ["1000", "76127"],
  ["10000", "5736396"],
  ["65535", "202288087"],
  ["65536", "202288087"],
  ["65537", "202353624"],
  ["321059", "4236201628"],
  ["321060", "4236201628"],
  ["321072", "4236201628"],
  ["321073", "4236522701"],
  ["1000000", "37550402023"],
  ["10000000", "3203324994356"],
  ["100000000", "279209790387276"],
  ["1000000000", "24739512092254535"],
  ["2000000000", "95673602693282040"],
  ["4000000000", "370412807102643725"],
  ["4294967295", "425649736193687430"],
  ["4294967296", "425649736193687430"],
  ["5000000000", "572840944428163514"],
  ["10000000000", "2220822432581729238"],
  ["19505444491", "8205714352685564257"],
  ["20705000000", "9222267677559724371"],
  ["29505444490", "18446744057541225032"],
  ["29505444491", "18446744087046669523"],
);

@large = grep { ($_->[0] <= 1000000) ||
                ($extra && ($usexs || $_->[0] <= 100000000))
              } @large;

plan tests => 1 + scalar(keys %sums) + scalar(@large);

{
  my @sum;
  my @exp;
  for (0..1000) {
    push @exp, vecsum( @{primes($_)} );
    push @sum, sum_primes($_);
  }
  is_deeply( \@sum, \@exp, "sum_primes for 0 to 1000" );
}
while (my($range, $expect) = each (%sums)) {
  my($low,$high) = $range =~ /(\d+) to (\d+)/;
  is( sum_primes($low,$high), $expect, "sum primes from $low to $high" );
}

foreach my $pair (@large) {
  my($n,$sum) = @$pair;
  is( "".sum_primes(0,$n), $sum, "sum_primes($n) = $sum" );
}
