#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/valuation remove_factors remove_factors_exp/;

#my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
#my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
#my $usegmp= Math::Prime::Util::prime_get_config->{'gmp'};
#my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
#$use64 = 0 if $use64 && 18446744073709550592 == ~0;

plan tests => 2;

subtest 'valuation' => sub {
  my @valuations = (
    [-4,2, 2],
    #[0,0, 0],  error
    #[1,0, 0],
    [0,2, undef],
    [1,2, 0],
    [96552,6, 3],
    [1879048192,2, 28],
    ["65520150907877741108803406077280119039314703968014509493068998974809747144832",2, 7],
  );
  foreach my $r (@valuations) {
    my($n, $k, $exp) = @$r;
    is(valuation($n, $k), $exp, "valuation($n,$k) = ".(defined($exp)?$exp:"<undef>"));
  }
};

subtest 'remove_factors' => sub {
  my @removes = (
    [ 72,  6,  2,  2],
    [ 72, 12,  6,  1],
    [ 72,  8,  9,  1],
    [  5,  6,  5,  0],
    [-72,  6, -2,  2],
    [-12,  3, -4,  1],
    ["184467440737095516160", 10, "18446744073709551616", 1],
    ["-184467440737095516160", 10, "-18446744073709551616", 1],
  );
  foreach my $r (@removes) {
    my($n, $k, $rem, $exp) = @$r;
    is("".remove_factors($n,$k), "$rem", "remove_factors($n,$k) = $rem");
    is_deeply([map { defined($_) ? "$_" : undef } remove_factors_exp($n,$k)],
              ["$rem","$exp"],
              "remove_factors_exp($n,$k) = ($rem,$exp)");
  }
  is(remove_factors(0,6), undef, "remove_factors(0,6) = undef");
  is_deeply([remove_factors_exp(0,6)], [undef,undef], "remove_factors_exp(0,6) = (undef,undef)");
  ok(!eval { remove_factors(10,1); 1 }, "remove_factors rejects k=1");
  ok(!eval { remove_factors_exp(10,0); 1 }, "remove_factors_exp rejects k=0");
};
