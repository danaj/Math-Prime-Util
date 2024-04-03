#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/cornacchia/;

#my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
#my $usexs  = Math::Prime::Util::prime_get_config->{'xs'};
#my $usegmp = Math::Prime::Util::prime_get_config->{'gmp'};

my @tests = (
  [0,0, [0,0]],
  [0,13, [undef]],
  [0,16, [4,0]],
  [1,113, [8,7]],
  [5,29, [3,2]],                 # Examples from NZMATH
  [7,4*29, [2,4]],
  [5,12829, [52,45]],            # Example from Buhler and Wagon
  [6,103, [7,3]],
  [1, 4*113, [16,14]],           # Some examples from Pari/GP
  [1, 4*109, [20,6]],
  [2,603, [21,9]],
  [24,"122838793181521", [10547339,694995]],
  [59551,100123456, [135,41]],   # Can find with sqrtmod
  [57564,100123456, [9934,5]],   # Can find with allsqrtmod
  [56892,100123456, [3016,40]],  # Can find with loop
);

plan tests => 0
            + scalar(@tests)
            ;

for my $data (@tests) {
  my($d,$n,$exp) = @$data;
  my @got = cornacchia($d,$n);
  is_deeply( \@got, $exp, defined $exp->[0]
             ? "Expected ($exp->[0],$exp->[1]) for x^2 + $d*y^2 = $n"
             : "Expected no solution for for x^2 + $d*y^2 = $n" );
}
