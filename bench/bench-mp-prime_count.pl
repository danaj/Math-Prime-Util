#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util;
use Math::Prime::Util::GMP;
use Math::Prime::Util::PP;
use Math::Primality;
use Benchmark qw/:all/;
my $count = shift || -2;

#my($n, $exp) = (100000,9592);
my($n, $exp) = (1000000,78498);
#my($n, $exp) = (10000000,664579);
cmpthese($count,{
  'MP'           =>sub { die unless $exp == Math::Primality::prime_count($n); },
  'MPU default'  =>sub { die unless $exp == Math::Prime::Util::prime_count($n); },
  'MPU XS Sieve' =>sub { die unless $exp == Math::Prime::Util::_XS_segment_pi($n); },
  #'MPU XS Lehmer'=>sub { die unless $exp == Math::Prime::Util::_XS_lehmer_pi($n); },
  'MPU PP Sieve' =>sub { die unless $exp == Math::Prime::Util::PP::_sieve_prime_count($n); },
  'MPU PP Lehmer'=>sub { die unless $exp == Math::Prime::Util::PP::_lehmer_pi($n); },
  'MPU GMP Trial'=>sub { die unless $exp == Math::Prime::Util::GMP::prime_count(2,$n); },
});
