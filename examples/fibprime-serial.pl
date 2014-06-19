#!/usr/bin/env perl
use strict;
use warnings;

# Overkill, but let's try to select a good bigint module.
my $bigint_class;
if      (eval { require Math::GMPz; 1; }) {
  $bigint_class = "Math::GMPz";
} elsif (eval { require Math::GMP; 1; }) {
  $bigint_class = "Math::GMP";
} else {
  require Math::BigInt;
  Math::BigInt->import(try=>"GMP,Pari");
  $bigint_class = "Math::BigInt";
}

use Math::Prime::Util ':all';
use Time::HiRes qw(gettimeofday tv_interval);
$| = 1;

my $time_start = [gettimeofday];
prime_precalc(1_000_000);


{
  my @fibstate;
  my $nth = 1;
  my $n = 0;
  while (1) {
    # Exploit knowledge that excepting k=4, all prime F_k have a prime k.
    my $k = ($nth <= 2) ?  2 + $nth  :  nth_prime($nth);
    $nth++;
    my $Fk = fib_n($k, \@fibstate);
    if (is_prob_prime($Fk)) {
      my $time_int = tv_interval($time_start);
      printf "%3d %7d %20.5f\n", ++$n, $k, $time_int;
    }
  }
}


sub fib_n {
  my ($n, $fibstate) = @_;
  @$fibstate = (1, $bigint_class->new(0), $bigint_class->new(1))
     unless defined $fibstate->[0];
  my ($curn, $a, $b) = @$fibstate;
  die "fib_n only increases" if $n < $curn;
  do { ($a, $b) = ($b, $a+$b); } for (1 .. $n-$curn);
  @$fibstate = ($n, $a, $b);
  $b;
}

