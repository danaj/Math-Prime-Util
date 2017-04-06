#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util;
use Math::Prime::Util::GMP;
use Math::Primality;
use Benchmark qw/:all/;
my $count = shift || -2;
srand(29);  # So we have repeatable results
Math::Prime::Util::srand(29);

test_at_digits($_, 1000) for (5, 15, 25, 50, 200);

sub test_at_digits {
  my($digits, $numbers) = @_;
  die "Digits must be > 0" unless $digits > 0;

  my $start = Math::Prime::Util::random_ndigit_prime($digits) - 3;
  my $end = $start;
  $end = Math::Prime::Util::GMP::next_prime($end) for 1 .. $numbers;

  print "next_prime x $numbers starting at $start\n";

  cmpthese($count,{
    'MP'       => sub { my $n = $start; $n = Math::Primality::next_prime($n) for 1..$numbers; die "MP ended with $n instead of $end" unless $n == $end; },
    'MPU'      => sub { my $n = $start; $n = Math::Prime::Util::next_prime($n) for 1..$numbers; die "MPU ended with $n instead of $end" unless $n == $end; },
    'MPU GMP'  => sub { my $n = $start; $n = Math::Prime::Util::GMP::next_prime($n) for 1..$numbers; die "MPU GMP ended with $n instead of $end" unless $n == $end; },
  });
}
