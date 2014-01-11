#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util qw/factor/;
# Compare to Math::Factor::XS, which uses trial division.
use Math::Factor::XS qw/prime_factors/;

use Benchmark qw/:all/;
use List::Util qw/min max reduce/;
my $count = shift || -2;
my $is64bit = (~0 > 4294967295);
my $maxdigits = ($is64bit) ? 20 : 10;  # Noting the range is limited for max.
my $semiprimes = 0;
my $howmany = 1000;

for my $d ( 3 .. $maxdigits ) {
  print "Factor $howmany $d-digit numbers\n";
  test_at_digits($d, $howmany);
}

sub test_at_digits {
  my $digits = shift;
  die "Digits has to be >= 1" unless $digits >= 1;
  die "Digits has to be <= $maxdigits" if $digits > $maxdigits;
  my $quantity = shift;

  my @rnd = ndigit_rand($digits, $quantity);
  my @smp = genrough($digits, $quantity);

  # verify (can be _really_ slow for 18+ digits)
  foreach my $p (@rnd, @smp) {
    next if $p < 2;
    verify_factor($p, [prime_factors($p)], [factor($p)], "Math::Prime::Util $Math::Prime::Util::VERSION");
  }

  #my $min_num = min @nums;
  #my $max_num = max @nums;
  #my $whatstr = "$digits-digit ", $semiprimes ? "semiprime" : "random";
  #print "factoring 1000 $digits-digit ",
  #      $semiprimes ? "semiprimes" : "random numbers",
  #      " ($min_num - $max_num)\n";

  my $lref = {
    "MPU   random"    => sub { my@a=factor($_) for @rnd },
    "MPU   nonsmooth" => sub { my@a=factor($_) for @smp },
    "MFXS  random"    => sub { my@a=prime_factors($_) for @rnd },
    "MFXS  nonsmooth" => sub { my@a=prime_factors($_) for @smp },
  };
  cmpthese($count, $lref);
}

sub verify_factor {
  my ($n, $aref1, $aref2, $name) = @_;

  return 1 if "@$aref1" eq "@$aref2";

  my @master = @$aref1;
  my @check  = @$aref2;
  die "Factor $n master fail!" unless $n == reduce { $a * $b } @master;
  die "Factor $n fail: $name" unless $#check == $#master;
  die "Factor $n fail: $name" unless $n == reduce { $a * $b } @check;
  for (0 .. $#master) {
    die "Factor $n fail: $name" unless $master[$_] == $check[$_];
  }
  1;
}

sub genrough {
  my ($digits, $num) = @_;

  my @min_factors_by_digit = (2,2,3,5,7,13,23,47,97);
  my $smallest_factor = $min_factors_by_digit[$digits];
  $smallest_factor = $min_factors_by_digit[-1] unless defined $smallest_factor;

  my @semiprimes;
  foreach my $i (1 .. $num) {
    my $n;
    my @facn;
    do {
      $n = ndigit_rand($digits, 1);
      @facn = Math::Prime::Util::trial_factor($n,$smallest_factor);
    } while scalar(@facn) > 1;
    push @semiprimes, $n;
  }
  return @semiprimes;
}

use Bytes::Random::Secure qw/random_string_from/;
sub ndigit_rand {
  my($digits, $howmany) = @_;
  die "digits must be > 0" if $digits < 1;
  $howmany = 1 unless defined $howmany;
  # TODO: need to skip things larger than ~0 for this module
  my @nums = map { random_string_from("123456789",1) . random_string_from("0123456789",$digits-1) } 1 .. $howmany;
  if (10**$digits > ~0) {  @nums = map { Math::BigInt->new($_) } @nums;  }
  else                  {  @nums = map { int($_) } @nums;                }
  return wantarray ? @nums : $nums[0];
}
