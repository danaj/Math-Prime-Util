#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util qw/factor/;
use Math::Factor::XS qw/prime_factors/;
use Benchmark qw/:all/;
use List::Util qw/min max reduce/;
my $count = shift || -2;
my $is64bit = (~0 > 4294967295);
my $maxdigits = ($is64bit) ? 20 : 10;  # Noting the range is limited for max.
my $semiprimes = 0;
my $howmany = 10000;

srand(29);

for my $d ( 3 .. $maxdigits ) {
  print "Factor $howmany $d-digit numbers\n";
  test_at_digits($d, $howmany);
}

sub test_at_digits {
  my $digits = shift;
  die "Digits has to be >= 1" unless $digits >= 1;
  die "Digits has to be <= $maxdigits" if $digits > $maxdigits;
  my $quantity = shift;

  my @rnd = genrand($digits, $quantity);
  my @smp = gensemi($digits, $quantity);

  # verify
  foreach my $p (@rnd, @smp) {
    my @mpxs = prime_factors($p);  push @mpxs, $p if $p < 2;

    verify_factor($p, \@mpxs, [factor($p)], "Math::Prime::Util $Math::Prime::Util::VERSION");
  }


  #my $min_num = min @nums;
  #my $max_num = max @nums;
  #my $whatstr = "$digits-digit ", $semiprimes ? "semiprime" : "random";
  #print "factoring 1000 $digits-digit ",
  #      $semiprimes ? "semiprimes" : "random numbers",
  #      " ($min_num - $max_num)\n";

  my $lref = {
    "MPU  rand"  => sub { factor($_) for @rnd },
    "MPU  semi"  => sub { factor($_) for @smp },
    "MFXS rand"  => sub { prime_factors($_) for @rnd },
    "MFXS semi"  => sub { prime_factors($_) for @smp },
  };
  cmpthese($count, $lref);
}

sub verify_factor {
  my $n = shift;
  my $aref_master = shift;
  my $aref_check = shift;
  my $name = shift;

  my @master = sort {$a<=>$b} @{$aref_master};
  my @check  = sort {$a<=>$b} @{$aref_check};
  die "Factor $n master fail!" unless $n == reduce { $a * $b } @master;
  die "Factor $n fail: $name" unless $#check == $#master;
  die "Factor $n fail: $name" unless $n == reduce { $a * $b } @check;
  for (0 .. $#master) {
    die "Factor $n fail: $name" unless $master[$_] == $check[$_];
  }
  1;
}

sub genrand {
  my $digits = shift;
  my $num = shift;

  my $base = ($digits == 1) ? 0 : int(10 ** ($digits-1));
  my $max = int(10 ** $digits);
  $max = ~0 if $max > ~0;
  my @nums = map { $base+int(rand($max-$base)) } (1 .. $num);
  return @nums;
}

sub gensemi {
  my $digits = shift;
  my $num = shift;

  my @min_factors_by_digit = (2,2,3,5,7,13,23,47,97);
  my $smallest_factor = $min_factors_by_digit[$digits];
  $smallest_factor = $min_factors_by_digit[-1] unless defined $smallest_factor;

  my $base = ($digits == 1) ? 0 : int(10 ** ($digits-1));
  my $max = int(10 ** $digits);
  $max = (~0-4) if $max > (~0-4);
  my @semiprimes;

  foreach my $i (1 .. $num) {
    my @factors;
    my $n;
    while (1) {
      $n = $base + int(rand($max-$base));
      $n += 1 if ($n%2) == 0;
      $n += 3 if ($n%3) == 0;
      @factors = Math::Prime::Util::factor($n);
      next if scalar @factors != 2;
      next if $factors[0] < $smallest_factor;
      next if $factors[1] < $smallest_factor;
      last if scalar @factors == 2;
    }
    die "ummm... $n != $factors[0] * $factors[1]\n" unless $n == $factors[0] * $factors[1];
    push @semiprimes, $n;
  }
  return @semiprimes;
}
