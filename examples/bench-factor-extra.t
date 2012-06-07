#!/usr/bin/perl -w

use strict;
use Math::Prime::Util;
use Benchmark qw/:all/;
use List::Util qw/min max/;
my $count = shift || -2;

srand(29);
my $rounds = 100;
my $sqrounds = 64*1024;
test_at_digits($_) for (3..10);


sub test_at_digits {
  my $digits = shift;

  die "Digits has to be >= 1" unless $digits >= 1;
  die "Digits has to be <= 10" if (~0 == 4294967295) && ($digits > 10);
  die "Digits has to be <= 19" if $digits > 19;

  my @nums = genrand($digits, 1000);
  #my @nums = gensemi($digits, 1000, 23);
  my $min_num = min @nums;
  my $max_num = max @nums;

  # Determine success rates
  my %nfactored;
  for (@nums) {
    $nfactored{'prho'} += int((scalar grep { $_ > 5 } Math::Prime::Util::prho_factor($_, $rounds))/2);
    $nfactored{'pbrent'} += int((scalar grep { $_ > 5 } Math::Prime::Util::pbrent_factor($_, $rounds))/2);
    $nfactored{'pminus1'} += int((scalar grep { $_ > 5 } Math::Prime::Util::pminus1_factor($_, $rounds))/2);
    $nfactored{'fermat'} += int((scalar grep { $_ > 5 } Math::Prime::Util::fermat_factor($_, $rounds))/2);
    $nfactored{'squfof'} += int((scalar grep { $_ > 5 } Math::Prime::Util::squfof_factor($_, $sqrounds))/2);
    #$nfactored{'trial'} += int((scalar grep { $_ > 5 } Math::Prime::Util::trial_factor($_))/2);
  }
  my $tfac = $nfactored{'fermat'};

  print "factoring 1000 random $digits-digit numbers ($min_num - $max_num)\n";
  print "Factorizations: ",
         join(", ", map { sprintf "%s %4.1f%%", $_, 100*$nfactored{$_}/$tfac }
                    grep { $_ ne 'fermat' }
                    sort {$nfactored{$a} <=> $nfactored{$b}} keys %nfactored),
         "\n";

  my $lref = {
    "prho"    => sub { Math::Prime::Util::prho_factor($_, $rounds) for @nums },
    "pbrent"  => sub { Math::Prime::Util::pbrent_factor($_, $rounds) for @nums },
    "pminus1" => sub { Math::Prime::Util::pminus1_factor($_, $rounds) for @nums },
    "fermat"  => sub { Math::Prime::Util::fermat_factor($_, $rounds) for @nums },
    "squfof"  => sub { Math::Prime::Util::squfof_factor($_, $sqrounds) for @nums },
    "trial"   => sub { Math::Prime::Util::trial_factor($_) for @nums },
  };
  delete $lref->{'fermat'} if $digits >= 9;
  cmpthese($count, $lref);
  print "\n";
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
  my $smallest_factor = shift;

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
