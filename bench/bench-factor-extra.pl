#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util qw/-nobigint/;
use Benchmark qw/:all/;
use List::Util qw/min max/;
use Config;
my $count = shift || -2;
my $is64bit = (~0 > 4294967295);
my $maxdigits = ($is64bit) ? 20 : 10;  # Noting the range is limited for max.

my $rgen = sub {
  my $range = shift;
  return 0 if $range <= 0;
  my $rbits = 0; { my $t = $range; while ($t) { $rbits++; $t >>= 1; } }
  while (1) {
    my $rbitsleft = $rbits;
    my $U = 0;
    while ($rbitsleft > 0) {
      my $usebits = ($rbitsleft > $Config{randbits}) ? $Config{randbits} : $rbitsleft;
      $U = ($U << $usebits) + int(rand(1 << $usebits));
      $rbitsleft -= $usebits;
    }
    return $U if $U <= $range;
  }
};

srand(29);
my $rounds = 400;
my $sqrounds = 256*1024;
my $rsqrounds = 32*1024;
my $p1smooth = 1000;
my $hrounds = 10000;
my $num_nums = 1000;
test_at_digits($_) for ( 3 .. $maxdigits );


sub test_at_digits {
  my $digits = shift;

  die "Digits has to be >= 1" unless $digits >= 1;
  die "Digits has to be <= $maxdigits" if $digits > $maxdigits;

  my @nums = genrand($digits, $num_nums);
  #my @nums = gensemi($digits, $num_nums, 23);
  my $min_num = min @nums;
  my $max_num = max @nums;

  # Determine success rates
  my %nfactored;
  my $tfac = 0;
  # Did we find any non-trivial factors?
  my $calc_nfacs = sub { ((scalar grep { $_ > 5 } @_) > 1) ? 1 : 0 };
  for (@nums) {
    $tfac += $calc_nfacs->(Math::Prime::Util::factor($_));
    $nfactored{'prho'} += $calc_nfacs->(Math::Prime::Util::prho_factor($_, $rounds));
    $nfactored{'pbrent'} += $calc_nfacs->(Math::Prime::Util::pbrent_factor($_, $rounds));
    $nfactored{'pminus1'} += $calc_nfacs->(Math::Prime::Util::pminus1_factor($_, $p1smooth));
    $nfactored{'pplus1'} += $calc_nfacs->(Math::Prime::Util::pplus1_factor($_, $p1smooth));
    $nfactored{'squfof'} += $calc_nfacs->(Math::Prime::Util::squfof_factor($_, $sqrounds));
    #$nfactored{'trial'} += $calc_nfacs->(Math::Prime::Util::trial_factor($_));
    #$nfactored{'fermat'} += $calc_nfacs->(Math::Prime::Util::fermat_factor($_, $rounds));
    $nfactored{'holf'} += $calc_nfacs->(Math::Prime::Util::holf_factor($_, $hrounds));
  }

  print "factoring $num_nums random $digits-digit numbers ($min_num - $max_num)\n";
  print "Factorizations: ",
         join(", ", map { sprintf "%s %4.1f%%", $_, 100*$nfactored{$_}/$tfac }
                    grep { $_ ne 'fermat' }
                    sort {$nfactored{$a} <=> $nfactored{$b}} keys %nfactored),
         "\n";

  my $lref = {
    "prho"    => sub { Math::Prime::Util::prho_factor($_, $rounds) for @nums },
    "pbrent"  => sub { Math::Prime::Util::pbrent_factor($_, $rounds) for @nums },
    "pminus1" => sub { Math::Prime::Util::pminus1_factor($_, $rounds) for @nums },
    "pplus1"  => sub { Math::Prime::Util::pplus1_factor($_, $rounds) for @nums},
    "fermat"  => sub { Math::Prime::Util::fermat_factor($_, $rounds) for @nums},
    "holf"    => sub { Math::Prime::Util::holf_factor($_, $hrounds) for @nums },
    "squfof"  => sub { Math::Prime::Util::squfof_factor($_, $sqrounds) for @nums },
    "trial"   => sub { Math::Prime::Util::trial_factor($_) for @nums },
  };
  delete $lref->{'fermat'} if $digits >= 9;
  delete $lref->{'holf'} if $digits >= 17;
  delete $lref->{'trial'} if $digits >= 15;
  cmpthese($count, $lref);
  print "\n";
}


sub genrand {
  my $digits = shift;
  my $num = shift;

  my $base = ($digits == 1) ? 0 : int(10 ** ($digits-1));
  my $max = int(10 ** $digits);
  $max = ~0 if $max > ~0;
  my @nums = map { $base + $rgen->($max-$base) } (1 .. $num);
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
      $n = $base + $rgen->($max-$base);
      $n += (1,0,5,4,3,2,1,0,3,2,1,0,1,0,3,2,1,0,1,0,3,2,1,0,5,4,3,2,1,0)[$n%30];
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
