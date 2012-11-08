#!/usr/bin/env perl
use strict;
use warnings;
$| = 1;  # fast pipes
srand(377);

use Math::Prime::Util qw/factor -nobigint/;
use Math::Factor::XS qw/prime_factors/;
use Math::Pari qw/factorint/;
use Benchmark qw/:all/;
use Data::Dumper;
use Config;
my $digits = shift || 15;
my $count = shift || -3;

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

my @min_factors_by_digit = (2,2,3,3,5,11,17,47,97);
my $smallest_factor_allowed = $min_factors_by_digit[$digits];
$smallest_factor_allowed = $min_factors_by_digit[-1] unless defined $smallest_factor_allowed;
my $numprimes = 200;

die "Digits has to be >= 2" unless $digits >= 2;
die "Digits has to be <= 10" if (~0 == 4294967295) && ($digits > 10);
die "Digits has to be <= 19" if $digits > 19;

my $skip_mfxs = ($digits > 17);

# Construct some semiprimes of the appropriate number of digits
# There are much cleverer ways of doing this, using randomly selected
# nth_primes, and so on, but this works well until we get lots of digits.
print "Generating $numprimes random $digits-digit semiprimes (min factor $smallest_factor_allowed) ";
my @semiprimes;
foreach my $i ( 1 .. $numprimes ) {
  my $base = int(10 ** ($digits-1));
  my $add = int(10 ** ($digits)) - $base;
  my @factors;
  my $n;
  while (1) {
    $n = $base + $rgen->($add);
    next if $n > (~0 - 4);
    $n += (1,0,5,4,3,2,1,0,3,2,1,0,1,0,3,2,1,0,1,0,3,2,1,0,5,4,3,2,1,0)[$n%30];
    @factors = factor($n);
    next if scalar @factors != 2;
    next if $factors[0] < $smallest_factor_allowed;
    next if $factors[1] < $smallest_factor_allowed;
    last if scalar @factors == 2;
  }
  die "ummm... $n != $factors[0] * $factors[1]\n" unless $n == $factors[0] * $factors[1];
  #print "$n == $factors[0] * $factors[1]\n";
  push @semiprimes, $n;
  print "." if ($i % ($numprimes/10)) == 0;
}
print "done.\n";

print "Verifying Math::Prime::Util $Math::Prime::Util::VERSION ...";
foreach my $sp (@semiprimes) {
  my @factors = factor($sp);
  die "wrong for $sp\n" unless ($#factors == 1) && ($factors[0] * $factors[1]) == $sp;
}
print "OK\n";
if (!$skip_mfxs) {
  print "Verifying Math::Factor::XS $Math::Factor::XS::VERSION ...";
  foreach my $sp (@semiprimes) {
    my @factors = prime_factors($sp);
    die "wrong for $sp\n" unless ($#factors == 1) && ($factors[0] * $factors[1]) == $sp;
  }
  print "OK\n";
} else {
  print "Math::Factor::XS is too slow for $digits digits.  Skipping.\n";
}
print "Verifying Math::Pari $Math::Pari::VERSION ...";
foreach my $sp (@semiprimes) {
  my @factors;
  my ($pn,$pc) = @{factorint($sp)};
  push @factors, (int($pn->[$_])) x $pc->[$_] for (0 .. $#{$pn});
  die "wrong for $sp\n" unless ($#factors == 1) && ($factors[0] * $factors[1]) == $sp;
}
print "OK\n";

my %compare = (
    'MPU'   => sub { factor($_) for @semiprimes; },
    'MFXS'  => sub { prime_factors($_) for @semiprimes; },
    'Pari'  => sub { foreach my $n (@semiprimes) { my @factors; my ($pn,$pc) = @{factorint($n)}; push @factors, (int($pn->[$_])) x $pc->[$_] for (0 .. $#{$pn}); } }
);
delete $compare{'MFXS'} if $skip_mfxs;

cmpthese($count, \%compare);
