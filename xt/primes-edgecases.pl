#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util ':all';
use Test::More;

my @primes = qw/2 3 5 7 11 13 17 19 23 29 31 37 41 43 47 53 59 61 67 71 73 79 83 89 97/;
my $end = 20;

plan tests => 4*(($end+1)*($end+2)/2)
            + 4*((101*102)/2)
            + 10;

diag "Checking small numbers";
foreach my $b (0 .. $end) {
  foreach my $e ($b .. $end) {
    my @p = grep { $_ >= $b && $_ <= $e } @primes;
    is_deeply( gen_primes($b,$e), \@p, "primes($b,$e)");
    is_deeply( gen_forprimes($b,$e), \@p, "forprimes {} $b,$e");
    is_deeply( gen_piterate($b,$e), \@p, "prime_iterator($b) while <= $e");
  }
}
SKIP: {
  skip "No OO iterator", (($end+1)*($end+2)/2) unless defined &Math::Prime::Util::prime_iterator_object;
  foreach my $b (0 .. $end) {
    foreach my $e ($b .. $end) {
      my @p = grep { $_ >= $b && $_ <= $e } @primes;
      is_deeply( gen_ooiterate($b,$e), \@p, "prime_iterator object $b to $e");
    }
  }
}

# TODO We should check boundaries around 1k*30, then segments around 256k*30 and 64k*30

my @lprimes = (~0 > 4294967295)
 ? (qw/18446744073709550671 18446744073709550681 18446744073709550717 18446744073709550719 18446744073709550771 18446744073709550773 18446744073709550791 18446744073709550873 18446744073709551113 18446744073709551163 18446744073709551191 18446744073709551253 18446744073709551263 18446744073709551293 18446744073709551337 18446744073709551359 18446744073709551427 18446744073709551437 18446744073709551521 18446744073709551533 18446744073709551557/)
 : (qw/4294966297 4294966337 4294966367 4294966373 4294966427 4294966441 4294966447 4294966477 4294966553 4294966583 4294966591 4294966619 4294966639 4294966651 4294966657 4294966661 4294966667 4294966769 4294966813 4294966829 4294966877 4294966909 4294966927 4294966943 4294966981 4294966997 4294967029 4294967087 4294967111 4294967143 4294967161 4294967189 4294967197 4294967231 4294967279 4294967291/);

diag "\nChecking numbers near end with iterator\n";
foreach my $bdelta (reverse 0 .. 100) {
  foreach my $edelta (reverse 0 .. $bdelta) {
    my ($b, $e) = (~0 - $bdelta, ~0 - $edelta);
    my @p = grep { $_ >= $b && $_ <= $e } @lprimes;
    is_deeply( gen_piterate($b,$e), \@p, "prime_iterator($b) while <= $e");
  }
}
SKIP: {
  skip "No OO iterator", ((101*102)/2) unless defined &Math::Prime::Util::prime_iterator_object;
  diag "\nChecking numbers near end with OO iterator\n";
  foreach my $bdelta (reverse 0 .. 100) {
    foreach my $edelta (reverse 0 .. $bdelta) {
      my ($b, $e) = (~0 - $bdelta, ~0 - $edelta);
      my @p = grep { $_ >= $b && $_ <= $e } @lprimes;
      is_deeply( gen_ooiterate($b,$e), \@p, "prime_iterator object $b to $e");
    }
  }
}

diag "\nChecking numbers near end with primes()\n";
foreach my $bdelta (reverse 0 .. 100) {
  foreach my $edelta (reverse 0 .. $bdelta) {
    my ($b, $e) = (~0 - $bdelta, ~0 - $edelta);
    my @p = grep { $_ >= $b && $_ <= $e } @lprimes;
    is_deeply( gen_primes($b,$e), \@p, "primes($b,$e)");
  }
}
diag "\nChecking numbers near end with forprimes.\n";
foreach my $bdelta (reverse 0 .. 100) {
  foreach my $edelta (reverse 0 .. $bdelta) {
    my ($b, $e) = (~0 - $bdelta, ~0 - $edelta);
    my @p = grep { $_ >= $b && $_ <= $e } @lprimes;
    is_deeply( gen_forprimes($b,$e), \@p, "forprimes {} $b,$e");
  }
}
diag "\nChecking numbers near end with segment primes().  Very slow.\n";
{
  my $b = $lprimes[-1] - 1;
  my $e = ~0;
  my @p = ($lprimes[-1]);
  diag "\n    Window around $lprimes[-1]\n";
  is_deeply( gen_primes({method => 'Segment'}, $b, $b), [], "primes($b,$b)");
  is_deeply( gen_primes({method => 'Segment'}, $b, $b+1), \@p, "primes($b,$b+1)");
  is_deeply( gen_primes({method => 'Segment'}, $b, $b+2), \@p, "primes($b,$b+2)");
  is_deeply( gen_primes({method => 'Segment'}, $b+1, $b+1), \@p, "primes($b+1,$b+1)");
  is_deeply( gen_primes({method => 'Segment'}, $b+1, $b+2), \@p, "primes($b+1,$b+2)");
  is_deeply( gen_primes({method => 'Segment'}, $b+2, $b+2), [], "primes($b+2,$b+2)");
  diag "\n    Window around $e\n";
  is_deeply( gen_primes({method => 'Segment'}, $e-2, $e-2), [], "primes($e-2,$e-2)");
  is_deeply( gen_primes({method => 'Segment'}, $e-2, $e), [], "primes($e-2,$e)");
  is_deeply( gen_primes({method => 'Segment'}, $e-1, $e), [], "primes($e-1,$e)");
  is_deeply( gen_primes({method => 'Segment'}, $e, $e), [], "primes($e,$e)");
}

#diag "\nChecking numbers near end with forprimes.  This will take a *very* long time.\n";
#foreach my $bdelta (reverse 0 .. 9) {
#  foreach my $edelta (reverse 0 .. $bdelta) {
#    my ($b, $e) = (~0 - $bdelta, ~0 - $edelta);
#    my @p = grep { $_ >= $b && $_ <= $e } @lprimes;
#    is_deeply( gen_forprimes($b,$e), \@p, "forprimes {} $b,$e");
#  }
#}

sub gen_primes {
  return primes(@_);
}
sub gen_forprimes {
  my($b, $e) = @_;
  my @p;
  forprimes { push @p, $_ } $b,$e;
  return \@p;
}
sub gen_piterate {
  my($b, $e) = @_;
  my @p;
  my $it = prime_iterator($b);
  my $n;
  while (1) {
    $n = $it->();
    last if $n > $e || $n == 0;
    push @p, $n;
  }
  return \@p;
}
sub gen_ooiterate {
  my($b, $e) = @_;
  my @p;
  my $it = prime_iterator_object($b);
  push @p, $it->iterate while $it->value <= $e;
  return \@p;
}
