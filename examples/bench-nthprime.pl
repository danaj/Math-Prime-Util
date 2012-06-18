#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util qw/nth_prime prime_precalc/;
use Benchmark qw/:all :hireswallclock/;
use Data::Dumper;

my $count = shift || -5;

#prime_precalc(1000000000);

srand(29);
my @darray;
push @darray, [gendigits($_,int(2700/($_*$_*$_)))]  for (2 .. 9);

my $sum;
foreach my $digits (3 .. 9) {
  my @digarray = @{$darray[$digits-2]};
  my $numitems = scalar @digarray;
  my $timing = cmpthese(
    $count,
    { "$digits" => sub { $sum += nth_prime($_) for @digarray }, },
    'none',
    );
  my $secondsper = $timing->[1]->[1];
  if ($timing->[0]->[1] eq 'Rate') {
    $secondsper =~ s/\/s$//;
    $secondsper = 1.0 / $secondsper;
  }
  $secondsper /= $numitems;
  my $timestr = (1.0 / $secondsper) . "/s per number";
  printf "%4d %2d-digit numbers: %s\n", $numitems, $digits, $timestr;
}

sub gendigits {
  my $digits = shift;
  die "Digits must be > 0" unless $digits > 0;
  my $num = shift;

  my $base = ($digits == 1) ? 0 : int(10 ** ($digits-1));
  my $max = int(10 ** $digits);
  $max = ~0 if $max > ~0;
  my @nums = map { $base+int(rand($max-$base)) } (1 .. $num);
  return @nums;
}
