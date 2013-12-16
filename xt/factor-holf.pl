#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util qw/is_prime/;
use List::Util qw/min max/;
my $count = shift || -2;
my $is64bit = (~0 > 4294967295);
my $maxdigits = ($is64bit) ? 20 : 10;  # Noting the range is limited for max.

my $hrounds = 64*1024*1024;
for (2 .. 100000000) {
  my @fs;
  my $s_fact = join(".",sort {$a<=>$b} Math::Prime::Util::factor($_));

  my @p_holf;
  push @fs, $_;
  while (@fs) {
    my $n = pop @fs;
    if (is_prime($n)) { push @p_holf, $n; }
    else              { push @fs, Math::Prime::Util::holf_factor($n); }
  }
  my $s_holf = join(".",sort {$a<=>$b} @p_holf);

  die "$_ $s_fact  holf $s_holf\n" unless $s_fact eq $s_holf;

  print "$_\n" if ($_ % 100000) == 0;
}
