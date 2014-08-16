#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util qw/is_prime factor holf_factor/;

my $hrounds = 512*1024*1024;
for (2 .. 1e10) {
  my @fs;
  my $s_fact = join(".",sort {$a<=>$b} factor($_));

  my @p_holf;
  push @fs, $_;
  while (@fs) {
    my $n = pop @fs;
    if (is_prime($n)) { push @p_holf, $n; }
    else              { my @f = holf_factor($n,$hrounds);
                        die "Could not factor $n\n" if scalar @f == 1;
                        push @fs, @f;  }
  }
  my $s_holf = join(".",sort {$a<=>$b} @p_holf);

  die "$_ $s_fact  holf $s_holf\n" unless $s_fact eq $s_holf;

  print "$_\n" if ($_ % 100000) == 0;
}
