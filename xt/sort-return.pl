#!/usr/bin/env perl
use strict;
use warnings;

use Math::Prime::Util qw/vecsort/;
use Math::Prime::Util::PP;


testit("sort direct", \&sort_direct);
testit("sort indirect", \&sort_indirect);
testit("sort workaround", \&sort_workaround);
print "\n";
testit("vecsort direct", \&vecsort_direct);
testit("vecsort indirect", \&vecsort_indirect);
testit("vecsort workaround", \&vecsort_workaround);
print "\n";
testit("PP vecsort direct", \&ppvecsort_direct);
testit("PP vecsort indirect", \&ppvecsort_indirect);
testit("PP vecsort workaround", \&ppvecsort_workaround);

sub testit {
  my($name, $func) = @_;
  my @v = (12,13,14,11);

  my @X = $func->(@v);
  my $x = $func->(@v);
  $x = "<undef>" unless defined $x;

  printf "%21s  %7s  %7s\n", $name, scalar(@X), $x;
}

sub sort_direct {
  return sort { $a<=> $b } @_;
}
sub sort_indirect {
  my @p = sort { $a<=> $b } @_;
  return @p;
}
sub sort_workaround {
  return scalar @_ unless wantarray;
  return sort { $a<=> $b } @_;
}

sub vecsort_direct {
  return vecsort(@_);
}
sub vecsort_indirect {
  my @p = vecsort(@_);
  return @p;
}
sub vecsort_workaround {
  return scalar @_ unless wantarray;
  return vecsort(@_);
}

sub ppvecsort_direct {
  return Math::Prime::Util::PP::vecsort(@_);
}
sub ppvecsort_indirect {
  my @p = Math::Prime::Util::PP::vecsort(@_);
  return @p;
}
sub ppvecsort_workaround {
  return scalar @_ unless wantarray;
  return Math::Prime::Util::PP::vecsort(@_);
}
