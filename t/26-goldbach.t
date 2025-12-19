#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/minimal_goldbach_pair
                         goldbach_pair_count
                         goldbach_pairs/;
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $usexs = Math::Prime::Util::prime_get_config->{'xs'};

my @gb = (
  [0, []],
  [1, []],
  [2, []],
  [3, []],
  [4, [2]],
  [5, [2]],
  [6, [3]],
  [7, [2]],
  [8, [3]],
  [9, [2]],
  [10, [3,5]],
  [11, []],
  [12, [5]],
  [13, [2]],
  [14, [3,7]],
  [15, [2]],
  [16, [3,5]],
  [17, []],
  [18, [5,7]],
  [19, [2]],
  [20, [3,7]],
  [21, [2]],
  [22, [3,5,11]],
  [23, []],
  [24, [5,7,11]],
  [25, [2]],
  [26, [3,7,13]],
  [27, []],
  [28, [5,11]],
  [29, []],
  [30, [7,11,13]],
  [31, [2]],
  [32, [3,13]],
  [50, [3,7,13,19]],
  [66, [5,7,13,19,23,29]],
  [130, [3,17,23,29,41,47,59]],
  [258, [7,17,19,29,31,47,59,61,67,79,101,107,109,127]],
  [4022, [3,19,79,103,199,229,283,313,331,349,379,409,439,463,523,631,661,691,709,751,769,853,859,1021,1051,1069,1171,1231,1291,1303,1429,1471,1483,1549,1741,1753,1783,1801,1861,1879,1933,1993,2011]],
);

my %bigmin = (
  60119912 => 1093,
  "15317795894" => 2017,
  "3325581707333960528" => 9781,
  "83778272185315920949659591651127238812" => 5441,  # Randomly chosen
);

plan tests => 1+scalar(keys %bigmin)      # minimal
            + 1                           # count
            + 1;                          # pairs


###### minimal

{
  my(@exp, @got);
  for my $g (@gb) {
    my($n,$L) = @$g;
    push @exp, (@$L > 0) ? $L->[0] : undef;
    push @got, minimal_goldbach_pair($n);
  }
  is_deeply(\@got, \@exp, "minimal_goldbach_pair for small inputs");
}
while (my($n, $p) = each (%bigmin)) {
  SKIP: {
    skip "skipping minimal_goldbach_pair($n) without EXTENDED_TESTING",1 unless $n < ~0 || $extra;
    is(minimal_goldbach_pair($n),$p,"minimal_goldbach_pair($n) = $p");
  }
}

###### count
{
  my(@exp, @got);
  for my $g (@gb) {
    my($n,$L) = @$g;
    push @exp, scalar @$L;
    push @got, goldbach_pair_count($n);
  }
  is_deeply(\@got, \@exp, "goldbach_pair_count for small inputs");
}

###### pairs
{
  my(@exp, @got);
  for my $g (@gb) {
    my($n,$L) = @$g;
    push @exp, $L;
    push @got, [goldbach_pairs($n)];
  }
  is_deeply(\@got, \@exp, "goldbach_pairs for small inputs");
}
