#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/lucky_numbers/;
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};

my @lucky = (qw/1 3 7 9 13 15 21 25 31 33 37 43 49 51 63 67 69 73 75 79 87 93 99 105 111 115 127 129 133 135 141 151 159 163 169 171 189 193 195/);

plan tests => 3 + 1*$extra;

# Simple test for correct set
is_deeply( lucky_numbers(200), \@lucky, "lucky_numbers(200)" );

{ # Verify we get the right set back for each n from 0 to 200
  my(@expect,@triangle);
  for my $n (0 .. 200) {
    push @triangle, @{lucky_numbers($n)}, ":";
    push @expect, (grep { $_ <= $n } @lucky), ":";
  }
  is_deeply( \@triangle, \@expect, "lucky numbers for each n from 0 to 200");
}

is( scalar(@{lucky_numbers(12345)}), 1349, "correct count for lucky_numbers(12345)" );

if ($extra) {
  is( scalar(@{lucky_numbers(145845)}), 12345, "correct count for lucky_numbers(145845)" );
  #is( scalar(@{lucky_numbers(1795453)}), 123456, "correct count for lucky_numbers(1795453)" );
}
