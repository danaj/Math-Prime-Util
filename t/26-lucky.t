#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/lucky_numbers nth_lucky is_lucky/;
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};

my @lucky = (qw/1 3 7 9 13 15 21 25 31 33 37 43 49 51 63 67 69 73 75 79 87 93 99 105 111 115 127 129 133 135 141 151 159 163 169 171 189 193 195/);
# 10 randomly selected samples up to 1e6
my %samples = (1144=>10239, 3378=>34639, 4996=>53595, 24406=>308119, 26427=>336201, 43449=>578395, 67655=>935911, 69526=>964101, 69985=>971011, 70970=>985983);
my $nsamples = scalar(keys %samples);

plan tests => 3 + 1*$extra                # lucky_numbers
            + 5 + (1+$nsamples)*$extra    # is_lucky
            + 3 + (0+$nsamples)*$extra;   # nth_lucky

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

###############################################################################

is_deeply( [grep { is_lucky($_) } 0..200], \@lucky, "is_lucky for 0 to 200" );

ok( is_lucky(42975), "42975 is a lucky number" );
ok( !is_lucky(513), "513 is not a lucky number" );
ok( !is_lucky(49023), "49023 is not a lucky number" );
ok( is_lucky(120001), "120001 is a lucky number" );

if ($extra) {
  ok( is_lucky(9999997), "9999997 is a lucky number" );
  while (my($n, $ln) = each (%samples)) {
    ok( is_lucky($ln), "is_lucky($ln) returns true" );
  }
}

###############################################################################

is( nth_lucky(0), undef, "nth_lucky(0) returns undef" );
is_deeply( [map { nth_lucky($_) } 1 .. 153], lucky_numbers(1000), "lucky numbers under 1000 with nth_lucky" );
is( nth_lucky(1 << 12), 42975, "42975 is the 2^12th lucky number" );
if ($extra) {
  while (my($n, $ln) = each (%samples)) {
    is( nth_lucky($n), $ln, "nth_lucky($n) = $ln" );
  }
}
