#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/lucky_numbers is_lucky
                         lucky_count lucky_count_approx
                         lucky_count_lower lucky_count_upper
                         nth_lucky nth_lucky_approx
                         nth_lucky_lower nth_lucky_upper/;
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $usexs = Math::Prime::Util::prime_get_config->{'xs'};

my @lucky = (qw/1 3 7 9 13 15 21 25 31 33 37 43 49 51 63 67 69 73 75 79 87 93 99 105 111 115 127 129 133 135 141 151 159 163 169 171 189 193 195/);
# 10 randomly selected samples up to 1e6
my %samples = (1144=>10239, 3378=>34639, 4996=>53595, 24406=>308119, 26427=>336201, 43449=>578395, 67655=>935911, 69526=>964101, 69985=>971011, 70970=>985983);
my $nsamples = scalar(keys %samples);

plan tests => 5 + 1*$extra                # lucky_numbers
            + 1                           # lucky_numbers ranges
            + 7 + (1+$nsamples)*$extra    # is_lucky
            + 3 + (0+$nsamples)*$extra    # nth_lucky
            + 1                           # lucky_count
            + 1                           # lucky_count ranges
            + 4                           # lucky_count bounds
            + 10;                         # nth_lucky bounds


###### lucky_numbers sieve

# Simple test for correct set
is_deeply( lucky_numbers(200), \@lucky, "lucky_numbers(200)" );

my $ln1k = lucky_numbers(1000);  # We will use this later

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
SKIP: {
  skip "skip mid-size lucky sieve for PP without EXTENDED_TESTING",2
    unless $usexs || $extra;
  my $l = lucky_numbers(350000);
  is( scalar(@$l), 27420, "Lucky numbers under 350k: 27420" );
  my $sum = 0;
  $sum += $_ for @$l;
  is( $sum, 4574808744, "Lucky numbers under 350k: correct sum" );
}

# ranges
SKIP: {
  skip "lucky_numbers should take ranges.  Skipping until it does.",1;
  my(@expect,@got);
  for my $lo (0 .. 40) {
    for my $hi (0 .. 40) {
      #push @expect, [$lo,$hi,grep { $_ >= $lo && $_ <= $hi } @lucky];
      #push @got,    [$lo,$hi,@{lucky_numbers($lo,$hi)}];
    }
  }
  is_deeply( \@got, \@expect, "all lucky numbers ranges 0 .. 40" );
}

###### is_lucky

is_deeply( [grep { is_lucky($_) } 0..200], \@lucky, "is_lucky for 0 to 200" );

ok(  is_lucky(42975), "42975 is a lucky number" );
ok( !is_lucky(513), "513 is not a lucky number" );
ok( !is_lucky(49023), "49023 is not a lucky number" );
ok(  is_lucky(120001), "120001 is a lucky number" );
ok( !is_lucky(1000047), "1000047 is not a lucky number" );
ok(  is_lucky(1000071), "1000047 is a lucky number" );

if ($extra) {
  ok( is_lucky(9999997), "9999997 is a lucky number" );
  while (my($n, $ln) = each (%samples)) {
    ok( is_lucky($ln), "is_lucky($ln) returns true" );
  }
}

###### lucky_count
{
  my @exact = (0);
  for my $l (@$ln1k) {
    push @exact, $exact[-1] while scalar(@exact) < $l;
    push @exact, $exact[-1] + 1;
  }
  is_deeply( [map { lucky_count($_) } 0 .. $#exact],
             \@exact,
             "lucky_count(0..$#exact)" );
}

{
  my(@expect,@got);
  for my $lo (0 .. 40) {
    for my $hi (0 .. 40) {
      push @expect, [$lo,$hi,scalar(grep { $_ >= $lo && $_ <= $hi } @lucky)];
      push @got,    [$lo,$hi,lucky_count($lo,$hi)];
    }
  }
  is_deeply( \@got, \@expect, "lucky_count ranges 0 .. 40" );
}

###### lucky_count_{upper,lower,approx}

is(check_count_bounds(513, 86), 1, "lucky count bounds for 513");
is(check_count_bounds(5964377, 376288), 1, "lucky count bounds for 5964377");

is_deeply( [map { check_count_bounds($_, lucky_count($_)) } 0..100],
           [map { 1 } 0..100],
           "lucky count bounds for small numbers" );

is_deeply( [map { check_count_bounds($samples{$_},$_) } keys %samples],
           [map { 1 } keys %samples],
           "lucky count bounds for small samples" );

###### nth_lucky

is( nth_lucky(0), undef, "nth_lucky(0) returns undef" );
is_deeply( [map { nth_lucky($_) } 1 .. 153], $ln1k, "lucky numbers under 1000 with nth_lucky" );
is( nth_lucky(1 << 12), 42975, "42975 is the 2^12th lucky number" );
if ($extra) {
  while (my($n, $ln) = each (%samples)) {
    is( nth_lucky($n), $ln, "nth_lucky($n) = $ln" );
  }
}

###### nth_lucky{upper,lower,approx}

is( nth_lucky_lower(0), undef, "nth_lucky_lower(0) returns undef" );
is( nth_lucky_upper(0), undef, "nth_lucky_upper(0) returns undef" );
is( nth_lucky_approx(0), undef, "nth_lucky_approx(0) returns undef" );

is(check_nth_bounds(86, 511), 1, "nth_lucky(86) bounds");
is(check_nth_bounds(123456, 1795453), 1, "nth_lucky(123456) bounds");
is(check_nth_bounds(5286238, 99999979), 1, "nth_lucky(5286238) bounds");
is(check_nth_bounds(46697909, 999999991), 1, "nth_lucky(46697909) bounds");
is(check_nth_bounds(1<<31, 55291335127), 1, "nth_lucky(2^31) bounds");

is_deeply( [map { check_nth_bounds($_, $ln1k->[$_-1]) } 1..100],
           [map { 1 } 1..100],
           "nth_lucky(1..100) bounds" );

is_deeply( [map { check_nth_bounds($_, $samples{$_}) } keys %samples],
           [map { 1 } keys %samples],
           "nth_lucky bounds for small samples" );



sub check_count_bounds {
  my($n, $count) = @_;

  return "bad lower count for $n" if lucky_count_lower($n) > $count;
  return "bad upper count for $n" if lucky_count_upper($n) < $count;
  my $approx = lucky_count_approx($n);
  return "approx count too low for $n"  if $approx < 0.9 * $count;
  return "approx count too high for $n" if $approx > 1.1 * $count;
  1;
}

sub check_nth_bounds {
  my($n, $nth) = @_;

  return "bad nth lower for $n" if nth_lucky_lower($n) > $nth;
  return "bad nth upper for $n" if nth_lucky_upper($n) < $nth;
  my $approx = nth_lucky_approx($n);
  return "approx nth too low for $n"  if $approx < 0.9 * $nth;
  return "approx nth too high for $n" if $approx > 1.1 * $nth;
  1;
}
