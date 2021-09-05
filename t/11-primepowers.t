#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/prime_powers is_prime_power
                         next_prime_power prev_prime_power
                         prime_power_count prime_power_count_approx
                         prime_power_count_lower prime_power_count_upper
                         nth_prime_power nth_prime_power_approx
                         nth_prime_power_lower nth_prime_power_upper
                         /;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $usexs = Math::Prime::Util::prime_get_config->{'xs'};

my @small = (2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25, 27, 29, 31, 32, 37, 41, 43, 47, 49, 53, 59, 61, 64, 67, 71, 73, 79, 81, 83, 89, 97);

my @A025528 = (0, 1, 2, 3, 4, 4, 5, 6, 7, 7, 8, 8, 9, 9, 9, 10, 11, 11, 12, 12, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 18, 18, 18, 18, 18, 19, 19, 19, 19, 20, 20, 21, 21, 21, 21, 22, 22, 23, 23, 23, 23, 24, 24, 24, 24, 24, 24, 25, 25, 26, 26, 26, 27, 27, 27, 28, 28, 28, 28, 29, 29, 30, 30);
my @A267712 = (7, 35, 193, 1280, 9700, 78734, 665134, 5762859, 50851223, 455062595, 4118082969, 37607992088, 346065767406, 3204942420923, 29844572385358, 279238346816392, 2623557174778438, 24739954338671299, 234057667428388198, 2220819603016308079);
$#A025528 = 40;
$#A267712 = ($usexs || $extra) ? 7 : 6;

my %samples = (72=>263, 89=>353, 311=>1831, 423=>2677, 814=>5857,
               1509=>12149, 4484=>42089, 9163=>93893, 10957=>114671,
               20942=>234599, 51526=>629819, 417867=>6071249,
               717421=>10843321);
$samples{1031932}   =   16002401 if $usexs || $extra;
$samples{2687492}   =   44442791 if $usexs || $extra;
$samples{8337143}   =  147948547 if $usexs || $extra;
$samples{21208173}  =  397416121 if $usexs || $extra;
$samples{47490104}  =  930475697 if $usexs || $extra;


plan tests =>   0   # is_prime_power (not tested here)
              + 3   # prime_powers
              + 3   # next_prime_power
              + 4   # prev_prime_power
              + 6   # prime_power_count (first set)
              + 2   # prime_power_count (second set)
              + 4   # prime_power_count bounds
              + 3   # nth_prime_power
              + $extra * scalar(keys %samples)
              + 9   # nth_prime_power bounds
              + 0;


###### is_prime_power

# See t/26-ispower.t

###### prime_powers

my $pp1k = prime_powers(1000);

is_deeply( [grep { $_ <= 300 } @$pp1k],
           [grep { is_prime_power($_) } 1..300],
           "prime_powers(300)" );

is_deeply( [map { prime_powers($_) } 1..100],
           [map { my $n=$_; [grep { $_ <= $n } @small] } 1..100],
           "prime_powers(1..100)" );

my $base1 = 1441897;
is_deeply( prime_powers($base1,$base1+1000),
           [map { $base1+$_ } (34,36,52,66,84,106,112,120,156,160,172,174,190,246,262,276,294,312,330,354,370,382,402,420,424,430,436,444,454,480,496,504,514,532,540,556,562,612,616,630,634,652,682,694,702,714,724,730,736,744,756,760,772,820,826,834,846,886,900,930,952,966,972,976,990,1000)],
           "prime_powers($base1, $base1 + 1000)" );

###### next_prime_power

is_deeply( [map { next_prime_power($_-1) } @small],
           \@small,
           "next_prime_power");

is_deeply( [map { next_prime_power(2**$_-1) } 1..9],
           [2, 4, 8, 16, 32, 64, 128, 256, 512],
           "next_prime_power(2^i+1)" );
is_deeply( [map { next_prime_power(2**$_) } 0..9],
           [2, 3, 5, 9, 17, 37, 67, 131, 257, 521],
           "next_prime_power(2^i)" );

###### prev_prime_power

is_deeply( [map { prev_prime_power($_) } 0..2], [undef,undef,undef], "prev_prime_power(0..2) = undef" );

is_deeply( [map { prev_prime_power($_+1) } @small],
           \@small,
           "prev_prime_power");

is_deeply( [map { prev_prime_power(2**$_+1) } 1..9],
           [2, 4, 8, 16, 32, 64, 128, 256, 512],
           "prev_prime_power(2^i+1)" );
is_deeply( [map { prev_prime_power(2**$_) } 1..9],
           [undef, 3, 7, 13, 31, 61, 127, 251, 509],
           "prev_prime_power(2^i)" );


###### prime_power_count

is(prime_power_count(0), 0, "prime_power_count(0) = 0");
is(prime_power_count(1), 0, "prime_power_count(1) = 0");
is_deeply( [map { prime_power_count(1+$_) } 0..$#A025528], \@A025528,  "prime_power_count(n) for 1..".scalar(@A025528) );
is_deeply( [map { prime_power_count(10**(1+$_)) } 0..$#A267712], \@A267712,  "prime_power_count(10^n) for 1..".scalar(@A267712) );

is(prime_power_count(12345678), 809830, "prime_power_count(12345678) = 809830");
is(prime_power_count(123456, 133332), 847, "prime_power_count(123456,133332) = 847");

{
  my $tlimit = 400;
  my @exact = (0);
  for my $l (@$pp1k) {
    last if $l > $tlimit;
    push @exact, $exact[-1] while scalar(@exact) < $l;
    push @exact, $exact[-1] + 1;
  }
  is_deeply( [map { prime_power_count($_) } 0 .. $#exact],
             \@exact,
             "prime_power_count(0..$#exact)" );
}
{
  my(@expect,@got);
  for my $lo (0 .. 30) {
    for my $hi (0 .. 30) {
      push @expect, [$lo,$hi,scalar(grep { $_ >= $lo && $_ <= $hi } @small)];
      push @got,    [$lo,$hi,prime_power_count($lo,$hi)];
    }
  }
  is_deeply( \@got, \@expect, "prime_power_count ranges 0 .. 40" );
}

###### prime_power_count_{upper,lower,approx}

is(check_count_bounds(513, 117), 1, "prime_power count bounds for 513");
is(check_count_bounds(5964377, 411055), 1, "prime_power count bounds for 5964377");

is_deeply( [map { check_count_bounds($_, prime_power_count($_)) } 0..100],
           [map { 1 } 0..100],
           "prime_power count bounds for small numbers" );

is_deeply( [map { check_count_bounds($samples{$_},$_) } keys %samples],
           [map { 1 } keys %samples],
           "prime_power count bounds for small samples" );

###### nth_prime_power

is( nth_prime_power(0), undef, "nth_prime_power(0) returns undef" );
is_deeply( [map { nth_prime_power($_) } 1 .. 100],
           [@$pp1k[0..100-1]],
           "first 100 prime powers with nth_prime_power" );
is( nth_prime_power(1 << 12), 37993, "37993 is the 2^12th prime power" );
if ($extra) {
  while (my($n, $npp) = each (%samples)) {
    is( nth_prime_power($n), $npp, "nth_prime_power($n) = $npp" );
  }
}

###### nth_prime_power{upper,lower,approx}

is( nth_prime_power_lower(0), undef, "nth_prime_power_lower(0) returns undef" );
is( nth_prime_power_upper(0), undef, "nth_prime_power_upper(0) returns undef" );
is( nth_prime_power_approx(0), undef, "nth_prime_power_approx(0) returns undef" );

is(check_nth_bounds(86, 343), 1, "nth_prime_power(86) bounds");
is(check_nth_bounds(123456, 1628909), 1, "nth_prime_power(123456) bounds");
is(check_nth_bounds(5286238, 91241503), 1, "nth_prime_power(5286238) bounds");
is(check_nth_bounds(46697909, 914119573), 1, "nth_prime_power(46697909) bounds");

is_deeply( [map { check_nth_bounds($_, $pp1k->[$_-1]) } 1..70],
           [map { 1 } 1..70],
           "nth_prime_power(1..70) bounds" );

is_deeply( [map { check_nth_bounds($_, $samples{$_}) } keys %samples],
           [map { 1 } keys %samples],
           "nth_prime_power bounds for small samples" );



sub check_count_bounds {
  my($n, $count) = @_;

  return "bad lower count for $n" if prime_power_count_lower($n) > $count;
  return "bad upper count for $n" if prime_power_count_upper($n) < $count;
  my $approx = prime_power_count_approx($n);
  return "approx count too low for $n"  if $approx < 0.9 * $count;
  return "approx count too high for $n" if $approx > 1.1 * $count;
  1;
}

sub check_nth_bounds {
  my($n, $nth) = @_;

  return "bad nth lower for $n" if nth_prime_power_lower($n) > $nth;
  return "bad nth upper for $n" if nth_prime_power_upper($n) < $nth;
  my $approx = nth_prime_power_approx($n);
  return "approx nth too low for $n"  if $approx < 0.9 * $nth;
  return "approx nth too high for $n" if $approx > 1.1 * $nth;
  1;
}
