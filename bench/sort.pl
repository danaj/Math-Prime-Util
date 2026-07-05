#!perl
use warnings;
use strict;
use feature 'say';
use Time::HiRes qw( gettimeofday tv_interval );

use Math::Prime::Util ":all";
use Sort::XS;
use Sort::Key::Radix;
use Sort::Key;
use Sort::DJB;
use List::MoreUtils;

my $narrays = 10000;

for my $len (10,100,1000,10000,100000,1000000) {
  $narrays /= 10 if $len == 1000;
  $narrays /= 10 if $len == 100000;
  $narrays /= 10 if $len == 1000000;

  my(@times) = (0) x 8;
  for (1..20) {
    my @ints;
    for (0..$narrays-1) {
      #$ints[$_] = [map { irand } 1..$len];          # 32-bit random
      $ints[$_] = [map { irand64 } 1..$len];         # 64-bit random
      #$ints[$_] = [map { 1..3 } 1..divint($len,3)]; # sawtooth
    }
    $times[0] += time_sort(\@ints);
    $times[1] += time_sortxsq(\@ints);
    $times[2] += time_sortkey(\@ints);
    $times[3] += time_sortkeyradix(\@ints);
    $times[4] += time_vecsort(\@ints);
    $times[5] += time_lmu(\@ints);
    $times[6] += time_djb(\@ints);
    $times[7] += time_vecsorti(\@ints);
  }
  show_res($times[0], $times[5], "LMU::qsort", 20*$narrays, $len);
  show_res($times[0], $times[0], "sort", 20*$narrays, $len);
  show_res($times[0], $times[2], "Sort::Key::usort", 20*$narrays, $len);
  show_res($times[0], $times[1], "Sort::XS::quick_sort", 20*$narrays, $len);
  show_res($times[0], $times[3], "Sort::Key::Radix::usort", 20*$narrays, $len);
  show_res($times[0], $times[6], "Sort::DJB", 20*$narrays, $len);
  show_res($times[0], $times[4], "vecsort", 20*$narrays, $len);
  show_res($times[0], $times[7], "vecsorti", 20*$narrays, $len);
  print "\n";
}

sub show_res {
  my($tsort, $tthis, $name, $narr, $len) = @_;
  $tthis = $tsort if $tsort == 0;
  my $us = $tthis * 1e6 / $narr;
  my $mult = sprintf "%4.1fx %8.1fuS", $tsort/$tthis, $us;

  say "$mult  $name  $len random 64-bit integers";
  return $tthis;
}

sub time_sort {
  my $ints = shift;
  my $t0 = [gettimeofday];
  for my $t (@$ints) { my @sorted = sort {$a<=>$b} @$t; }
  return tv_interval($t0);
}

sub time_vecsort {
  my $ints = shift;
  my $t0 = [gettimeofday];
  for my $t (@$ints) { my @sorted = vecsort($t); }
  return tv_interval($t0);
}

sub time_sortxsq {
  my $ints = shift;
  my $t0 = [gettimeofday];
  for my $t (@$ints) { my @sorted = Sort::XS::quick_sort($t); }
  return tv_interval($t0);
}
sub time_sortxs {
  my $ints = shift;
  my $t0 = [gettimeofday];
  for my $t (@$ints) { my @sorted = Sort::XS::xsort(list => $t, algorithm => 'quick', type => 'integer'); }
  return tv_interval($t0);
}
sub time_sortkey {
  my $ints = shift;
  my $t0 = [gettimeofday];
  for my $t (@$ints) { my @sorted = Sort::Key::usort(@$t); }
  return tv_interval($t0);
}
sub time_sortkeyradix {
  my $ints = shift;
  my $t0 = [gettimeofday];
  for my $t (@$ints) { my @sorted = Sort::Key::Radix::usort(@$t); }
  return tv_interval($t0);
}
sub time_lmu {
  my $ints = shift;
  my @nints = map { [@$_] } @$ints;
  my $t0 = [gettimeofday];
  for my $t (@nints) { List::MoreUtils::qsort {$a<=>$b} @$t; }
  return tv_interval($t0);
}
sub time_djb {
  my $ints = shift;
  my @nints = map { [@$_] } @$ints;
  my $t0 = [gettimeofday];
  for my $t (@nints) { Sort::DJB::sort_uint64($t); }
  return tv_interval($t0);
}
sub time_vecsorti {
  my $ints = shift;
  my @nints = map { [@$_] } @$ints;
  my $t0 = [gettimeofday];
  for my $t (@nints) { vecsorti($t); }
  return tv_interval($t0);
}
