#!/usr/bin/env perl

# Verify prime gaps, version 1.0
# Dana Jacobsen, 2014
#
# This is an alternative to T.R. Nicely's cglp4 program from:
#      https://faculty.lynchburg.edu/~nicely/#Downloads
# This runs 2-4x faster on my machines.  If cglp4 can use PFGW, then it will
# cross over speed around 3000 digits, and PFGW is much faster at 10k+.
#
# It will use the extra-strong BPSW test plus a Frobenius-Underwood test
# for the endpoints so is more stringent about endpoint testing (cglp4 uses
# the strong BPSW test).
#
# The gaps are in one of the formats:
#    <gapsize>  <P1-expression>
#    <gapsize>  <merit>  <P1-expression>
#    <merit>  <gapsize>  PRP#### = <P1-expression>
#
# This program will DIE if an invalid gap is found.  I believe this is
# preferable to printing a 0 result in a list which may be thousands of
# lines long, and hence missed.  If the gaps have been properly supplied,
# this should never come up.

use warnings;
use strict;
use Math::BigInt lib=>"GMP";
use Math::Prime::Util qw/:all/;
use Math::Prime::Util::GMP;     # Ensure we're using this
use Time::HiRes qw(gettimeofday tv_interval);
$|=1;

# TODO: Use a command line argument
my $use_pfgw = 0;
#my $pfgw_exec = "/users/jacobsen/src/pfgw-3.7.10/pfgw64";
my $pfgw_exec = "pfgw64";
my $pfgw_thresh = 2400;  # PFGW faster only for this many digits

my $fstart = [gettimeofday];
my $procn = 0;
while (<>) {
  chomp;
  next if /^#/ || /^\s*$/;
  my($mer, $gap, $expr);
  if (/^\s*(\d+)  (\S+)  (\S+)$/) {
    ($mer, $gap, $expr) = ($2, $1, $3);
  } elsif (/^\s*(\S+)\s+(\d+)\s+PRP\d+ = (.*)/) {
    ($mer, $gap, $expr) = ($1, $2, $3);
  } elsif (/^(\d+)  (\S+)$/) {
    ($gap, $expr) = ($1, $2);
  } else {
    warn "skipping $_\n";
    next;
  }
  $procn++;
  my $start = [gettimeofday];
  $expr =~ s/^1\*//;

  my $orig_expr = $expr;
  my $n = numerate($expr);
  my $end = $n + $gap;
  my $dstr = length($n) . "D";
  my $dstr2 = length($end) . "D";
  my $log2n = int(length($n) * 3.322);  # approx

  printf "G=%7d %10.2fs Checking P1 ($dstr)...\r", $gap, tv_interval($start);
  die "beg of '$expr' is not prime" unless test($n);
  printf "G=%7d %10.2fs Checking P2 ($dstr2)...  \r", $gap, tv_interval($start);
  die "end of '$expr' is not prime" unless test($end);
  my $next;

  # To avoid all the overhead of timing and printing, for very small
  # gaps we can just call next_prime which will check all the interior
  # points.  The only downside is that we're losing some manual control.
  if (0 && $gap < 15000 && $log2n < 800) {
    printf "G=%7d %10.2fs Checking P1 ($dstr) interval...   \r", $gap, tv_interval($start);
    $next = next_prime($n);
  } else {
    my $depth = int( 1.2 * $log2n * $log2n * log($log2n) );
    printf "G=%7d %10.2fs Sieving to $depth ...%s \r", $gap, tv_interval($start), " " x 30;
    my @list = sieve_range($n+1, $gap-1, $depth);
    my $gapstart = [gettimeofday];
    my $ntests = scalar(@list);
    my $i = 0;
    my $nexti = 1;
    printf "G=%7d %10.2fs Checking P1 ($dstr) + %d...     \r", $gap, tv_interval($start), $list[0]-$n;
    foreach my $rgap (@list) {
      my $pgap = $rgap + 1;  # We sieved from $n+1
      die "Interior point $expr + $pgap is prime\n" if testint($n+$pgap);
      $i++;
      if ($i >= $nexti) {
        my $startint = tv_interval($start);
        my $gaptime = tv_interval($gapstart);
        my $est = $startint + ($ntests-$i) * $gaptime/$i;
        printf "G=%7d %10.2fs (est %.2fs) Checking P1 ($dstr) + $pgap...   \r",  $gap, $startint, $est;
        my $display_intervals = int(0.4 / ($gaptime/$i));
        #$display_intervals = 256 if $display_intervals > 256;
        $nexti = $i + $display_intervals;
      }
    }
    $next = $end;
  }
  if ($next == $end) {
    printf "G=%7d P1=%-40sOK BPSW+FU=1 (%.3fs)\n",
      $gap, $expr, tv_interval($start);
  } else {
    die "gap $gap for $expr should be ", $next-$n, "\n";
  }
}
printf "\n Errors=0.  OK=%d.  T=%.3f.\n", $procn, tv_interval($fstart);

sub numerate {
  my $expr = shift;
  $expr =~ s/\b(\d+)#/primorial($1)/g;
  $expr =~ s/\^/**/;
  $expr =~ s/(\d+)/ Math::BigInt->new("$1") /g;
  my $n = eval $expr;
  die "Cannot eval: $expr\n" if !defined $n;
  return $n;
}

sub test {
  my $n = shift;
  return is_bpsw_prime($n) && is_frobenius_underwood_pseudoprime($n);
}

sub testint {
  my $n = shift;
  if ($use_pfgw && length($n) >= $pfgw_thresh) {
    return 0 if system("$pfgw_exec -k -Cquiet -f0 -u0 -q\"$n\" >/dev/null 2>1");
  }
  return is_bpsw_prime($n) && is_frobenius_underwood_pseudoprime($n);
}
