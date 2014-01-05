#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec::Functions;
use FindBin;
use Time::HiRes qw(gettimeofday tv_interval);
use bigint;
use Math::NumSeq;
$|++; #flush the output buffer after every write() or print() function
my $use64;
BEGIN { no bigint; $use64 = (~0 > 4294967295); }

compare('Primes',
        10000000,
        "$FindBin::Bin/../bin/primes.pl 1 LASTNUM",
        q/perl -MMath::NumSeq::Primes -e 'my $seq = Math::NumSeq::Primes->new; while (1) { my $v = ($seq->next)[1]; last if $v > LASTNUM; print "$v\n"; }'/);

compare('Twin',
        10000000,
        "$FindBin::Bin/../bin/primes.pl --twin 1 LASTNUM",
        q/perl -MMath::NumSeq::TwinPrimes -e 'my $seq = Math::NumSeq::TwinPrimes->new; while (1) { my $v = ($seq->next)[1]; last if $v > LASTNUM; print "$v\n"; }'/);

compare('Sophie Germain',
        10000000,
        "$FindBin::Bin/../bin/primes.pl --sophie 1 LASTNUM",
        q/perl -MMath::NumSeq::SophieGermainPrimes -e 'my $seq = Math::NumSeq::SophieGermainPrimes->new; while (1) { my $v = ($seq->next)[1]; last if $v > LASTNUM; print "$v\n"; }'/);

# Why Math::Prime::Util::is_prime instead of Math::Prime::XS::is_prime?
#   1) it's much faster for the palindrome tests
#   2) it supports bignums, which is required for Fib, Euclid, Lucas, etc.

compare('Palindromic',
        $use64   ?  '10**11'  :  '10**10',
        "$FindBin::Bin/../bin/primes.pl --palin 1 LASTNUM",
        q/perl -MMath::Prime::Util=is_prime -MMath::NumSeq::Palindromes -e 'my $seq = Math::NumSeq::Palindromes->new; while (1) { my $v = ($seq->next)[1]; last if $v > LASTNUM; print "$v\n" if is_prime($v); }'/);

# Sadly Math::NumSeq::LucasNumbers uses OEIS 204 (1,3) instead of OEIS 32 (-1,2)
# and neither package offers a way to adjust.
#compare('Lucas',
#        '10**100',
#        "$FindBin::Bin/../bin/primes.pl --lucas 1 LASTNUM",
#        q/perl -MMath::Prime::Util=is_prime -MMath::NumSeq::LucasNumbers -e 'my $seq = Math::NumSeq::LucasNumbers->new; while (1) { my $v = ($seq->next)[1]; last if $v > LASTNUM; print "$v\n" if is_prime($v); }'/);

compare('Fibonacci',
        '10**100',
        "$FindBin::Bin/../bin/primes.pl --fib 1 LASTNUM",
        q/perl -MMath::Prime::Util=is_prime -MMath::NumSeq::Fibonacci -e 'my $seq = Math::NumSeq::Fibonacci->new; while (1) { my $v = ($seq->next)[1]; last if $v > LASTNUM; print "$v\n" if is_prime($v); }'/);

compare('Euclid',
        '10**200',
        "$FindBin::Bin/../bin/primes.pl --euclid 1 LASTNUM",
        q/perl -MMath::Prime::Util=is_prime -MMath::NumSeq::Primorials -e 'my $seq = Math::NumSeq::Primorials->new; while (1) { my $v = ($seq->next)[1] + 1; last if $v > LASTNUM; print "$v\n" if is_prime($v); }'/);

compare('Lucky',
        '100000',
        "$FindBin::Bin/../bin/primes.pl --lucky 1 LASTNUM",
        q/perl -MMath::Prime::Util=is_prime -MMath::NumSeq::LuckyNumbers -e 'my $seq = Math::NumSeq::LuckyNumbers->new; while (1) { my $v = ($seq->next)[1]; last if $v > LASTNUM; print "$v\n" if is_prime($v); }'/);


sub compare {
  my($name, $end, $command_scr, $command_mns) = @_;
  no bigint;
  $command_scr =~ s/LASTNUM/$end/;
  $command_mns =~ s/LASTNUM/$end/;

  printf "%15s to %8s", $name, $end;

  my $start_scr = [gettimeofday];
  my @scr = split /\s+/, qx/$command_scr/;
  my $seconds_scr = tv_interval($start_scr);

  printf " (%7d).  primes.pl %6.2fs", scalar @scr, $seconds_scr;

  my $start_mns = [gettimeofday];
  my @mns = split /\s+/, qx/$command_mns/;
  my $seconds_mns = tv_interval($start_mns);

  printf "  Math::NumSeq %6.2fs\n", $seconds_mns;

  die "$name:  primes.pl generated ", scalar @scr, " results.  MNS generated ", scalar @mns, " results." if scalar @scr != scalar @mns;

  foreach my $i (0 .. $#scr) {
    die "$name prime $i not equal:\n  primes.pl: $scr[$i]\n  MNumSeq:  $mns[$i]\n"
        if $scr[$i] != $mns[$i];
  }
}
