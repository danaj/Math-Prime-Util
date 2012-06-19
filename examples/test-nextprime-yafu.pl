#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util qw/next_prime/;
use File::Temp qw/tempfile/;
use autodie;
my $maxdigits = (~0 <= 4294967295) ? 10 : 20;
$| = 1;  # fast pipes
my $num = shift || 10000;
my $yafu_fname = "yafu_batchfile_$$.txt";
$SIG{'INT'} = \&gotsig;

foreach my $digits (4 .. $maxdigits) {
  printf "%2d-digit numbers", $digits;
  my @narray = gendigits($digits, $num);
  print ".";
  my @mpuarray = mpu_next_primes(@narray);
  print ".";
  die "mpu_next_primes didn't get enough numbers" unless $#mpuarray == $#narray;
  my @yafuarray = yafu_next_primes(@narray);
  die "yafunext_primes didn't get enough numbers" unless $#yafuarray == $#narray;
  print ".";
  foreach my $n (@narray) {
    my $mpu = shift @mpuarray;
    my $yafu = shift @yafuarray;
    die "next_prime($n):  MPU: $mpu  YAFU: $yafu\n" unless $mpu == $yafu;
  }
  print ".";
  print "OK\n";
}

sub gendigits {
  my $digits = shift;
  die "Digits must be > 0" unless $digits > 0;
  my $howmany = shift;

  my $base = ($digits == 1) ? 0 : int(10 ** ($digits-1));
  my $max = int(10 ** $digits);
  $max = ~0 if $max > ~0;
  my @nums = map { $base+int(rand($max-$base)) } (1 .. $howmany);
  return @nums;
}

sub mpu_next_primes {
  my @nparray;
  push @nparray, next_prime($_) for @_;
  @nparray;
}

sub yafu_next_primes {
  my @nparray;
  # Yafu 1.31 seems to go out of its way to make it hard to process more than
  # one number at a time.  The batchfile system will infinite loop if the data
  # file isn't in the current directory.
  # It does its darndest to see if you're on a terminal or not, and if not it
  # just cuts you off after one number.  So any sort of tempfile or pipe stuff
  # just plain doesn't work.  Faking it using IO::*tty* would probably work.

  #my $fh = File::Temp->new;   # .... autodie
  #print $fh, "$_\n" for @_;
  #$fh->flush;

  # Shudder.  Read comments above about why I have to do this.
  open(my $fh, '>', $yafu_fname);
  print $fh "$_\n" for @_;
  close $fh;

  open my $yafu, "yafu \"nextprime(\@)\" -batchfile $yafu_fname |";
  while (<$yafu>) {
    next unless /ans = (\d+)/;
    push @nparray, $1;
  }
  close($yafu);
  @nparray;
}
sub gotsig { my $sig = shift; die "Die because SIG$sig\n"; }
END {
  unlink $yafu_fname if -e $yafu_fname;
  # YAFU leaves stuff around
  unlink "__tmpbatchfile" if -e "__tmpbatchfile";
  unlink "session.log" if -e "session.log";
}
