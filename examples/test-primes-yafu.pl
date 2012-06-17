#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util qw/primes/;
use File::Temp qw/tempfile/;
use autodie;
use Text::Diff;
my $maxdigits = (~0 <= 4294967295) ? 10 : 18;
$| = 1;  # fast pipes
my $num = 5000;
my $interval = 8000;
my $yafu_fname = "yafu_batchfile_$$.txt";
$SIG{'INT'} = \&gotsig;


# Note -- yafu 1.31 will not sieve 19 digit numbers.  E.g.:
#   primes(8631424695497106432,8631424695497114432,0)
# gives:
#   input too high


foreach my $digits (3 .. $maxdigits) {
  printf "%5d %2d-digit numbers", $num, $digits;
  my @narray = gendigits($digits, $num);
  print ".";
  my @mpuarray = mpu_primes(@narray);
  print ".";
  #die "mpu_next_primes didn't get enough numbers" unless $#mpuarray-1 == $#narray;
  my @yafuarray = yafu_primes(@narray);
  #die "yafunext_primes didn't get enough numbers" unless $#yafuarray-1 == $#narray;
  print ".";
  if ($#mpuarray != $#yafuarray) {
    die "MPU got $#mpuarray primes, YAFU got $#yafuarray\n";
  }
  foreach my $n (@narray) {
    my @mpu = @{shift @mpuarray};
    my @yafu = @{shift @yafuarray};
    die "mpu array is for the wrong n?" unless $n == shift @mpu;
    die "yafu array is for the wrong n?" unless $n == shift @yafu;
    my $diff = diff \@mpu, \@yafu, { STYLE => 'Table' };
    die "primes($n,$n+$interval):\n$diff\n" if length($diff) > 0;
  }
  print ".";
  print "OK\n";
  $num = int($num * 0.75); # reduce as we go
}

sub gendigits {
  my $digits = shift;
  die "Digits must be > 0" unless $digits > 0;
  my $howmany = shift;

  my $base = ($digits == 1) ? 0 : int(10 ** ($digits-1));
  my $max = int(10 ** $digits);
  #$max = ~0 if $max > ~0;
  $max = ~0-$interval if $max > (~0-$interval); # special for us
  my @nums = map { $base+int(rand($max-$base)) } (1 .. $howmany);
  return @nums;
}

sub mpu_primes {
  my @piarray;
  push @piarray, [$_, @{primes($_, $_+$interval)}] for @_;
  @piarray;
}

sub yafu_primes {
  my @ns = @_;
  my @piarray;
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
  print $fh "$_,", $_+$interval, ",0\n" for @ns;
  close $fh;

  open my $yafu, "yafu \"primes(\@)\" -pscreen -batchfile $yafu_fname |";
  my @curprimes;
  while (<$yafu>) {
    chomp;
    if (/^\d+/) {
      push @curprimes, split(/ /);
    } elsif (/ans = (\d+)/) {
      foreach my $p (@curprimes) { die "Entry is '$p'" unless $p =~ /^\d+$/; }
      push @piarray, [shift @ns, @curprimes];
      @curprimes = ();
    }
  }
  close($yafu);
  @piarray;
}
sub gotsig { my $sig = shift; die "Die because SIG$sig\n"; }
END {
  unlink $yafu_fname if -e $yafu_fname;
  # YAFU leaves stuff around
  unlink "__tmpbatchfile" if -e "__tmpbatchfile";
  unlink "session.log" if -e "session.log";
}
