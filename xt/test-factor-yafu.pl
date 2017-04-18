#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util qw/factor urandomm/;
use File::Temp qw/tempfile/;
use Math::BigInt try => 'GMP,Pari';
use Config;
use autodie;
use Text::Diff;
my $maxdigits = 50;
$| = 1;  # fast pipes
my $num = 10000;
my $yafu_fname = "yafu_batchfile_$$.txt";
$SIG{'INT'} = \&gotsig;

{ # Test from 2 to 10000
  print "    2 -  1000"; test_array(    2 ..  1000);
  print " 1001 -  5000"; test_array( 1001 ..  5000);
  print " 5001 - 10000"; test_array( 5001 .. 10000);
}

foreach my $digits (5 .. $maxdigits) {
  printf "%5d %2d-digit numbers", $num, $digits;
  my @narray = gendigits($digits, $num);
  test_array(@narray);
  $num = int($num * 0.9) + 1; # reduce as we go
}

sub test_array {
  my @narray = @_;
  print ".";
  my @mpuarray = mpu_factors(@narray);
  print ".";
  my @yafuarray = yafu_factors(@narray);
  print ".";
  if ($#mpuarray != $#yafuarray) {
    die "MPU got $#mpuarray factors, YAFU got $#yafuarray\n";
  }
  foreach my $n (@narray) {
    my @mpu = @{shift @mpuarray};
    my @yafu = @{shift @yafuarray};
    die "mpu array is for the wrong n?" unless $n == shift @mpu;
    die "yafu array is for the wrong n?" unless $n == shift @yafu;
    my $diff = diff \@mpu, \@yafu, { STYLE => 'Table' };
    die "factor($n):\n$diff\n" if length($diff) > 0;
  }
  print ".";
  print "OK\n";
}

sub gendigits {
  my $digits = shift;
  die "Digits must be > 0" unless $digits > 0;
  my $howmany = shift;
  my ($base, $max);

  if ( 10**$digits < ~0) {
    $base = ($digits == 1) ? 0 : int(10 ** ($digits-1));
    $max = int(10 ** $digits);
    $max = ~0 if $max > ~0;
  } else {
    $base = Math::BigInt->new(10)->bpow($digits-1);
    $max = Math::BigInt->new(10)->bpow($digits) - 1;
  }
  my @nums = map { $base + urandomm($max-$base) } (1 .. $howmany);
  return @nums;
}

sub mpu_factors {
  my @piarray;
  push @piarray, [$_, factor($_)] for @_;
  @piarray;
}

sub yafu_factors {
  my @ns = @_;
  my @piarray;

  #my $fh = File::Temp->new;   # .... autodie
  #print $fh, "$_\n" for @_;
  #$fh->flush;

  # Shudder.  Yafu must have a file in the current directory.
  open(my $fh, '>', $yafu_fname);
  print $fh "$_\n" for @ns;
  close $fh;

  open my $yafu, "yafu \"factor(\@)\" -batchfile $yafu_fname |";
  my @curfactors;
  while (<$yafu>) {
    chomp;
    if (/^P(RP)?\d+ = (\d+)/) {
      push @curfactors, $2;
    } elsif (/^C\d+ = (\d+)/) {
      # Yafu didn't factor this one completely.  Sneakily do it ourselves.
      push @curfactors, factor( Math::BigInt->new("$1") );
    } elsif (/ans = (\d+)/ || /^1$/) {
      push @piarray, [shift @ns, sort {$a<=>$b} @curfactors];
      @curfactors = ();
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
  unlink "factor.log" if -e "factor.log";
}
