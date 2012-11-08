#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util qw/factor/;
use Math::Pari qw/factorint/;
use File::Temp qw/tempfile/;
use Math::BigInt try => 'GMP,Pari';
use Config;
use autodie;
use Text::Diff;
use Time::HiRes qw(gettimeofday tv_interval);
my $maxdigits = 50;
$| = 1;  # fast pipes
srand(87431);
my $num = 1000;

my $rgen = sub {
  my $range = shift;
  return 0 if $range <= 0;
  my $rbits = 0; { my $t = $range; while ($t) { $rbits++; $t >>= 1; } }
  while (1) {
    my $rbitsleft = $rbits;
    my $U = $range - $range;  # 0 or bigint 0
    while ($rbitsleft > 0) {
      my $usebits = ($rbitsleft > $Config{randbits}) ? $Config{randbits} : $rbitsleft;
      $U = ($U << $usebits) + int(rand(1 << $usebits));
      $rbitsleft -= $usebits;
    }
    return $U if $U <= $range;
  }
};

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
  my($start, $mpusec, $gnusec, $parisec, $diff);

  print ".";
  $start = [gettimeofday];
  my @mpuarray = mpu_factors(@narray);
  $mpusec = tv_interval($start);

  print ".";
  $start = [gettimeofday];
  my @gnuarray = gnu_factors(@narray);
  $gnusec = tv_interval($start);

  print ".";
  $start = [gettimeofday];
  my @pariarray = pari_factors(@narray);
  $parisec = tv_interval($start);

  print ".";
  die "MPU got ", scalar @mpuarray, " factors.  GNU factor got ",
      scalar @gnuarray, "\n" unless $#mpuarray == $#gnuarray;
  die "MPU got ", scalar @mpuarray, " factors.  Pari factor got ",
      scalar @pariarray, "\n" unless $#mpuarray == $#pariarray;
  foreach my $n (@narray) {
    my @mpu  = @{shift @mpuarray};
    my @gnu  = @{shift @gnuarray};
    my @pari = @{shift @pariarray};
    die "mpu array is for the wrong n?" unless $n == shift @mpu;
    die "gnu array is for the wrong n?" unless $n == shift @gnu;
    die "pari array is for the wrong n?" unless $n == shift @pari;
    $diff = diff \@mpu, \@gnu, { STYLE => 'Table' };
    die "factor($n): MPU/GNU\n$diff\n" if length($diff) > 0;
    my $diff = diff \@mpu, \@pari, { STYLE => 'Table' };
    die "factor($n): MPU/Pari\n$diff\n" if length($diff) > 0;
  }
  print ".";
  # We should ignore the small digits, since we're comparing direct
  # Perl functions with multiple command line invocations.  It really
  # doesn't make sense until we're over 1ms per number.
  printf "OK  MPU:%7.2f ms  GNU:%7.2f ms  Pari:%7.2f ms\n",
         (($mpusec*1000) / scalar @narray),
         (($gnusec*1000) / scalar @narray),
         (($parisec*1000) / scalar @narray);
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
  my @nums = map { $base + $rgen->($max-$base) } (1 .. $howmany);
  return @nums;
}

sub mpu_factors {
  my @piarray;
  push @piarray, [$_, factor($_)] for @_;
  @piarray;
}

sub gnu_factors {
  my @ns = @_;
  my @piarray;
  my $numpercommand = int( (4000-30)/(length($ns[-1])+1) );

  while (@ns) {
    my $cs = '/usr/bin/factor';
    foreach my $n (1 .. $numpercommand) {
      last unless @ns;
      $cs .= " " . shift @ns;
    }
    my $fout = qx{$cs};
    my @flines = split(/\n/, $fout);
    foreach my $fline (@flines) {
      $fline =~ s/^(\d+): //;
      push @piarray, [$1, split(/ /, $fline)];
    }
  }
  @piarray;
}

sub pari_factors {
  my @piarray;
  foreach my $n (@_) {
    my @factors;
    my ($pn,$pc) = @{factorint($n)};
    # Map the Math::Pari objects returned into Math::BigInts, because Pari will
    # throw a hissy fit later when we try to compare them to anything else.
    push @piarray, [ $n, map { (Math::BigInt->new($pn->[$_])) x $pc->[$_] } (0 .. $#$pn) ];
  }
  @piarray;
}
