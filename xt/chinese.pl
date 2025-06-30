#!/usr/bin/env perl
use strict;
use warnings;
$| = 1;  # fast pipes

use Math::Prime::Util qw/chinese lcm chinese2 powint urandomb urandomm/;
use Math::ModInt qw(mod);
use Math::ModInt::ChineseRemainder qw(cr_combine);

my $limit = shift || 10_000;
my $printmod = int(($limit+77) / 78);


cmp_chinese($limit, powint(2,16), 2);
cmp_chinese($limit, powint(2,32), 2);
cmp_chinese($limit, powint(2,16), 4);
cmp_chinese($limit, powint(2,40), 3);
cmp_chinese($limit, powint(2,31), 13);
print "\nDone\n";


sub rpairs {
  my($lim, $num) = @_;
  my @p;
  for (1..$num) {
    my $mod = 1+urandomm($lim);  # random modulo between 1 and $lim inclusive
    push @p, [ urandomm($mod), $mod ]
  }
  @p;
}

sub printpairs {
  my(@rm) = @_;
  print "\ninput: ";
  print join(" ", map { "[$_->[0] $_->[1]]" } @rm), "\n";
  1;
}

sub cmp_chinese {
  my($ntests, $lim, $num) = @_;

  my $size = $lim > 2**30 ? "large" : "small";
  print "Running $limit random tests with $num $size inputs...\n";
  for my $n (1 .. $ntests) {
    print '.' unless $n % $printmod;

    my @rm = rpairs($lim, $num);
    #printpairs(@rm);

    my $mic = cr_combine( map { mod($_->[0],$_->[1]) } @rm );
    if ($mic->is_undefined) {
      my $mpu_res = chinese(@rm);
      printpairs(@rm) && die "MIC: undef  MPU: $mpu_res\n" if defined $mpu_res;
      next;
    }
    my $mic_res = $mic->residue;
    my $mic_mod = $mic->modulus;

    my($mpu_res,$mpu_mod) = chinese2(@rm);
    printpairs(@rm) && die "MIC: $mic_res $mic_mod  MPU: undef\n" if !defined $mpu_res;
    printpairs(@rm) && die "MIC: $mic_res $mic_mod   MPU: $mpu_res  $mpu_mod" if $mpu_res != $mic_res || $mpu_mod != $mic_mod;
  }
  print "\n";
}
