#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/bernfrac bernreal harmfrac harmreal/;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};

my @A000367 = (qw/1 1 -1 1 -1 5 -691 7 -3617 43867 -174611 854513 -236364091 8553103 -23749461029 8615841276005 -7709321041217 2577687858367 -26315271553053477373 2929993913841559 -261082718496449122051 1520097643918070802691 -27833269579301024235023 596451111593912163277961 -5609403368997817686249127547 495057205241079648212477525 -801165718135489957347924991853 29149963634884862421418123812691 -2479392929313226753685415739663229 84483613348880041862046775994036021 -1215233140483755572040304994079820246041491/);
my @A002445 = (qw/1 6 30 42 30 66 2730 6 510 798 330 138 2730 6 870 14322 510 6 1919190 6 13530 1806 690 282 46410 66 1590 798 870 354 56786730 6 510 64722 30 4686 140100870 6 30 3318 230010 498 3404310 6 61410 272118 1410 6 4501770 6 33330 4326 1590 642 209191710 1518 1671270 42/);
my @A001008 = (qw/1 3 11 25 137 49 363 761 7129 7381 83711 86021 1145993 1171733 1195757 2436559 42142223 14274301 275295799 55835135 18858053 19093197 444316699 1347822955 34052522467 34395742267 312536252003 315404588903 9227046511387/);
my @A002805 = (qw/1 2 6 12 60 20 140 280 2520 2520 27720 27720 360360 360360 360360 720720 12252240 4084080 77597520 15519504 5173168 5173168 118982864 356948592 8923714800 8923714800 80313433200 80313433200 2329089562800/);

plan tests => 1 + 1 + 1;

subtest 'bernfrac (Bernoulli numbers)', sub {
  $#A000367 = 20 unless $extra;
  $#A002445 = 20 unless $extra;
  my @num = map { (bernfrac(2*$_))[0] }  0 .. $#A000367;
  my @den = map { (bernfrac(2*$_))[1] }  0 .. $#A002445;
  is_deeply( \@num, \@A000367, "B_2n numerators 0 .. $#A000367" );
  is_deeply( \@den, \@A002445, "B_2n denominators 0 .. $#A002445" );
  SKIP: {
    skip "bernfrac(502) only in EXTENDED_TESTING",1 unless $extra;
    my($num,$den) = bernfrac(502);
    my $sum = 0;
    $sum += $_ for split(//, $num);
    is($sum-$den, 242,"sumdigits(bernfrac(502) numerator) - denominator = 242");
  }
};

subtest 'harmfrac (Harmonic numbers)', sub {
  $#A001008 = 20 unless $extra;
  $#A002805 = 20 unless $extra;
  my @num = map { (harmfrac(1+$_))[0] }  0 .. $#A001008;
  my @den = map { (harmfrac(1+$_))[1] }  0 .. $#A002805;
  is_deeply( \@num, \@A001008, "H_n numerators 0 .. $#A001008" );
  is_deeply( \@den, \@A002805, "H_n denominators 0 .. $#A002805" );
};

subtest 'bernreal and harmreal', sub {
  my @bern_reals = (1,1/2,1/6,0,-1/30,0,1/42,0,-1/30,0,5/66,0,-691/2730,0,7/6,0,-3617/510,0,43867/798,0,-174611/330,0,854513/138,0,-236364091/2730);
  my $lbr = $#bern_reals;
  is_deeply([map {is_closeto(bernreal($_),$bern_reals[$_],1e-8)} 0..$lbr],
            [map { 1 } 0..$lbr],
            "bernreal(0..$lbr) within tolerance");
  my @harm_reals = (0/1,1/1,3/2,11/6,50/24,274/120,1764/720,13068/5040,109584/40320,1026576/362880,10628640/3628800,120543840/39916800,1486442880/479001600,19802759040/6227020800,283465647360/87178291200,4339163001600/1307674368000,70734282393600/20922789888000,1223405590579200/355687428096000,22376988058521600/6402373705728000,431565146817638400/121645100408832000,8752948036761600000/2432902008176640000);
  my $lhr = $#harm_reals;
  is_deeply([map {is_closeto(harmreal($_),$harm_reals[$_],1e-8)} 0..$lhr],
            [map { 1 } 0..$lhr],
            "harmreal(0..$lhr) within tolerance");
  SKIP: {
    skip "bernreal(46) and harmreal(46) with EXTENDED_TESTING",2 unless $extra;
    like( bernreal(46), qr/2115074863808199160560.145/, "bernreal(46)" );
    like( harmreal(46), qr/4.416687245986104750714329/, "harmreal(46)" );
  }
};

sub is_closeto {
  my($got,$exp,$tol) = @_;
  return 0 + (abs($got-$exp) <= $tol);
}
sub cmp_closeto {
  my($got,$exp,$tol,$mess) = @_;
  cmp_ok(abs($got - $exp), '<=', $tol, $mess);
}
