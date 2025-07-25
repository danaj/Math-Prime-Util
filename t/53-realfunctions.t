#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/ExponentialIntegral LogarithmicIntegral
                         RiemannR RiemannZeta LambertW/;

my $infinity = 20**20**20;

# We can check using Pari/GP:  -real(eint1(-x))
my %eivals = (
         -10 =>  -0.00000415696892968532438,
        -0.5 =>  -0.55977359477616,
        -0.1 =>  -1.8229239584193906660809,
      -0.001 =>  -6.33153936413615,
    -0.00001 => -10.9357198000436956,
 -0.00000001 => -17.843465089050832587,
 0.693147180559945 => 1.0451637801174927848446,           # log2
         1   =>  1.8951178163559367554665,
         1.5 =>  3.3012854491297978379574,
         2   =>  4.9542343560018901633795,
         5   =>  40.185275355803177455091,
         10  =>  2492.2289762418777591384,
         12  =>  14959.532666397528852292,
         20  =>  25615652.664056588820481,
         40  =>  6039718263611241.5783592,
         41  =>  16006649143245041.110700,
         79  =>  2.61362206325045575150640392249037e+32,
);

my %livals = (  # In pari these are:  -eint1(-log($n))
              0 =>  0,
           1.01 => -4.0229586739299358695031,
              2 =>  1.0451637801174927848446,
             10 =>  6.1655995047872979375230,
             24 =>  11.200315795232698830550,
           1000 =>  177.60965799015222668764,
         100000 =>  9629.8090010507982050343,
      100000000 =>  5762209.3754480314675691,
     4294967295 =>  203284081.95454158906409,
    10000000000 =>  455055614.58662307560953,
   100000000000 =>  4118066400.6216115150394,
);

# Values from T. R. Nicely for comparison
my %rvals = (
           1.01 =>  1.0060697180622924796117,
              2 =>  1.5410090161871318832885,
             10 =>  4.5645831410050902398658,
           1000 =>  168.35944628116734806491,
        1000000 =>  78527.399429127704858870,
       10000000 =>  664667.44756474776798535,
     4294967295 =>  203280697.51326064541983,
    10000000000 =>  455050683.30684692446315,
18446744073709551615 => 4.25656284014012122706963685602e17,
);

my %rzvals = (
            2   =>  0.6449340668482264364724151666,
            2.5 =>  0.3414872572509171797567696934,
            4.5 =>  0.0547075107614542640229672890,
            7   =>  0.0083492773819228268397975498,
            8.5 =>  0.0028592508824156277133439825,
           20.6 =>  0.0000006293391573578212882457,
);

my %lamvals = (
            -0.3678794411714423215955237701614608674458 => -0.99999995824889,  # Ideally this would be -1
            -.1 => -0.11183255915896296483356945682026584227264536229126586332968,
            0 => 0,
            0.3678794411714423215955237701614608674458 => 0.278464542761073795109358739022980155439470898229676526861772,
            1 => 0.567143290409783872999968662210355549753815787186512508135131,
            10 => 1.7455280027406993830743012648753899115,
            10000 => 7.2318460380933727064756185001412538839,
            100000000000 => 22.227122734961075624690200512898589272,
            18446744073709551615 => 40.656266572498926634921823566267328254,
);


plan tests => 3 + 6 + 1
              + scalar(keys(%eivals))
              + scalar(keys(%livals))
              + scalar(keys(%rvals))
              + scalar(keys(%rzvals))
              + scalar(keys(%lamvals))
              + 1
              ;

eval { LogarithmicIntegral(-1); };
like($@, qr/invalid/i, "li(-1) is invalid");
eval { RiemannR(0); };
like($@, qr/invalid/i, "R(0) is invalid");
eval { RiemannR(-1); };
like($@, qr/invalid/i, "R(-1) is invalid");

cmp_ok( ExponentialIntegral(0),         '<=',-$infinity, "Ei(0) is -infinity");
cmp_ok( ExponentialIntegral(-$infinity),'==', 0,         "Ei(-inf) is 0" );
cmp_ok( ExponentialIntegral($infinity), '>=', $infinity, "Ei(inf) is infinity");

cmp_ok( LogarithmicIntegral(0),         '==', 0,         "li(0) is 0");
cmp_ok( LogarithmicIntegral(1),         '<=',-$infinity, "li(1) is -infinity");
cmp_ok( LogarithmicIntegral($infinity), '>=', $infinity, "li(inf) is infinity");

# Example used in Math::Cephes
cmp_closeto( ExponentialIntegral(2.2), 5.732614700, 1e-06, "Ei(2.2)");

while (my($n, $ein) = each (%eivals)) {
  cmp_closeto( ExponentialIntegral($n), $ein, 0.00000001 * abs($ein), "Ei($n) ~= $ein");
}
while (my($n, $lin) = each (%livals)) {
  cmp_closeto( LogarithmicIntegral($n), $lin, 0.00000001 * abs($lin), "li($n) ~= $lin");
}

while (my($n, $rin) = each (%rvals)) {
  cmp_closeto( RiemannR($n), $rin, 0.00000001 * abs($rin), "R($n) ~= $rin");
}

while (my($n, $zin) = each (%rzvals)) {
  cmp_closeto( RiemannZeta($n), $zin, 0.00000001 * abs($zin), "Zeta($n) ~= $zin");
}
while (my($n, $lin) = each (%lamvals)) {
  # Machines with long double will be a little different near -1/e
  cmp_closeto( LambertW($n), $lin, 0.0000001 * abs($lin), "LambertW($n) ~= $lin");
}


# Put at end to not hit bug in pre-0.53 MPU::GMP
{
  my($n,$ein) = (170, 4.00120321792254728767739056606721e+71);
  cmp_closeto( ExponentialIntegral($n), $ein, 0.00000001 * abs($ein), "Ei($n) ~= $ein");
}



sub cmp_closeto {
  my $got = shift;
  my $expect = shift;
  my $tolerance = shift;
  my $message = shift;
  cmp_ok( abs($got - $expect), '<=', $tolerance, $message );
}
