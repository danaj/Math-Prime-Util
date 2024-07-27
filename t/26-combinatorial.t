#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/binomial factorial subfactorial factorialmod fubini
                         falling_factorial rising_factorial
                         forcomb forperm forderange formultiperm
                         numtoperm permtonum randperm shuffle
                         csrand
                        /;
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $usexs  = Math::Prime::Util::prime_get_config->{'xs'};
my $usegmp = Math::Prime::Util::prime_get_config->{'gmp'};

use Math::BigInt try => "GMP,Pari";

my %perms = (
 0 => [[]],
 1 => [[0]],
 2 => [[0,1],[1,0]],
 3 => [[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,0,1],[2,1,0]],
 4 => [[0,1,2,3],[0,1,3,2],[0,2,1,3],[0,2,3,1],[0,3,1,2],[0,3,2,1],[1,0,2,3],[1,0,3,2],[1,2,0,3],[1,2,3,0],[1,3,0,2],[1,3,2,0],[2,0,1,3],[2,0,3,1],[2,1,0,3],[2,1,3,0],[2,3,0,1],[2,3,1,0],[3,0,1,2],[3,0,2,1],[3,1,0,2],[3,1,2,0],[3,2,0,1],[3,2,1,0]],
);

my @binomials = (
 [ 35,16, 4059928950 ],             # We can do this natively even in 32-bit
 [ 40,19, "131282408400" ],         # We can do this in 64-bit
 [ 67,31, "11923179284862717872" ], # ...and this
 [ 228,12, "30689926618143230620" ],# But the result of this is too big.
 [ 177,78, "3314450882216440395106465322941753788648564665022000" ],
 [ 61,17, "536830054536825" ],
 [ -11,22, 64512240 ],
 [ -12,23, -286097760 ],
 [ -23,-26, -2300 ],     # Kronenburg extension
 [ -12,-23, -705432 ],   # same
 [  12,-23, 0 ],
 [  12,-12, 0 ],
 [ -12,0, 1 ],
);

# TODO: Add a bunch of combs here:  "5,3" => [[..],[..],[..]],

plan tests => 2                        # factorial
            + 1 + 1*$extra + 2         # factorialmod
            + 3                        # subfactorial
            + 1                        # fubini
            + 6                        # binomial
            + 3                        # binomial (Loeb)
            + scalar(@binomials)       # binomial
            + 2                        # falling_factorial
            + 2                        # rising_factorial
            + 7 + 4                    # Combinations
            + scalar(keys(%perms)) + 1 # Permutations
            + 4                        # Multiset Permutations
            + 5                        # Derangements
            + 5 + 5 + 1                # numtoperm, permtonum
            + 5                        # randperm
            + 5                        # shuffle
            + 11                       # randperm using csrand
            + 3                        # shuffle using csrand
            ;

###### factorial
sub fact { my $n = Math::BigInt->new("$_[0]"); $n->bfac; }
{
  my @result = map { factorial($_) } 0 .. 100;
  my @expect = map { fact($_) } 0 .. 100;
  is_deeply( \@result, \@expect, "Factorials 0 to 100" );
}
ok(!eval { factorial(-5); }, "factorial(-5) gives error");

###### factorialmod
{
  my @result = map { my $m=$_; map { factorialmod($_,$m) } 0..$m-1; } 1 .. 40;
  my @expect = map { my $m=$_; map { factorial($_) % $m; } 0..$m-1; } 1 .. 40;
  is_deeply( \@result, \@expect, "factorialmod n! mod m for m 1 to 50, n 0 to m" );
}
if ($extra) {
  # Very slow in PP without GMP
  is( factorialmod(5000001,"8000036000054000027"), "4179720539133404343", "factorialmod with large n and large composite non-square-free m" );
}
SKIP: {
  skip "medium size factorialmods in PP",2 unless $usexs || $extra;
  is( factorialmod(1000000000,1000000008), 0, "1000000000! mod 1000000008 is zero" );
  is( factorialmod(50000,10000019), 8482159, "50000! mod 10000019" );
}

###### subfactorial
is_deeply( [ map { subfactorial($_) } 0..23 ],
           [qw/1 0 1 2 9 44 265 1854 14833 133496 1334961 14684570 176214841 2290792932 32071101049 481066515734 7697064251745 130850092279664 2355301661033953 44750731559645106 895014631192902121 18795307255050944540 413496759611120779881 9510425471055777937262/],
           "subfactoral(n) for 0..23" );
is(subfactorial(110), "5842828822584214646127804296800556812003401310647230252823417994828330749128488139372248218138294820842482275906806317309680576864190217329860297083368061950972635498019334565561", "subfactorial(110)");
ok(!eval { subfactorial(-5); }, "subfactorial(-5) gives error");

###### fubini
is_deeply( [ map { fubini($_) } 0..23 ],
           [qw/1 1 3 13 75 541 4683 47293 545835 7087261 102247563 1622632573 28091567595 526858348381 10641342970443 230283190977853 5315654681981355 130370767029135901 3385534663256845323 92801587319328411133 2677687796244384203115 81124824998504073881821 2574844419803190384544203 85438451336745709294580413/],
           "fubini(n) for 0..23" );

###### binomial
{
  my @bin; for my $n (0..10) { for my $k (0..10) { push @bin, [$n,$k]; } }
  is_deeply( [map { binomial($_->[0],$_->[1]) } @bin],
             [qw/1 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 1 2 1 0 0 0 0 0 0 0 0 1 3 3 1 0 0 0 0 0 0 0 1 4 6 4 1 0 0 0 0 0 0 1 5 10 10 5 1 0 0 0 0 0 1 6 15 20 15 6 1 0 0 0 0 1 7 21 35 35 21 7 1 0 0 0 1 8 28 56 70 56 28 8 1 0 0 1 9 36 84 126 126 84 36 9 1 0 1 10 45 120 210 252 210 120 45 10 1/],
             "binomial(0..10,0..10)" );
}
{
  my @bin; for my $n (-10..-1) { for my $k (0..10) { push @bin, [$n,$k]; } }
  is_deeply( [map { binomial($_->[0],$_->[1]) } @bin],
             [qw/1 -10 55 -220 715 -2002 5005 -11440 24310 -48620 92378 1 -9 45 -165 495 -1287 3003 -6435 12870 -24310 43758 1 -8 36 -120 330 -792 1716 -3432 6435 -11440 19448 1 -7 28 -84 210 -462 924 -1716 3003 -5005 8008 1 -6 21 -56 126 -252 462 -792 1287 -2002 3003 1 -5 15 -35 70 -126 210 -330 495 -715 1001 1 -4 10 -20 35 -56 84 -120 165 -220 286 1 -3 6 -10 15 -21 28 -36 45 -55 66 1 -2 3 -4 5 -6 7 -8 9 -10 11 1 -1 1 -1 1 -1 1 -1 1 -1 1/],
             "binomial(-10..-1,0..10)" );
}
{
  my @bin; for my $n (0..10) { for my $k (-10..-1) { push @bin, [$n,$k]; } }
  is_deeply( [map { binomial($_->[0],$_->[1]) } @bin],
             [qw/0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0/],
             "binomial(0..10,-10..-1)" );
}
{
  my @bin; for my $n (-10..-1) { for my $k (-10..-1) { push @bin, [$n,$k]; } }
  is_deeply( [map { binomial($_->[0],$_->[1]) } @bin],
             [qw/1 0 0 0 0 0 0 0 0 0 -9 1 0 0 0 0 0 0 0 0 36 -8 1 0 0 0 0 0 0 0 -84 28 -7 1 0 0 0 0 0 0 126 -56 21 -6 1 0 0 0 0 0 -126 70 -35 15 -5 1 0 0 0 0 84 -56 35 -20 10 -4 1 0 0 0 -36 28 -21 15 -10 6 -3 1 0 0 9 -8 7 -6 5 -4 3 -2 1 0 -1 1 -1 1 -1 1 -1 1 -1 1/],
             "binomial(-10..-1,-10..-1)" );
}

is_deeply( [map { binomial(13, $_) } -15 .. 15],
           [qw/0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 13 78 286 715 1287 1716 1716 1287 715 286 78 13 1 0 0/],
           "binomial(13,n) for n in -15 .. 15" );
is_deeply( [map { binomial(-13, $_) } -15 .. 15],
           [qw/91 -13 1 0 0 0 0 0 0 0 0 0 0 0 0 1 -13 91 -455 1820 -6188 18564 -50388 125970 -293930 646646 -1352078 2704156 -5200300 9657700 -17383860/],
           "binomial(-13,n) for n in -15 .. 15" );

# Binomials from Loeb 1995
{ my $n;
  is_deeply( [map { $n=$_; map { binomial($n,$_) } 0..$n } 0..7],
             [qw/1 1 1 1 2 1 1 3 3 1 1 4 6 4 1 1 5 10 10 5 1 1 6 15 20 15 6 1 1 7 21 35 35 21 7 1/],
             "binomial: Loeb 1995 Region 1 (positive n, positive k)" );
  is_deeply( [map { $n=$_; map { binomial($n,$_) } 0..6 } -5 .. -1],
             [qw/1 -5 15 -35 70 -126 210 1 -4 10 -20 35 -56 84 1 -3 6 -10 15 -21 28 1 -2 3 -4 5 -6 7 1 -1 1 -1 1 -1 1/],
             "binomial: Loeb 1995 Region 2 (negative n, positive k)" );
  is_deeply( [map { $n=$_; map { binomial($n,$_) } -6..$n } -6 .. -1],
             [qw/1 -5 1 10 -4 1 -10 6 -3 1 5 -4 3 -2 1 -1 1 -1 1 -1 1/],
             "binomial: Loeb 1995 Region 3 (negative n, negative k)" );
  # Region 4 is positive n and k > n.  We are always 0.
  # Region 5 is positive n and negative k.  We are always 0.
  # Region 6 is negative n and n < k < 0.  We are always 0.
}

# Selected binomials
foreach my $r (@binomials) {
  my($n, $k, $exp) = @$r;
  is( binomial($n,$k), $exp, "binomial($n,$k)) = $exp" );
}

###### falling_factorial
{
  my @s;
  for my $k (0..10) { push @s, map { falling_factorial($_,$k) } -10..10; }
  is_deeply( \@s,
             [qw/1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10 110 90 72 56 42 30 20 12 6 2 0 0 2 6 12 20 30 42 56 72 90 -1320 -990 -720 -504 -336 -210 -120 -60 -24 -6 0 0 0 6 24 60 120 210 336 504 720 17160 11880 7920 5040 3024 1680 840 360 120 24 0 0 0 0 24 120 360 840 1680 3024 5040 -240240 -154440 -95040 -55440 -30240 -15120 -6720 -2520 -720 -120 0 0 0 0 0 120 720 2520 6720 15120 30240 3603600 2162160 1235520 665280 332640 151200 60480 20160 5040 720 0 0 0 0 0 0 720 5040 20160 60480 151200 -57657600 -32432400 -17297280 -8648640 -3991680 -1663200 -604800 -181440 -40320 -5040 0 0 0 0 0 0 0 5040 40320 181440 604800 980179200 518918400 259459200 121080960 51891840 19958400 6652800 1814400 362880 40320 0 0 0 0 0 0 0 0 40320 362880 1814400 -17643225600 -8821612800 -4151347200 -1816214400 -726485760 -259459200 -79833600 -19958400 -3628800 -362880 0 0 0 0 0 0 0 0 0 362880 3628800 335221286400 158789030400 70572902400 29059430400 10897286400 3632428800 1037836800 239500800 39916800 3628800 0 0 0 0 0 0 0 0 0 0 3628800/],
  "falling_factorial(-10..10, 0..10)" );
}
is_deeply( [map { falling_factorial($_->[0],$_->[1]) } ([515,7],[516,7],[568,7],[89,10],[103,101])],
           [qw/9222879462222182400 9349716704335257600 18378924259448108160 18452514066426316800 49514503582430902037733576272908866745450829110572462415026402773499383329208111416070720536941769246326758192988646046611441067207574945792000000000000000000000000/],
           "falling_factorial selected values");
###### rising_factorial
{
  my @s;
  for my $k (0..10) { push @s, map { rising_factorial($_,$k) } -10..10; }
  is_deeply( \@s,
             [qw/1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10 90 72 56 42 30 20 12 6 2 0 0 2 6 12 20 30 42 56 72 90 110 -720 -504 -336 -210 -120 -60 -24 -6 0 0 0 6 24 60 120 210 336 504 720 990 1320 5040 3024 1680 840 360 120 24 0 0 0 0 24 120 360 840 1680 3024 5040 7920 11880 17160 -30240 -15120 -6720 -2520 -720 -120 0 0 0 0 0 120 720 2520 6720 15120 30240 55440 95040 154440 240240 151200 60480 20160 5040 720 0 0 0 0 0 0 720 5040 20160 60480 151200 332640 665280 1235520 2162160 3603600 -604800 -181440 -40320 -5040 0 0 0 0 0 0 0 5040 40320 181440 604800 1663200 3991680 8648640 17297280 32432400 57657600 1814400 362880 40320 0 0 0 0 0 0 0 0 40320 362880 1814400 6652800 19958400 51891840 121080960 259459200 518918400 980179200 -3628800 -362880 0 0 0 0 0 0 0 0 0 362880 3628800 19958400 79833600 259459200 726485760 1816214400 4151347200 8821612800 17643225600 3628800 0 0 0 0 0 0 0 0 0 0 3628800 39916800 239500800 1037836800 3632428800 10897286400 29059430400 70572902400 158789030400 335221286400/],
  "rising_factorial(-10..10, 0..10)" );
}
is_deeply( [map { rising_factorial($_->[0],$_->[1]) } ([509,7],[510,7],[562,7],[80,10],[103,101])],
           [qw/9222879462222182400 9349716704335257600 18378924259448108160 18452514066426316800 6760937240727169751346751449031021029092236987417146776093364751481076175432048515956305908925637116481562056123160956910787676051553407749205364947724300581490631820332063331242347041889126973440000000000000000000000000/],
           "rising_factorial selected values");


###### forcomb
{ my @p = (); forcomb { push @p, [@_] } 0;
  is_deeply( [@p], [[]], "forcomb 0" ); }
{ my @p = (); forcomb { push @p, [@_] } 1;
  is_deeply( [@p], [[],[0]], "forcomb 1" ); }
{ my @p = (); forcomb { push @p, [@_] } 0,0;
  is_deeply( [@p], [[]], "forcomb 0,0" ); }
{ my @p = (); forcomb { push @p, [@_] } 5,0;
  is_deeply( [@p], [[]], "forcomb 5,0" ); }
{ my @p = (); forcomb { push @p, [@_] } 5,6;
  is_deeply( [@p], [], "forcomb 5,6" ); }
{ my @p = (); forcomb { push @p, [@_] } 5,5;
  is_deeply( [@p], [[0,1,2,3,4]], "forcomb 5,5" ); }
{ my @p = (); forcomb { push @p, [@_] } 3;
  is_deeply( [@p], [[],[0],[1],[2],[0,1],[0,2],[1,2],[0,1,2]], "forcomb 3 (power set)" ); }

{ my @data = (qw/apple bread curry/);
  my @p = (); forcomb { push @p, [@data[@_]] } @data,2;
  my @e = ([qw/apple bread/],[qw/apple curry/],[qw/bread curry/]);
  is_deeply( \@p,\@e, "forcomb 3,2" ); }
{ my @data = (qw/ant bee cat dog/);
  my @p = (); forcomb { push @p, [@data[@_]] } @data,3;
  my @e = ([qw/ant bee cat/],[qw/ant bee dog/],[qw/ant cat dog/],[qw/bee cat dog/]);
  is_deeply( \@p,\@e, "forcomb 4,3" ); }

{ my $b = binomial(20,15);
  my $s = 0; forcomb { $s++ } 20,15;
  is($b, 15504, "binomial(20,15) is 15504");
  is($s, $b, "forcomb 20,15 yields binomial(20,15) combinations"); }


###### forperm
while (my($n, $expect) = each (%perms)) {
  my @p = (); forperm { push @p, [@_] } $n;
  is_deeply( \@p, $expect, "forperm $n" );
}

{ my $s = 0; forperm { $s++ } 7;
  is($s, factorial(7), "forperm 7 yields factorial(7) permutations"); }

###### formultiperm
{ my @p; formultiperm { push @p, [@_] } [];
  is_deeply(\@p, [], "formultiperm []"); }

{ my @p; formultiperm { push @p, [@_] } [1,2,2];
  is_deeply(\@p, [ [1,2,2], [2,1,2], [2,2,1] ], "formultiperm 1,2,2"); }

{ my @p; formultiperm { push @p, [@_] } [qw/a a b b/];
  is_deeply(\@p, [map{[split(//,$_)]} qw/aabb abab abba baab baba bbaa/], "formultiperm a,a,b,b"); }

{ my @p; formultiperm { push @p, join("",@_) } [qw/a a b b/];
  is_deeply(\@p, [qw/aabb abab abba baab baba bbaa/], "formultiperm aabb"); }

###### forderange
{ my @p; forderange { push @p, [@_] } 0;
  is_deeply(\@p, [[]], "forderange 0"); }
{ my @p; forderange { push @p, [@_] } 1;
  is_deeply(\@p, [], "forderange 1"); }
{ my @p; forderange { push @p, [@_] } 2;
  is_deeply(\@p, [[1,0]], "forderange 2"); }
{ my @p; forderange { push @p, [@_] } 3;
  is_deeply(\@p, [[1,2,0],[2,0,1]], "forderange 3"); }
{ my $n=0; forderange { $n++ } 7; is($n, 1854, "forderange 7 count"); }

###### numtoperm / permtonum

is_deeply([numtoperm(0,0)],[],"numtoperm(0,0)");
is_deeply([numtoperm(1,0)],[0],"numtoperm(1,0)");
is_deeply([numtoperm(1,1)],[0],"numtoperm(1,1)");
is_deeply([numtoperm(5,15)],[0,3,2,4,1],"numtoperm(5,15)");
is_deeply([numtoperm(24,987654321)],[0,1,2,3,4,5,6,7,8,9,10,13,11,21,14,20,17,15,12,22,18,19,23,16],"numtoperm(24,987654321)");

is(permtonum([]),0,"permtonum([])");
is(permtonum([0]),0,"permtonum([0])");
is(permtonum([6,3,4,2,5,0,1]),4768,"permtonum([6,3,4,2,5,0,1])");
is(permtonum([reverse(0..14),15..19]),"1790774578500738480","permtonum( 20 )");
is(permtonum([reverse(0..12),reverse(13..25)]),"193228515634198442606207999","permtonum( 26 )");

is(permtonum([numtoperm(14,8467582)]),8467582,"permtonum(numtoperm)");

###### randperm
# TODO: better randperm tests

is(@{[randperm(0)]},0,"randperm(0)");
is(@{[randperm(1)]},1,"randperm(1)");
is(@{[randperm(100,4)]},4,"randperm(100,4)");
{ my @p = 1..100;
  my @s = @p[randperm(0+@p)];
  isnt("@s", "@p", "randperm shuffle has shuffled input");
  my @ss = sort { $a<=>$b } @s;
  is("@ss", "@p", "randperm shuffle contains original data");
}

###### shuffle

is_deeply([shuffle()],[],"shuffle with no args");
is_deeply([shuffle("a")],["a"],"shuffle with one arg");
{ my @p = 1..100;
  my @s = shuffle(@p);
  is(0+@s,0+@p,"argument count is the same for 100 elem shuffle");
  isnt("@s", "@p", "shuffle has shuffled input");
  my @ss = sort { $a<=>$b } @s;
  is("@ss", "@p", "shuffle contains original data");
}

###### randperm, making sure we hit different algorithms
# TODO: Make versions for the PP code
SKIP: {
  skip "Skip set randperm tests with PP", 11 unless $usexs;
  csrand(15);  # arbitrarily chosen
  is_deeply([randperm(4,8)], [3,2,1,0], "randperm(4,8)");
  is_deeply([randperm(4,0)], [], "randperm(4,0)");
  is_deeply([randperm(42,1)], [35], "randperm(42,1)");
  is_deeply([randperm(2,2)], [1,0], "randperm(2,2)");
  is_deeply([randperm(2,2)], [0,1], "randperm(2,2)");
  is_deeply([randperm(1024,2)], [251,774], "randperm(1024,2)");
  is_deeply([randperm(75,6)], [74,19,0,24,33,11], "randperm(75,6)");
  is_deeply([randperm(123456789,37)], [2562439,85749631,40232376,82949734,102392910,100235558,85230001,27316639,65518150,64663313,73135111,90372409,112792441,89573412,46169851,106373363,112630472,20654195,111100168,48606846,88079756,109850048,61974647,123348106,121858399,60997453,105081822,84147357,81592378,11253701,24115274,70091937,22151205,69290886,52258596,118497340,10911837], "randperm(123456789,37)");
  is_deeply([randperm(54321,10)], [25785,17023,6559,15532,23397,32841,48735,15208,46324,2659], "randperm(54321,10)");
  is_deeply([randperm(30,12)], [15,7,11,4,26,2,13,25,6,17,14,19], "randperm(30,12)");
  is_deeply([randperm(16)], [3,10,7,1,5,0,8,14,2,9,13,6,4,12,15,11], "randperm(16)");
}

###### shuffling
SKIP: {
  skip "Skip set shuffle tests with PP", 3 unless $usexs;
  csrand(20); # arbitrary
  is_deeply([shuffle(qw<a b c d e f>)], [qw<d f c a b e>], "shuffle(a b c d e f)");
  is_deeply([shuffle(10..20)], [11,17,14,20,13,15,10,19,12,18,16], "shuffle(10..20)");
  is_deeply([shuffle(qw<a b b b b b c>)], [qw<b b b b a c b>], "shuffle(a b b b b b c)");
}
