#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/binomial factorial factorialmod
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
 [ 0,0, 1 ],
 [ 0,1, 0 ],
 [ 1,0, 1 ],
 [ 1,1, 1 ],
 [ 1,2, 0 ],
 [ 13,13, 1 ],
 [ 13,14, 0 ],
 [ 35,16, 4059928950 ],             # We can do this natively even in 32-bit
 [ 40,19, "131282408400" ],         # We can do this in 64-bit
 [ 67,31, "11923179284862717872" ], # ...and this
 [ 228,12, "30689926618143230620" ],# But the result of this is too big.
 [ 177,78, "3314450882216440395106465322941753788648564665022000" ],
 [ -10,5, -2002 ],
 [ -11,22, 64512240 ],
 [ -12,23, -286097760 ],
 [ -23,-26, -2300 ],     # Kronenburg extension
 [ -12,-23, -705432 ],   # same
 [  12,-23, 0 ],
 [  12,-12, 0 ],
 [ -12,0, 1 ],
 [  0,-1, 0 ],
);

# TODO: Add a bunch of combs here:  "5,3" => [[..],[..],[..]],

plan tests => 1                        # Factorial
            + 1 + 1*$extra + 2         # Factorialmod
            + 2 + scalar(@binomials)   # Binomial
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

###### binomial
foreach my $r (@binomials) {
  my($n, $k, $exp) = @$r;
  is( binomial($n,$k), $exp, "binomial($n,$k)) = $exp" );
}
is_deeply( [map { binomial(10, $_) } -15 .. 15],
           [qw/0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 10 45 120 210 252 210 120 45 10 1 0 0 0 0 0/],
           "binomial(10,n) for n in -15 .. 15" );
is_deeply( [map { binomial(-10, $_) } -15 .. 15],
           [qw/-2002 715 -220 55 -10 1 0 0 0 0 0 0 0 0 0 1 -10 55 -220 715 -2002 5005 -11440 24310 -48620 92378 -167960 293930 -497420 817190 -1307504/],
           "binomial(-10,n) for n in -15 .. 15" );

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
