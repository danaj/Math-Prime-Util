#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/forcomb forperm forderange formultiperm
                         numtoperm permtonum/;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $usexs  = Math::Prime::Util::prime_get_config->{'xs'};
my $usegmp = Math::Prime::Util::prime_get_config->{'gmp'};

plan tests => 0
            + 1  # Combinations
            + 1  # Permutations
            + 1  # Multiset Permutations
            + 1  # Derangements
            + 1  # numtoperm, permtonum
            ;

subtest 'forcomb', sub {
  my @t1data = ( [ 0, "[]" ],
                 [ 1, "[] [0]" ],
                 [ 2, "[] [0] [1] [0 1]" ],
                 [ 3, "[] [0] [1] [2] [0 1] [0 2] [1 2] [0 1 2]" ],
               );
  my @t2data = ( [ 0,0, "[]" ],
                 [ 5,6, "" ],
                 [ 5,5, "[0 1 2 3 4]" ],
                 [ 5,4, "[0 1 2 3] [0 1 2 4] [0 1 3 4] [0 2 3 4] [1 2 3 4]" ],
                 [ 5,3, "[0 1 2] [0 1 3] [0 1 4] [0 2 3] [0 2 4] [0 3 4] [1 2 3] [1 2 4] [1 3 4] [2 3 4]" ],
                 [ 5,2, "[0 1] [0 2] [0 3] [0 4] [1 2] [1 3] [1 4] [2 3] [2 4] [3 4]" ],
                 [ 5,1, "[0] [1] [2] [3] [4]" ],
                 [ 5,0, "[]" ],
                 [ 4,4, "[0 1 2 3]" ],
                 [ 4,3, "[0 1 2] [0 1 3] [0 2 3] [1 2 3]" ],
                 [ 4,2, "[0 1] [0 2] [0 3] [1 2] [1 3] [2 3]" ],
                 [ 4,1, "[0] [1] [2] [3]" ],
                 [ 4,0, "[]" ],
               );

  for my $test (@t1data) {
    my($e,$exp) = @$test;
    my @p;  forcomb { push @p, "[@_]" } $e;
    is("@p", $exp, "forcomb {} $e");
  }
  for my $test (@t2data) {
    my($s,$e,$exp) = @$test;
    my @p;  forcomb { push @p, "[@_]" } $s,$e;
    is("@p", $exp, "forcomb {} $s,$e");
  }
  {
    my @data = (qw/apple bread curry/);
    my @p;
    forcomb { push @p, [@data[@_]] } @data,2;
    is_deeply(\@p, [[qw/apple bread/],[qw/apple curry/],[qw/bread curry/]], "forcomb {} 3,2");
  }
  {
    my @data = (qw/ant bee cat dog/);
    my @p;
    forcomb { push @p, [@data[@_]] } @data,3;
    is_deeply(\@p, [[qw/ant bee cat/],[qw/ant bee dog/],[qw/ant cat dog/],[qw/bee cat dog/]], "forcomb 4,3");
  }
  {
    my $s = 0;
    forcomb { $s++ } 20,15;
    is($s, 15504, "forcomb { } 20,15 yields binomial(20,15) combinations");
  }
};


subtest 'forperm', sub {
  my %perms = (
   0 => [[]],
   1 => [[0]],
   2 => [[0,1],[1,0]],
   3 => [[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,0,1],[2,1,0]],
   4 => [[0,1,2,3],[0,1,3,2],[0,2,1,3],[0,2,3,1],[0,3,1,2],[0,3,2,1],[1,0,2,3],[1,0,3,2],[1,2,0,3],[1,2,3,0],[1,3,0,2],[1,3,2,0],[2,0,1,3],[2,0,3,1],[2,1,0,3],[2,1,3,0],[2,3,0,1],[2,3,1,0],[3,0,1,2],[3,0,2,1],[3,1,0,2],[3,1,2,0],[3,2,0,1],[3,2,1,0]],
  );
  while (my($n, $expect) = each (%perms)) {
    my @p; forperm { push @p, [@_] } $n;
    is_deeply(\@p, $expect, "forperm $n");
  }
  {
    my $s = 0;
    forperm { $s++ } 7;
    is($s, 5040, "forperm 7 yields factorial(7) permutations");
  }
};

subtest 'formultiperm', sub {
  { my @p; formultiperm { push @p, [@_] } [];
    is_deeply(\@p, [], "formultiperm []"); }

  { my @p; formultiperm { push @p, [@_] } [1,2,2];
    is_deeply(\@p, [ [1,2,2], [2,1,2], [2,2,1] ], "formultiperm 1,2,2"); }

  { my @p; formultiperm { push @p, [@_] } [qw/a a b b/];
    is_deeply(\@p, [map{[split(//,$_)]} qw/aabb abab abba baab baba bbaa/], "formultiperm a,a,b,b"); }

  { my @p; formultiperm { push @p, join("",@_) } [qw/a a b b/];
    is_deeply(\@p, [qw/aabb abab abba baab baba bbaa/], "formultiperm aabb"); }
};

subtest 'forderange', sub {
  { my @p; forderange { push @p, [@_] } 0;
    is_deeply(\@p, [[]], "forderange 0"); }
  { my @p; forderange { push @p, [@_] } 1;
    is_deeply(\@p, [], "forderange 1"); }
  { my @p; forderange { push @p, [@_] } 2;
    is_deeply(\@p, [[1,0]], "forderange 2"); }
  { my @p; forderange { push @p, [@_] } 3;
    is_deeply(\@p, [[1,2,0],[2,0,1]], "forderange 3"); }
  { my $n=0; forderange { $n++ } 7; is($n, 1854, "forderange 7 count"); }
};

subtest 'numtoperm / permtonum', sub {
  is_deeply([numtoperm(0,0)],[],"numtoperm(0,0)");
  is_deeply([numtoperm(1,0)],[0],"numtoperm(1,0)");
  is_deeply([numtoperm(1,1)],[0],"numtoperm(1,1)");
  is_deeply([numtoperm(5,15)],[0,3,2,4,1],"numtoperm(5,15)");
  is_deeply([numtoperm(24,987654321)],[0,1,2,3,4,5,6,7,8,9,10,13,11,21,14,20,17,15,12,22,18,19,23,16],"numtoperm(24,987654321)");

  is(permtonum([]),0,"permtonum([])");
  is(permtonum([0]),0,"permtonum([0])");
  is(permtonum([6,3,4,2,5,0,1]),4768,"permtonum([6,3,4,2,5,0,1])");
  is(permtonum([reverse(0..14),15..19]),"1790774578500738480","permtonum( 20 )");
  is("".permtonum([reverse(0..12),reverse(13..25)]),"193228515634198442606207999","permtonum( 26 )");

  is(permtonum([numtoperm(14,8467582)]),8467582,"permtonum(numtoperm)");
};
