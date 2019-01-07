#!/usr/bin/env perl
use strict;
use warnings;

BEGIN {
  unless ($ENV{RELEASE_TESTING}) {
    require Test::More;
    Test::More::plan(skip_all => 'these tests are for release candidate testing');
  }
}

#---------------------------------------------------------------------


use Test::More;
eval "use Test::Spellunker";
plan skip_all => "Test::Spellunker required for testing POD spelling" if $@;

add_stopwords(qw/bigint bigints
                 bignum bignums
                 quadmath
                 pseudoprime pseudoprimes
                 primorial primorials
                 semiprime semiprimes
                 precalculated premultiplier
                 benchmarking hardcoded online
                 unoptimized unusably
                 coprime summatory
                 RiemannR LambertW
                 csrand srand irand irand64 drand urandomb urandomm
                 forprimes forcomposites foroddcomposites fordivisors
                 forpart forcomp forcomb forperm forderange formultiperm forsetproduct
                 forsemiprimes forfactored forsquarefree
                 lastfor
                 numtoperm permtonum randperm
                 totient moebius mertens liouville kronecker znorder znprimroot znlog
                 gcd lcm gcdext chinese invmod sqrtmod addmod mulmod powmod divmod
                 bernfrac bernreal harmfrac harmreal stirling hclassno
                 vecsum vecprod vecmin vecmax vecreduce vecextract
                 vecall vecany vecnone vecnotall vecfirst vecfirstidx
                 sqrtint logint rootint powint addint mulint divint modint divrem tdivrem
                 factorialmod
                 todigits todigitstring fromdigits sumdigits hammingweight
                 lucasu lucasv
                 pp/);

all_pod_files_spelling_ok();
