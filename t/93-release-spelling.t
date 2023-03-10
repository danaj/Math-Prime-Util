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

add_stopwords(qw/-th
                 bigint bigints
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
                 forsemiprimes forfactored forsquarefree foralmostprimes
                 lastfor
                 numtoperm permtonum randperm
                 totient moebius mertens liouville kronecker znorder znprimroot znlog
                 sumliouville sumtotient
                 gcd lcm gcdext chinese
                 sqrtmod allsqrtmod rootmod allrootmod
                 invmod addmod submod mulmod powmod divmod
                 binomialmod factorialmod
                 bernfrac bernreal harmfrac harmreal stirling hclassno
                 vecsum vecprod vecmin vecmax vecreduce vecextract vecequal
                 vecall vecany vecnone vecnotall vecfirst vecfirstidx vecmex vecpmex
                 sqrtint logint rootint powint addint subint mulint divint modint negint absint divrem tdivrem fdivrem
                 qnr
                 todigits todigitstring fromdigits sumdigits hammingweight
                 tozeckendorf fromzeckendorf
                 lucasu lucasv lucasuv lucasumod lucasvmod lucasuvmod
                 lshiftint rshiftint rashiftint signint cmpint
                 biquadrate powerfree k-powerfree k-free
                 pp/);

all_pod_files_spelling_ok();
