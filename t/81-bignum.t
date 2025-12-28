#!/usr/bin/env perl
use strict;
use warnings;

# If you're not using ancient perl 5.6.2 with super early releases of bigint,
# then you can define bigint up here and not have to quote every number.

# Note: In 5.16.0 (and perhaps others?), using labels like "SKIP:" will create
# a small memory leak.  So running the test suite through valgrind will show
# some small leaks in this test, which has nothing to do with the module.

my $extra = 0+(defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING});
my $use64 = ~0 > 4294967295;
my $broken64 = (18446744073709550592 == ~0);

use Test::More;

# Performance issues here are in:
#   primality        very little we can do about it, we must test
#   factoring        also must test though maybe can make faster?

if ($broken64) {
  # Perl from before 2002, built for 64-bit.  Not supported.
  plan skip_all => "Broken 64-bit Perl, skipping all tests";
} else {
  plan tests =>  0
               + 1   # basic int
               + 1   # basic mod
               + 1   # other mod
               + 1   # gcd and lcm
               + 1   # gcdext and chinese
               + 1   # primality
               + 1   # primes, twin primes, semiprimes, almost primes, etc.
               + 1   # next/prev prime
               + 1   # prime_iterator
               + 1   # pseudoprimes
               + 1   # primecount and lower/upper/approx
               + 1   # factoring
               + 1   # znorder znprimroot znlog
               + 1   # divisor_sum
               + 5   # moebius, euler_phi, kronecker, valuation, etc.
               + 1   # jordan_totient
               + 1   # liouville
               + 1   # ispower
               + 1   # random primes
               + 1   # Perrin PsP
               + 1   # vecequal
               + 1;  # $_ didn't get changed
}

# Using GMP makes these tests run about 2x faster on some machines
#use bigint try => 'GMP,GMPz,LTM,Pari'; #  <-------- large numbers ahead!  > 2^64
use bigint;
use Math::BigFloat;

use Math::Prime::Util qw/
  prime_set_config
  negint absint signint cmpint
  addint subint add1int sub1int
  mulint powint lshiftint rshiftint rashiftint
  divint modint cdivint divrem fdivrem cdivrem tdivrem
  sqrtint rootint logint
  negmod addmod submod mulmod muladdmod mulsubmod
  divmod powmod invmod sqrtmod rootmod allsqrtmod allrootmod
  is_congruent is_qr qnr is_primitive_root
  factorialmod binomialmod
  lucasumod lucasvmod
  is_prob_prime
  prime_count_lower
  prime_count_upper
  prime_count_approx
  nth_prime_lower
  nth_prime_upper
  nth_prime_approx
  factor
  factor_exp
  divisors
  kronecker
  moebius
  euler_phi
  carmichael_lambda
  jordan_totient
  divisor_sum
  znorder
  znprimroot
  znlog
  liouville
  gcd
  lcm
  gcdext
  chinese
  is_power
  pn_primorial
  ExponentialIntegral
  LogarithmicIntegral
  RiemannR
  primes
  twin_primes
  prime_count
  nth_prime
  is_prime
  is_provable_prime
  next_prime
  prev_prime
  prime_iterator
  is_strong_pseudoprime
  random_prime
  random_ndigit_prime
  random_nbit_prime
  random_safe_prime
  random_strong_prime
  random_maurer_prime
  miller_rabin_random
  is_perrin_pseudoprime
  is_bpsw_prime
  verify_prime
  valuation
  vecequal
/;
# TODO:  is_strong_lucas_pseudoprime
#        ExponentialIntegral
#        LogarithmicIntegral
#        RiemannR

my $usegmp = Math::Prime::Util::prime_get_config->{gmp};
my $bignumver = $bigint::VERSION;
my $bigintver = $Math::BigInt::VERSION;
my $bigintlib = Math::BigInt->config()->{lib};
   $bigintlib =~ s/^Math::BigInt:://;
my $mpugmpver = $usegmp ? $Math::Prime::Util::GMP::VERSION : "<none>";
diag "BigInt $bignumver/$bigintver, lib: $bigintlib.  MPU::GMP $mpugmpver\n";

# Turn off use of BRS - ECM tries to use this.
# prime_set_config( irand => sub { int(rand(4294967296)) } );


###############################################################################

$_ = 'this should not change';

subtest 'arithmetic ops', sub {
  # The test numbers used are randomly chosen, nothing special about them
  my $n = 29326327445836963809762302067795652280;
  my $negn = negint($n);
  is($negn, -29326327445836963809762302067795652280, "negint(n)");
  is(negint($negn), $n, "negint(negint(n))");
  is(absint($negn), $n, "absint(negint(n))");
  is(signint($n), 1, "signint(n)");
  is(signint($negn), -1, "signint(negint(n))");
  is(cmpint($n,$n),0,"cmpint(n,n) = 0");
  is(cmpint($n,$negn),1,"cmpint(n,-n) = 1");
  is(cmpint($negn,$n),-1,"cmpint(-n,n) = -1");

  my $n2 = 58652654891673927619524604135591304560;
  my $n3 = 87978982337510891429286906203386956840;
  is(addint($n,$n),$n2,"addint(n,n) = 2n");
  is(addint($n,$n2),$n3,"addint(2*n,n) = 3n");
  is(addint($negn,$n),"0","addint(-n,n) = 0");
  is(subint($n3,$n2),$n,"subint(3n,2n) = n");
  is(subint($n,$negn),$n2,"subint(n,-n) = 2n");
  is(add1int(19437056586806593268),19437056586806593269,"add1int(n) = n+1");
  is(sub1int(19437056586806593268),19437056586806593267,"sub1int(n) = n-1");

  my $nsqr = 860033481460450377515108424533203603187404404374382196284947588070669198400;
  my $ncub = 25221623491692321487348209824650128345904008945462108558207852225847625066929594776910302455742481841632732352000;
  is(mulint(2,$n),$n2,"mulint(2,n) = 2n");
  is(mulint($n,3),$n3,"mulint(n,3) = 3n");
  is(mulint($n,-3),negint($n3),"mulint(n,-3) = -3n");
  is(mulint($n,$n),$nsqr,"mulint(n,n) = n^2");
  is(powint($n,0),"1","powint(n,0) = 1");
  is(powint($n,1),$n,"powint(n,1) = n");
  is(powint($n,2),$nsqr,"powint(n,2) = n^2");
  is(powint($n,3),$ncub,"powint(n,3) = n^3");
  is(powint($negn,0),"1","powint(-n,0) = 1");
  is(powint($negn,1),$negn,"powint(-n,1) = -n");
  is(powint($negn,2),$nsqr,"powint(-n,2) = n^2");
  is(powint($negn,3),negint($ncub),"powint(-n,3) = -n^3");

  is(lshiftint($n),$n2,"lshiftint(n) = 2n");
  is(lshiftint($n3,14),1441447646617778445177436671236291900866560,"lshiftint(3n,14) = 3n * 2^14");
  is(rshiftint($n2),$n,"rshiftint(2n) = n");
  is(rshiftint($n,7),229111933170601279763767984904653533,"rshiftint(n,7) = n / 2^7");
  is(rshiftint($negn,7),-229111933170601279763767984904653533,"rshiftint(-n,7) = -(n >> 7)");
  is(rashiftint($negn,7),-229111933170601279763767984904653534,"rashiftint(-n,7) = (fdivrem(n,2**7))[0] [Python right shift]");

  is(divint($n2,2),$n,"divint(2n,2) = n");
  is(divint(negint($n3),3),$negn,"divint(-3n,3) = -n");
  is(modint($n3,$n),"0","mod(3n,n) = 0");
  is(modint($n,1),"0","mod(n,1) = 0");
  is(modint($n,29),"5","mod(n,29)");
  is(modint($negn,37),"30","mod(-n,37)");
  is(divint($negn,511),-57390073279524391017147362167897559,"divint(-n,511)");
  is(cdivint($negn,511),-57390073279524391017147362167897558,"cdivint(-n,511)");

  my $m = -891149066616685318069719273514027361488599736308638170770;
  my $d = 11081041228439474651;
  is_deeply([divrem($m,$d)],
            [-80421058657335617430773462550803594698,5116452384810829628],
            "divrem(-m,d)");
  is_deeply([fdivrem($m,$d)],
            [-80421058657335617430773462550803594698,5116452384810829628],
            "fdivrem(-m,d)");
  is_deeply([cdivrem($m,$d)],
            [-80421058657335617430773462550803594697,-5964588843628645023],
            "cdivrem(-m,d)");
  is_deeply([tdivrem($m,$d)],
            [-80421058657335617430773462550803594697,-5964588843628645023],
            "tdivrem(-m,d)");

  is(sqrtint($nsqr),$n,"sqrtint(n^2) = n");
  is(sqrtint($ncub),158813171656800310523242313199648286862899164761029420165,"sqrtint(n^3) = n^(3/2)");

  {
    my($r,$v,$e);
    $r = rootint($ncub, 7, \$v);
    $e = 11412903582751672;
    is_deeply([$r,$v], [$e,powint($e,7)], "rootint(n^3,7,\\\$r)");

    $r = logint($n, 18, \$v);
    $e = 29;
    is_deeply([$r,$v], [$e,powint(18,$e)], "logint(n,18,\\\$r)");
  }
};

subtest 'basic mod ops', sub {
  # Randomly chosen integers
  my $A  = 400189273193594088451869468713181466450681802689782908252400244075;
  my $B  = 990429319976107635118049598881942833305583;
  my $C  = 130716934118556849942;
  my $Ca = 130716934118556849953;  # a composite
  my $D  = 95045564778113305157;   # a prime
  my $Da = 95045564778113305482;   # a non-squarefree composite

  is(negmod($A,$D),12736484928456701970,"negmod");
  is(addmod($A,$B,$C),15656969117801514554,"addmod");
  is(submod($A,$B,$C),74721010532486795616,"submod");
  is(mulmod($A,$B,$C),4443724024613552983,"mulmod");
  is(muladdmod($A,$B,$D,$C),99489288802726858140,"muladdmod");
  is(mulsubmod($A,$B,$D,$C),40115093365057097768,"mulsubmod");
  is(divmod($A,$B,$C),71001221878922773879,"divmod");
  is(powmod($A,$B,$C),80999932714432421083,"powmod");
  is(invmod($A,$C),33944411285878109917,"invmod");

  is(sqrtmod($C,$D),8765086989354512719,"sqrtmod [prime modulus]");
  my @all1 = (13353174708808101300,19189233687021677589,111527700431535172364,117363759409748748653);
  my @all2 = (2960938399489571940,22255231847764978494,72790332930348326988,92084626378623733542);
  ok(is_one_of(sqrtmod($D,$Ca),@all1), "sqrtmod [composite modulus]");
  ok(is_one_of(sqrtmod($C,$Da),@all2), "sqrtmod [composite modulus]");

  is(rootmod($C,5,$D),81379717162684976939,"rootmod(C,5,D) [prime modulus]");
  my @all3 =(1653087340990329197,15231025166435087123,28808962991879845049,42386900817324602975,55964838642769360901,69542776468214118827,83120714293658876753);
  ok(is_one_of(rootmod($D,35,$Da),@all3),"rootmod(C,35,Da) [composite modulus]");

  is_deeply([allsqrtmod($C,$D)],[8765086989354512719,86280477788758792438],"allsqrtmod");
  is_deeply([allrootmod($D,35,$Da)], \@all3, "allrootmod");
};

subtest 'other mod ops', sub {
  my $A  = 400189273193594088451869468713181466450681802689782908252400244075;
  my $B  = 990429319976107635118049598881942833305583;
  my $C  = 130716934118556849940;

  is(is_congruent($A,$B,$C),0,"is_congruent");
  is(is_congruent(addint($A,20208865721792522128),$B,$C),1,"is_congruent");

  is(is_qr($B,$C),0,"is_qr");
  is(is_qr(addint($B,142),$C),1,"is_qr");
  is(qnr($C), 2, "qnr");
  is(qnr(130716934118556861071), 37, "qnr");

  my $N = 72574551534231909331741171093173785967490646405143;  # 7^59
  is(is_primitive_root(3,$N),1,"is_primitive_root");
  is(is_primitive_root(7,$N),0,"is_primitive_root");

  is(factorialmod(777, $C), 27183400742691213580, "factorialmod");
  is(factorialmod(36893488147432436565, $C), "0", "factorialmod");

  is(binomialmod($C,73,$B),8314965924716374520974024285356149070687,"binomialmod");

  # lucasumod lucasvmod lucasuvmod done in t/25-lucas_sequences.t
  is(lucasumod(17,1,7777777778888999,$B),363910679494422198934676182699576192782510,"lucasumod with Q=1");
  is(lucasvmod(17,1,7777777778888999,$B),68398397794084099557151356283528547714286,"lucasvmod with Q=1");
};

subtest 'gcd and lcm', sub {
  is( gcd(921166566073002915606255698642,1168315374100658224561074758384,951943731056111403092536868444), 14, "gcd(a,b,c)" );
  is( gcd(1214969109355385138343690512057521757303400673155500334102084,1112036111724848964580068879654799564977409491290450115714228), 42996, "gcd(a,b)" );
  is( gcd(745845206184162095041321,61540282492897317017092677682588744425929751009997907259657808323805386381007), 1, "gcd of two primes = 1" );
  is( "".lcm(9999999998987,10000000001011), 99999999999979999998975857, "lcm(p1,p2)" );
  is( "".lcm(892478777297173184633,892478777297173184633), 892478777297173184633, "lcm(p1,p1)" );
  is( "".lcm(23498324,32497832409432,328732487324,328973248732,3487234897324), 1124956497899814324967019145509298020838481660295598696, "lcm(a,b,c,d,e)" );
};

subtest 'gcdext and chinese', sub {
  # Work around a Math::BigInt::Lite issue
  is_deeply( [gcdext(803028077895224634710139483024654235947,101394830246542359478030280778952246347)], [7687627417944666569835322284775478836, -60884570288210047004733169112173096587, "3"], "gcdext(a,b)" );
  is( chinese([26,17179869209],[17,34359738421]), 103079215280, "chinese([26,17179869209],[17,34359738421] = 103079215280" );
};

###############################################################################

subtest 'primality', sub {
  my @primes = (100000982717289000001);
  foreach my $n (@primes) {
    ok( is_prime($n), "$n is prime" );
    ok( is_prob_prime($n), "$n is probably prime");
  }

  my @composites = (
    36893488147419103233,   # div 3
    36893488147419103249,   # div 7
    36893488147419103261,   # div 61
    36893488147419103253,   # no small factors
    21652684502221,         # small pseudoprime
    1195068768795265792518361315725116351898245581,  # big pseudoprime
  );
  foreach my $n (@composites) {
    #ok( !is_prime($n), "$n is not prime" );
    ok( !is_prob_prime($n), "$n is not probably prime");
  }

  # Primes where n-1 is easy to factor, so we finish quickly.
  my @proveprimes = qw/
    65635624165761929287 1162566711635022452267983
    77123077103005189615466924501 3991617775553178702574451996736229
    273952953553395851092382714516720001799
  /;
  foreach my $n (@proveprimes) {
    #ok( is_prime($n), "$n is prime" );
    SKIP: {
      skip "Large proof on 32-bit machine without EXTENDED_TESTING.", 1
        if !$use64 && !$extra && $n > 2**66;
      skip "Large proof without GMP or EXTENDED_TESTING.", 1
        if !$usegmp && !$extra && $n > 2**66;
      ok( is_provable_prime($n), "$n is provably prime" );
    }
  }
};

###############################################################################

subtest 'range primes', sub {
  # Used to be (2**66, 2**66+100), but old ActiveState chokes for some reason.
  is_deeply( primes(73786976294838206464, 73786976294838206564),
             [73786976294838206473,73786976294838206549],
             "primes( 2^66, 2^66 + 100 )" );

  is_deeply( twin_primes(18446744073709558000, 18446744073709559000),
             [18446744073709558601,18446744073709558727],
             "twin_primes( 18446744073709558000, +1000)" );
};

###############################################################################

subtest 'next and prev', sub {
  is( next_prime(777777777777777777777777), 777777777777777777777787, "next_prime(777777777777777777777777)");
  is( prev_prime(777777777777777777777777), 777777777777777777777767, "prev_prime(777777777777777777777777)");
};

###############################################################################

subtest 'prime iterator', sub {
  my $it = prime_iterator(10**24+910);
  is_deeply( [map { $it->() } 1..3],
             [1000000000000000000000921,1000000000000000000000931,1000000000000000000000949],
             "iterator 3 primes starting at 10^24+910" );
};

###############################################################################

subtest 'pseudoprimes', sub {
  # composite pseudoprimes to various small prime bases
  my %pseudoprimes = (
   '75792980677' => [ qw/2/ ],
   '21652684502221' => [ qw/2 7 37 61 9375/ ],
   '3825123056546413051' => [ qw/2 3 5 7 11 13 17 19 23 29 31 325 9375/ ],
   '318665857834031151167461' => [ qw/2 3 5 7 11 13 17 19 23 29 31 37 325 9375/ ],
   '3317044064679887385961981' => [ qw/2 3 5 7 11 13 17 19 23 29 31 37 73 325 9375/ ],
   '6003094289670105800312596501' => [ qw/2 3 5 7 11 13 17 19 23 29 31 37 61 325 9375/ ],
   '59276361075595573263446330101' => [ qw/2 3 5 7 11 13 17 19 23 29 31 37 325 9375/ ],
   '564132928021909221014087501701' => [ qw/2 3 5 7 11 13 17 19 23 29 31 37 325 9375/ ],
   #'1543267864443420616877677640751301' => [ qw/2 3 5 7 11 13 17 19 23 29 31 37 61 325 9375/ ],
  );

  while (my($psrp, $baseref) = each (%pseudoprimes)) {
    my $baselist = join(",", @$baseref);
    my @expmr = map { (0!=1) } @$baseref;
    my @gotmr = map { is_strong_pseudoprime($psrp, $_) } @$baseref;
    is_deeply(\@gotmr, \@expmr, "$psrp is a strong pseudoprime to bases $baselist");
  }
};

###############################################################################

subtest 'prime counts', sub {
  # Testing prime_count only on a small range for time reasons
  is( prime_count(877777777777777777777752, 877777777777777777777872), 2, "prime_count(87..7752, 87..7872)");

  # test bounds
  check_pcbounds(31415926535897932384, 716115441142294636, '2e-8', '2e-8');
  if ($extra) {
    check_pcbounds(314159265358979323846, 6803848951392700268, '5e-9', '5e-9');
    check_pcbounds(31415926535897932384626433, 544551456607147153724423, '3e-6', '3e-11');
    # pi(10^23) = 1925320391606803968923
    check_pcbounds(10**23, 1925320391607837268776, '5e-10', '5e-10');
  }
};

###############################################################################

subtest 'factoring', sub {
  my %factors = (
    '1234567890' => [2, 3, 3, 5, 3607, 3803],
    '190128090927491' => [61, 73, 196291, 217517],
    '23489223467134234890234680' => [2, 2, 2, 5, 4073, 4283, "33662485846146713"],
    #'7674353466844111807691499613711' => [11783, 12239, 18869, 22277, 37861, 55163, 60617],
  );
  my %allfactors = (
    '23489223467134234890234680' => [qw/1 2 4 5 8 10 20 40 4073 4283 8146 8566 16292 17132 20365 21415 32584 34264 40730 42830 81460 85660 162920 171320 17444659 34889318 69778636 87223295 139557272 174446590 348893180 697786360 33662485846146713 67324971692293426 134649943384586852 168312429230733565 269299886769173704 336624858461467130 673249716922934260 1346499433845868520 137107304851355562049 144176426879046371779 274214609702711124098 288352853758092743558 548429219405422248196 576705707516185487116 685536524256777810245 720882134395231858895 1096858438810844496392 1153411415032370974232 1371073048513555620490 1441764268790463717790 2742146097027111240980 2883528537580927435580 5484292194054222481960 5767057075161854871160 587230586678355872255867 1174461173356711744511734 2348922346713423489023468 2936152933391779361279335 4697844693426846978046936 5872305866783558722558670 11744611733567117445117340 23489223467134234890234680/],
  );

  while (my($n, $factors) = each(%factors)) {
    is_deeply( [factor($n)], $factors, "factor($n)" );
    is_deeply( [factor_exp($n)], [linear_to_exp(@$factors)], "factor_exp($n)" );
  }
  while (my($n, $allfactors) = each(%allfactors)) {
    is_deeply( [divisors($n)], $allfactors, "divisors($n)" );
  }
};

###############################################################################

subtest 'znorder znprimroot znlog', sub {
  # Calc/FastCalc are slugs with this function, so tone things down.
  #is( znorder(82734587234,927208363107752634625923555185111613055040823736157),
  #    4360156780036190093445833597286118936800,
  #    "znorder" );
  is("".znorder(8267,927208363107752634625925),2838011904800209433220,"znorder(8267,927208363107752634625925)");
  is(znorder(902,827208363107752634625947),undef,"znorder(902,827208363107752634625947");

  is( znprimroot(333822190384002421914469856494764513809), 3, "znprimroot(333822190384002421914469856494764513809)" );

  is( znlog(232752345212475230211680, 23847293847923847239847098123812075234, 804842536444911030681947), 13, "znlog(b,g,p): find k where b^k = g mod p" );
};

###############################################################################

subtest 'divisor sum', sub {
  # Done wrong, the following will have a bunch of extra zeros.
  my $hundredfac = Math::BigInt->new(100)->bfac;
  is( "".divisor_sum($hundredfac), 774026292208877355243820142464115597282472420387824628823543695735957009720184359087194959566149232506852422409529601312686157396490982598473425595924480000000, "Divisor sum of 100!" );
  # These should yield bigint results.
  # Quoted 0 to prevent error in perl 5.8.2 + bigint 0.23 (0 turns into NaN)
  is( divisor_sum(pn_primorial(27),"0"), 134217728, "Divisor count(103#)" );
  is( "".divisor_sum(pn_primorial(27),1), 123801167235014219383860918985791897600000, "Divisor sum(103#)" );
  is( "".divisor_sum(pn_primorial(27),2), 872887488619258559049272439859735080160421720974947767918289356800000000000000000, "sigma_2(103#)" );
  if ($extra) {
    is( "".divisor_sum(pn_primorial(71),"0"), 2361183241434822606848, "Divisor count(353#)" );
  }
};

###############################################################################

{
  my $n = 618970019642690137449562110;
  is( moebius($n), -1, "moebius($n)" );
  is( "".euler_phi($n), 145857122964987051805507584, "euler_phi($n)" );
  is( carmichael_lambda($n), 3271601336256, "carmichael_lambda($n)" );
  is( kronecker(878944444444444447324234,216539985579699669610468715172511426009), -1, "kronecker(..., ...)" );
  is( valuation(6**10000-1,5), 5, "valuation(6^10000,5) = 5" );
}

subtest 'jordan totient', sub {
  my $n = 2188536338969724335815;
  is( "".jordan_totient(5,$n), 50191437962401621299442160569980637495164094052874237077859453863016726060686241335185072669341771886684800, "jordan_totient(5,$n)" );
  $n = 2188536338969724335807;
  is( divisor_sum( $n, sub { my $d=shift; $d**5 * moebius($n/$d); }), 50207524710890617788554288878260755791080217791665431423557510096680804997771551711694188532723268222129800, "jordan totient using divisor_sum and moebius" );

  if ($extra) {
    my $n = 48981631802481400359696467;
    is( jordan_totient(5,$n), 281946200770875813001683560563488308767928594805846855593191749929654015729263525162226378019837608857421063724603387506651820000, "jordan_totient(5,$n)" );
    is( divisor_sum( $n, sub { my $d=shift; $d**5 * moebius($n/$d); }), "281946200770875813001683560563488308767928594805846855593191749929654015729263525162226378019837608857421063724603387506651820000", "jordan totient using divisor_sum and moebius" );
  }
};

subtest 'liouville', sub {
  is( liouville(  560812147176208202656339069),"-1", "liouville(a x b x c) = -1" );
  is( liouville(10571644062695614514374497899), "1", "liouville(a x b x c x d) = 1" );
};
###############################################################################

subtest 'is_power', sub {
  is( is_power(18475335773296164196), "0", "ispower(18475335773296164196) == 0" );
  is( is_power(3089265681159475043336839581081873360674602365963130114355701114591322241990483812812582393906477998611814245513881), 14, "ispower(150607571^14) == 14" );
  my @negpowers = (qw/0 0 0 3 0 5 3 7 0 9 5 11 3 13 7 15 0 17 9 19 5 21 11 23 3 25 13 27 7 29 15 31/);
  push @negpowers, (qw/0 33 17 35 9 37 19 39 5 41 21 43 11 45 23 47 3 49 25 51 13 53 27 55 7 57 29 59 15 61 31 63 0 65 33 67 17 69 35 71 9 73 37 75 19 77 39 79 5 81 41 83 21 85 43 87 11 89 45 91 23 93 47 95 3 97 49 99 25 101 51 103 13 105 53 107 27 109 55 111 7 113 57 115 29 117 59 119 15 121 61 123 31 125 63 127 0 129 65 131 33 133 67 135 17 137 69 139 35 141 71 143 9 145 73 147 37 149 75/) if $extra;
  # Work around bug in Math::BigInt::Pari and Perl pre-5.18.
  if ($bigintlib eq 'Pari' && $] < "5.018") {
    is_deeply( [map { is_power("".-7 ** $_) } int("0") .. $#negpowers], \@negpowers, "-7 ^ i for 0 .. $#negpowers" );
    is_deeply( [map { my $r; my $p=is_power("".-7 ** $_, "0", \$r); $p ? (0+$r) ** $p : -7 ** $_; } int("0") .. $#negpowers], [map { -7 ** $_ } int("0") .. $#negpowers], "correct root from is_power for -7^i for 0 .. $#negpowers" );
  } else {
    is_deeply( [map { is_power(-7 ** $_) } int("0") .. $#negpowers], \@negpowers, "-7 ^ i for 0 .. $#negpowers" );
    is_deeply( [map { my $r; my $p=is_power(-7 ** $_, "0", \$r); $p ? (1*$r) ** $p : -7 ** $_; } int("0") .. $#negpowers], [map { -7 ** $_ } int("0") .. $#negpowers], "correct root from is_power for -7^i for 0 .. $#negpowers" );
  }
};

###############################################################################

# We have a separate test file for random primes (t/16-randomprime.t)
# These tests really should just all be moved there.
subtest 'random primes', sub {
  my $randprime;

  SKIP: {  # \/ \/ \/   skipped without $extra
  skip "Skipping random prime tests without EXTENDED_TESTING", 26 unless $extra;

  $randprime = random_prime(147573952590750158861, 340282366920939067930896100764782952647);
  cmp_ok( $randprime, '>=', 147573952590750158861, "random range prime isn't too small");
  cmp_ok( $randprime, '<=', 340282366920939067930896100764782952647, "random range prime isn't too big");
  ok( is_prime($randprime), "random range prime is prime");

  $randprime = random_ndigit_prime(25);
  cmp_ok( $randprime, '>', 10**24, "random 25-digit prime is not too small");
  cmp_ok( $randprime, '<', 10**25, "random 25-digit prime is not too big");
  ok( is_prime($randprime), "random 25-digit prime is just right");

  $randprime = random_nbit_prime(80);
  cmp_ok( $randprime, '>', 2**79, "random 80-bit prime is not too small");
  cmp_ok( $randprime, '<', 2**80, "random 80-bit prime is not too big");
  ok( is_prime($randprime), "random 80-bit prime is just right");

  # This routine is very slow without GMP
  $randprime = random_safe_prime(100);
  cmp_ok( $randprime, '>', 2**99, "random 100-bit safe prime is not too small");
  cmp_ok( $randprime, '<', 2**100, "random 100-bit safe prime is not too big");
  ok( is_prime($randprime), "random 100-bit safe prime is just right");

  $randprime = random_strong_prime(180);
  cmp_ok( $randprime,'>',2**179,"random 180-bit strong prime is not too small");
  cmp_ok( $randprime,'<',2**180,"random 180-bit strong prime is not too big");
  ok( is_prime($randprime), "random 180-bit strong prime is just right");

  $randprime = random_maurer_prime(80);
  cmp_ok( $randprime,'>',2**79,"random 80-bit Maurer prime is not too small");
  cmp_ok( $randprime,'<',2**80,"random 80-bit Maurer prime is not too big");
  ok( is_prime($randprime), "random 80-bit Maurer prime is just right");

  $randprime = random_nbit_prime(80);
  is( miller_rabin_random( $randprime, 20 ), 1, "80-bit prime passes Miller-Rabin with 20 random bases" );
  do { $randprime += 2 } while is_prime($randprime);
  is( miller_rabin_random( $randprime, 40 ), "0", "80-bit composite fails Miller-Rabin with 40 random bases" );

  # Failure and shortcuts for MRR:
  ok(!eval { miller_rabin_random(undef,4); },   "MRR(undef,4)");
  ok(!eval { miller_rabin_random(10007,-4); },   "MRR(10007,-4)");
  # Note use of 1-1 : bigint on perl 5.6 and 5.8 is totally borked
  is(miller_rabin_random(10007, 1-1), 1, "MRR(n,0) = 1");
  is(miller_rabin_random(61, 17), 1, "MRR(61,17) = 1");
  is(miller_rabin_random(62, 17), 1-1, "MRR(62,17) = 0");
  is(miller_rabin_random(1009), 1, "MRR(1009) = 1");   # runs one random base

  } # ^^^ skipped without $extra
};

###############################################################################

subtest 'Perrin primality test', sub {
  SKIP: {
  skip "Perrin pseudoprime tests without EXTENDED_TESTING.", 2 unless $extra;
  my $perrinpsp = "1872702918368901354491086980308187833191468631072304770659547218657051750499825897279325406141660412842572655186363032039901203993254366727915836984799032960354882761038920216623610400227219443050113697104123375722324640843102690830473074828429679607154504449403902608511103291058038852618235905156930862492532896467422733403061010774542590301998535381232230279731082501";
  # It's fast with a *new* version of the GMP code (that has the test).
  is( is_perrin_pseudoprime($perrinpsp), 1, "18727...2501 is a Perrin PRP" );
  is( is_bpsw_prime($perrinpsp), 1-1, "18727...2501 is not a BPSW prime" );
  }
};

###############################################################################

subtest 'vecequal', sub {
  my $ten = Math::BigInt->new("10");
  my $six = Math::BigInt->new("6");
  is( vecequal([$ten,20],[$ten,20]), 1, "vecequal with Math::BigInt" );
  is( vecequal([$ten,20],[10,20]), 1, "vecequal with Math::BigInt and scalar" );
  is( vecequal([$ten,$six],[$ten,$six]), 1, "vecequal with equal Math::BigInt" );
  is( vecequal([$ten,20],[$six,20]), "0", "vecequal with unequal Math::BigInt" );

  ok(!eval { vecequal([$ten,{}],[$ten,{}]); }, "vecequal with hash should error");
};

###############################################################################

is( $_, 'this should not change', "Nobody clobbered \$_" );


sub check_pcbounds {
  my ($n, $expn, $percent, $percentrh) = @_;
  $percent   = Math::BigFloat->new($percent);
  $percentrh = Math::BigFloat->new($percentrh);

  my $pcap = prime_count_approx($n);
  is( "$pcap", "$expn", "PC approx($n)" );

  my $pclo = prime_count_lower($n);
  my $pcup = prime_count_upper($n);
  prime_set_config(assume_rh=>1);
  my $pclo_rh = prime_count_lower($n);
  my $pcup_rh = prime_count_upper($n);
  prime_set_config(assume_rh => undef);

  #diag "lower:    " . $pclo->bstr() . "  " . ($pcap-$pclo)->bstr;
  #diag "rh lower: " . $pclo_rh->bstr() . "  " . ($pcap-$pclo_rh)->bstr;
  #diag "approx:   " . $pcap->bstr();
  #diag "rh upper: " . $pcup_rh->bstr() . "  " . ($pcup_rh-$pcap)->bstr;
  #diag "upper:    " . $pcup->bstr() . "  " . ($pcup-$pcap)->bstr;
  # lower:    544534406675337676203117  17049931809477521306
  # rh lower: 544551456594152957592704  12994196131719
  # approx:   544551456607147153724423
  # rh upper: 544551456620339152603564  13191998879141
  # upper:    544586259732074697890498  34803124927544166075

  ok( $pclo <= $pclo_rh && $pclo_rh <= $pcap &&
      $pcap <= $pcup_rh && $pcup_rh <= $pcup,
      "prime count bounds for $n are in the right order");

  my $pcapf = Math::BigFloat->new($pcap);
#diag "" . ($pcapf - $pclo_rh)/($pcapf) . "  " . $percentrh/100.0 . "";
  cmp_ok( ($pcapf - $pclo_rh)/$pcapf, '<=', $percentrh , "PC lower with RH");
  cmp_ok( ($pcup_rh - $pcapf)/$pcapf, '<=', $percentrh , "PC upper with RH");
  cmp_ok( ($pcapf - $pclo)/$pcapf,    '<=', $percent   , "PC lower");
  cmp_ok( ($pcup - $pcapf)/$pcapf,    '<=', $percent   , "PC upper");
}

###############################################################################

sub linear_to_exp {   # Convert factor() output to factor_exp() output
  my %exponents;
  my @factors = grep { !$exponents{$_}++ } @_;
  return (map { [$_, $exponents{$_}] } @factors);
}

sub is_one_of {
  my($n, @list) = @_;
  if (defined $n) {
    for (@list) {
      return 1 if defined $_ && "$n" eq $_;
    }
  } else {
    for (@list) {
      return 1 if !defined $_;
    }
  }
  0;
}
