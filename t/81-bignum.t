#!/usr/bin/env perl
use strict;
use warnings;

# If you're not using ancient perl 5.6.2 with super early releases of bigint,
# then you can define bigint up here and not have to quote every number.

# Note: In 5.16.0 (and perhaps others?), using labels like "SKIP:" will create
# a small memory leak.  So running the test suite through valgrind will show
# some small leaks in this test, which has nothing to do with the module.

my $extra = defined $ENV{RELEASE_TESTING} && $ENV{RELEASE_TESTING};
my $broken64 = (18446744073709550592 == ~0);

use Test::More;

my @primes = qw/
777777777777777777777767 777777777777777777777787
877777777777777777777753 877777777777777777777871
87777777777777777777777795577 890745785790123461234805903467891234681243
618970019642690137449562111
/;
push @primes, "531137992816767098689588206552468627329593117727031923199444138200403559860852242739162502265229285668889329486246501015346579337652707239409519978766587351943831270835393219031728127"  if $extra;

my @composites = qw/
777777777777777777777777 877777777777777777777777
87777777777777777777777795475 890745785790123461234805903467891234681234
/;

# pseudoprimes to various small prime bases
# These must not be themselves prime, as we're going to primality test them.
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
my $num_pseudoprime_tests = 0;
foreach my $psrp (keys %pseudoprimes) {
  push @composites, $psrp;
  $num_pseudoprime_tests += scalar @{$pseudoprimes{$psrp}};
}

my %factors = (
  '1234567890' => [2, 3, 3, 5, 3607, 3803],
  '190128090927491' => [61, 73, 196291, 217517],
  '23489223467134234890234680' => [2, 2, 2, 5, 4073, 4283, "33662485846146713"],
  #'7674353466844111807691499613711' => [11783, 12239, 18869, 22277, 37861, 55163, 60617],
);

my %allfactors = (
  #'7674353466844111807691499613711' => [qw/11783 12239 18869 22277 37861 55163 60617 144212137 222333427 230937691 262489891 272648203 420344713 446116163 463380779 649985629 675139957 714250111 714399209 741891463 843429497 1040870647 1143782173 1228866151 1350364909 2088526343 2295020237 3343815571 2721138813053 3212613775949 4952921753279 5144598942407 5460015718957 7955174113331 8417765879647 8741707108529 8743531918951 9938129763151 10322733613783 12264578833601 12739215848633 13477185344459 13998750015347 14479729857233 15040092822089 15911349722747 15914671178893 16527116121251 23187475403219 24609105899569 25480035467921 25561473911977 27042223452571 28088752680643 39400178873093 39408403566067 40924958773469 43304736851953 46526101343011 51126165819649 63094456009199 74490179475167 126600201333631 60618809338381681 103025036600999633 121632770171205089 150106180344442639 164947271430833701 177217413722674687 187522570500896219 194739009256700533 194779660558471427 273218022676129477 283791511459997341 300231257918513143 301190847104824991 311850154091885119 330969772836016469 464349219218967461 482218789227785227 482319451245094013 510259714326562199 530006674331052767 548217052124698613 569432954337111629 602419611852924167 625733143466684111 743441975156391817 772213047096586561 877717784755892761 877901006241274559 911685306596568913 964699622850956981 1405555196516926123 1491730172314174073 1549459864122309809 2388819198964283339 2820272685109297787 2295088740360468824341 3343915379533148669003 3674530365664682357177 5683170094020942755179 6245068643642794753561 6709628500954186324507 7373013629467938879913 9098986333939079448263 10344307556540938128697 10742387967627371501879 10744630415386959327601 11367055656052826107123 11806958684072862490459 16561656880558940507309 17202590050170658819397 18257285578953176479447 28147456621396150583437 29236758176123863786021 33231273048642855824221 34517317393052695615093 53215725295327339942903 126603980184504541757122583 139121394174430538725078397 202698118561160872868954851 344496721589267486990685443 406717550842339912432640819 627040891154842046547226049 651307261889511313561189817/],
  '23489223467134234890234680' => [qw/2 4 5 8 10 20 40 4073 4283 8146 8566 16292 17132 20365 21415 32584 34264 40730 42830 81460 85660 162920 171320 17444659 34889318 69778636 87223295 139557272 174446590 348893180 697786360 33662485846146713 67324971692293426 134649943384586852 168312429230733565 269299886769173704 336624858461467130 673249716922934260 1346499433845868520 137107304851355562049 144176426879046371779 274214609702711124098 288352853758092743558 548429219405422248196 576705707516185487116 685536524256777810245 720882134395231858895 1096858438810844496392 1153411415032370974232 1371073048513555620490 1441764268790463717790 2742146097027111240980 2883528537580927435580 5484292194054222481960 5767057075161854871160 587230586678355872255867 1174461173356711744511734 2348922346713423489023468 2936152933391779361279335 4697844693426846978046936 5872305866783558722558670 11744611733567117445117340/],
);

plan tests =>  0
             + 2*(@primes + @composites)
             + 1   # primes
             + 2   # next/prev prime
             + 1   # primecount large base small range
             + scalar(keys %pseudoprimes)
             + 6   # PC lower, upper, approx
             + 6*2*$extra # more PC tests
             + scalar(keys %factors)
             + scalar(keys %allfactors)
             + 2   # moebius, euler_phi
             + 12  # random primes
             + 0;

# Using GMP makes these tests run about 2x faster on some machines
use bigint try => 'GMP';   #  <--------------- large numbers ahead!  > 2^64
use Math::BigFloat;

use Math::Prime::Util qw/
  prime_set_config
  is_prob_prime
  prime_count_lower
  prime_count_upper
  prime_count_approx
  nth_prime_lower
  nth_prime_upper
  nth_prime_approx
  factor
  all_factors
  moebius
  euler_phi
  ExponentialIntegral
  LogarithmicIntegral
  RiemannR
  primes
  prime_count
  nth_prime
  is_prime
  next_prime
  prev_prime
  is_strong_pseudoprime
  random_prime
  random_ndigit_prime
  random_nbit_prime
  random_maurer_prime
/;
# TODO:  is_strong_lucas_pseudoprime
#        ExponentialIntegral
#        LogarithmicIntegral
#        RiemannR

# See if we're testing bignum with PP or GMP.
diag "Math::Prime::Util::GMP not being used.\n"
  unless Math::Prime::Util::prime_get_config->{gmp};


###############################################################################

foreach my $n (@primes) {
  ok( is_prime($n), "$n is prime" );
  ok( is_prob_prime($n), "$n is probably prime");
}
foreach my $n (@composites) {
  ok( !is_prime($n), "$n is not prime" );
  ok( !is_prob_prime($n), "$n is not probably prime");
}

###############################################################################

is_deeply( primes(2**66, 2**66+100), [73786976294838206473,73786976294838206549], "primes( 2^66, 2^66 + 100 )" );

###############################################################################

is( next_prime(777777777777777777777777), 777777777777777777777787, "next_prime(777777777777777777777777)");
is( prev_prime(777777777777777777777777), 777777777777777777777767, "prev_prime(777777777777777777777777)");

###############################################################################

# Testing prime_count only on a small range -- it would take a very long time
# otherwise.

is( prime_count(877777777777777777777752, 877777777777777777777872), 2, "prime_count(87..7752, 87..7872)");

###############################################################################

# Testing nth_prime would be far too time consuming.

###############################################################################

while (my($psrp, $baseref) = each (%pseudoprimes)) {
  SKIP: {
    skip "Your 64-bit Perl is broken, skipping pseudoprime tests for $psrp", 1 if $broken64 && $psrp == 3825123056546413051;
    my $baselist = join(",", @$baseref);
    my @expmr = map { 1 } @$baseref;
    my @gotmr = map { is_strong_pseudoprime($psrp, $_) } @$baseref;
    is_deeply(\@gotmr, \@expmr, "$psrp is a strong pseudoprime to bases $baselist");
  }
}

###############################################################################

check_pcbounds(31415926535897932384, 716115441142294636, '8e-5', '2e-8');
if ($extra) {
  check_pcbounds(314159265358979323846, 6803848951392700268, '7e-5', '5e-9');
  check_pcbounds(31415926535897932384626433, 544551456607147153724423, '4e-5', '3e-11');
}

###############################################################################

SKIP: {
  skip "Your 64-bit Perl is broken, skipping bignum factoring tests", scalar(keys %factors) + scalar(keys %allfactors) if $broken64;
  while (my($n, $factors) = each(%factors)) {
    is_deeply( [factor($n)], $factors, "factor($n)" );
  }
  while (my($n, $allfactors) = each(%allfactors)) {
    is_deeply( [all_factors($n)], $allfactors, "all_factors($n)" );
  }
}

###############################################################################

SKIP: {
  skip "Your 64-bit Perl is broken, skipping moebius and euler_phi tests", 2 if $broken64;
  is( moebius(618970019642690137449562110), -1, "moebius(618970019642690137449562110)" );
  is( euler_phi(618970019642690137449562110), 145857122964987051805507584, "euler_phi(618970019642690137449562110)" );
}

###############################################################################

my $randprime;

$randprime = random_prime(147573952590750158861, 340282366920939067930896100764782952647);
cmp_ok( $randprime, '>=', 147573952590750158861, "random range prime isn't too small");
cmp_ok( $randprime, '<=', 340282366920939067930896100764782952647, "random range prime isn't too big");
ok( is_prime($randprime), "random range prime is prime");

$randprime = random_ndigit_prime(25);
cmp_ok( $randprime, '>', 10**24, "random 25-digit prime isn't too small");
cmp_ok( $randprime, '<', 10**25, "random 25-digit prime isn't too big");
ok( is_prime($randprime), "random 25-digit prime is prime");

$randprime = random_nbit_prime(80);
cmp_ok( $randprime, '>', 2**79, "random 80-bit prime isn't too small");
cmp_ok( $randprime, '<', 2**80, "random 80-bit prime isn't too big");
ok( is_prime($randprime), "random 80-bit prime is prime");

SKIP: {
  skip "Your 64-bit Perl is broken, skipping maurer prime", 3 if $broken64;
  $randprime = random_maurer_prime(80);
  cmp_ok( $randprime, '>', 2**79, "random 80-bit Maurer prime isn't too small");
  cmp_ok( $randprime, '<', 2**80, "random 80-bit Maurer prime isn't too big");
  ok( is_prime($randprime), "random 80-bit Maurer prime is prime");
}

###############################################################################


sub check_pcbounds {
  my ($n, $expn, $percent, $percentrh) = @_;
  $percent   = Math::BigFloat->new($percent);
  $percentrh = Math::BigFloat->new($percentrh);

  my $pcap = prime_count_approx($n);
  is( $pcap, "$expn", "PC approx($n)" );

  my $pclo = prime_count_lower($n);
  my $pcup = prime_count_upper($n);
  prime_set_config(assume_rh=>1);
  my $pclo_rh = prime_count_lower($n);
  my $pcup_rh = prime_count_upper($n);
  prime_set_config(assume_rh=>0);

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
