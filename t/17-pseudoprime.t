#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_prime
                         is_pseudoprime
                         is_euler_pseudoprime
                         is_euler_plumb_pseudoprime
                         is_strong_pseudoprime
                         is_lucas_pseudoprime
                         is_strong_lucas_pseudoprime
                         is_extra_strong_lucas_pseudoprime
                         is_almost_extra_strong_lucas_pseudoprime
                         is_frobenius_underwood_pseudoprime
                         is_frobenius_khashin_pseudoprime
                         is_perrin_pseudoprime
                         is_catalan_pseudoprime
                         is_frobenius_pseudoprime
                         lucasumod kronecker powmod submod modint
                         vecsample irand urandomb/;

my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
my $usegmp =Math::Prime::Util::prime_get_config->{'gmp'};
my $extra = 0+(defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING});

# small primes
my @sp = qw/2 3 5 7 11 13 17 19 23 29 31 37 41 43 47 53 59 61 67 71 73 79 83 89 97/;

plan tests => 1    # invalid inputs
            + 1    # simple
            + 1    # strong pseudoprimes
            + 1    # pseudoprimes
            + 1    # Lucas pseudoprimes
            + 1    # other pseudoprimes
            + 1    # Perrin pseudoprimes
            + 1    # Catalan pseudoprimes
            + 1;   # Frobenius pseudoprimes

subtest 'invalid inputs should croak', sub {
  ok(!eval { is_strong_pseudoprime(2047,0); }, "MR base 0 fails");
  ok(!eval { is_strong_pseudoprime(2047,1); }, "MR base 1 fails");
};

subtest 'basic functionality', sub {
  is( is_strong_pseudoprime(0, 2), 0, "MR with 0 shortcut composite");
  is( is_strong_pseudoprime(1, 2), 0, "MR with 0 shortcut composite");
  is( is_strong_pseudoprime(2, 2), 1, "MR with 2 shortcut prime");
  is( is_strong_pseudoprime(3, 2), 1, "MR with 3 shortcut prime");

  is_deeply([map{is_pseudoprime($_)} 162193,452051],[1,1],"is_pseudoprime(n) = is_pseudoprime(n,2)");
  my @b235 = (2,3,5);
  is_deeply([map{is_pseudoprime($_,2,3,5)} 18595801,22066201],
            [map{is_pseudoprime($_,@b235)} 18595801,22066201],"is_pseudoprime(n,\@baselist)");
  my @bnull = ();
  is_deeply([map{is_pseudoprime($_,2)     } 40165093,66437841],
            [map{is_pseudoprime($_,@bnull)} 40165093,66437841],"is_pseudoprime(n,())");
};


subtest 'pseudoprimes (Fermat test)', sub {
  my @psp = (
    [qw/2 341 561 645 1105 1387 1729 1905 2047 2465 2701 2821 3277 4033 4369 4371 4681 5461 6601 7957 8321 8481 8911 10261 10585 11305 12801 13741 13747 13981 14491 15709 15841 16705 18705 18721 19951 23001 23377 25761 29341/],
    [qw/3 91 121 286 671 703 949 1105 1541 1729 1891 2465 2665 2701 2821 3281 3367 3751 4961 5551 6601 7381 8401 8911 10585 11011 12403 14383 15203 15457 15841 16471 16531 18721 19345 23521 24046 24661 24727 28009 29161/],
  );
  for my $pdata (@psp) {
    my($base, @vals) = @$pdata;
    my @fails = grep { !is_pseudoprime($_,$base) } @vals;
    is_deeply(\@fails, [], "Small PSP-$base");
  }
  ok(is_pseudoprime(4398117272641,3),"Large PSP-3") if $use64;
  ok(is_pseudoprime(143168581, 2, 3, 5, 7, 11), "143168581 is a Fermat pseudoprime to bases 2,3,5,7,11");
};


subtest 'strong pseudoprimes (Miller-Rabin test)', sub {
  # Using a different codebase to get reference values:
  #   perl -E '$|=1; use Math::Primality ":all"; for (2 .. 1000000) { print "$_ " if $_&1 && is_strong_pseudoprime($_,17) && !is_prime($_); } print "\n"'
  #
  # Compare to:
  #   perl -E '$|=1; use ntheory ":all"; foroddcomposites { print "$_ " if is_strong_pseudoprime($_,17); } 1000000; print "\n";'

  # strong pseudoprimes from 2-100k for various bases (base is the first value)
  my @spsp = (
    [qw/2 2047 3277 4033 4681 8321 15841 29341 42799 49141 52633 65281 74665 80581 85489 88357 90751 1194649/ ],
    [qw/3 121 703 1891 3281 8401 8911 10585 12403 16531 18721 19345 23521 31621 44287 47197 55969 63139 74593 79003 82513 87913 88573 97567/ ],
    [qw/5 781 1541 5461 5611 7813 13021 14981 15751 24211 25351 29539 38081 40501 44801 53971 79381/ ],
    [qw/7 25 325 703 2101 2353 4525 11041 14089 20197 29857 29891 39331 49241 58825 64681 76627 78937 79381 87673 88399 88831/ ],
   [qw/11 133 793 2047 4577 5041 12403 13333 14521 17711 23377 43213 43739 47611 48283 49601 50737 50997 56057 58969 68137 74089 85879 86347 87913 88831/ ],
   [qw/13 85 1099 5149 7107 8911 9637 13019 14491 17803 19757 20881 22177 23521 26521 35371 44173 45629 54097 56033 57205 75241 83333 85285 86347/ ],
   [qw/17 9 91 145 781 1111 2821 4033 4187 5365 5833 6697 7171 15805 19729 21781 22791 24211 26245 31621 33001 33227 34441 35371 38081 42127 49771 71071 74665 77293 78881 88831 96433 97921 98671/ ],
   [qw/19 9 49 169 343 1849 2353 2701 4033 4681 6541 6697 7957 9997 12403 13213 13747 15251 16531 18769 19729 24761 30589 31621 31861 32477 41003 49771 63139 64681 65161 66421 68257 73555 96049/ ],
   [qw/23 169 265 553 1271 2701 4033 4371 4681 6533 6541 7957 8321 8651 8911 9805 14981 18721 25201 31861 34133 44173 47611 47783 50737 57401 62849 82513 96049/ ],
   [qw/29 15 91 341 469 871 2257 4371 4411 5149 6097 8401 11581 12431 15577 16471 19093 25681 28009 29539 31417 33001 48133 49141 54913 79003/ ],
   [qw/31 15 49 133 481 931 6241 8911 9131 10963 11041 14191 17767 29341 56033 58969 68251 79003 83333 87061 88183/ ],
   [qw/37 9 451 469 589 685 817 1333 3781 8905 9271 18631 19517 20591 25327 34237 45551 46981 47587 48133 59563 61337 68101 68251 73633 79381 79501 83333 84151 96727/ ],
   [qw/61 217 341 1261 2701 3661 6541 6697 7613 13213 16213 22177 23653 23959 31417 50117 61777 63139 67721 76301 77421 79381 80041/ ],
   [qw/73 205 259 533 1441 1921 2665 3439 5257 15457 23281 24617 26797 27787 28939 34219 39481 44671 45629 64681 67069 76429 79501 93521/ ],
   [qw/325 341 343 697 1141 2059 2149 3097 3537 4033 4681 4941 5833 6517 7987 8911 12403 12913 15043 16021 20017 22261 23221 24649 24929 31841 35371 38503 43213 44173 47197 50041 55909 56033 58969 59089 61337 65441 68823 72641 76793 78409 85879/ ],
   [qw/9375 11521 14689 17893 18361 20591 28093 32809 37969 44287 60701 70801 79957 88357 88831 94249 96247 99547/ ],
   [qw/28178 28179 29381 30353 34441 35371 37051 38503 43387 50557 51491 57553 79003 82801 83333 87249 88507 97921 99811/ ],
   [qw/75088 75089 79381 81317 91001 100101 111361 114211 136927 148289 169641 176661 191407 195649/ ],
   [qw/450775 465991 468931 485357 505441 536851 556421 578771 585631 586249 606361 631651 638731 641683 645679/ ],
   [qw/642735 653251 653333 663181 676651 714653 759277 794683 805141 844097 872191 874171 894671/ ],
   [qw/9780504 9780505 9784915 9826489 9882457 9974791 10017517 10018081 10084177 10188481 10247357 10267951 10392241 10427209 10511201/ ],
   [qw/203659041 204172939 204456793 206407057 206976001 207373483 209301121 210339397 211867969 212146507 212337217 212355793 214400629 214539841 215161459/ ],
   [qw/553174392 553174393 553945231 554494951 554892787 555429169 557058133 557163157 557165209 558966793 559407061 560291719 561008251 563947141/ ],
   [qw/1005905886 1005905887 1007713171 1008793699 1010415421 1010487061 1010836369 1012732873 1015269391 1016250247 1018405741 1020182041/ ],
   [qw/1340600841 1345289261 1345582981 1347743101 1348964401 1350371821 1353332417 1355646961 1357500901 1361675929 1364378203 1366346521 1367104639/ ],
   [qw/1795265022 1795265023 1797174457 1797741901 1804469753 1807751977 1808043283 1808205701 1813675681 1816462201 1817936371 1819050257/ ],
   [qw/3046413974 3046413975 3048698683 3051199817 3068572849 3069705673 3070556233 3079010071 3089940811 3090723901 3109299161 3110951251 3113625601/ ],
   [qw/3613982119 3626488471 3630467017 3643480501 3651840727 3653628247 3654142177 3672033223 3672036061 3675774019 3687246109 3690036017 3720856369/ ],
  );
  for my $pdata (@spsp) {
    my($base, @vals) = @$pdata;
    my @fails = grep { !is_strong_pseudoprime($_,$base) } @vals;
    is_deeply(\@fails, [], "Small SPSP-$base");
  }

  ok( is_strong_pseudoprime(1099558795087,3),"Large SPSP-3") if $use64;

  # Verify MR for bases >= n
  is( is_strong_pseudoprime(  3,    3), 1, "spsp(  3,    3)");
  is( is_strong_pseudoprime( 11,   11), 1, "spsp( 11,   11)");
  is( is_strong_pseudoprime( 89, 5785), 1, "spsp( 89, 5785)");
  is( is_strong_pseudoprime(257, 6168), 1, "spsp(257, 6168)");
  is( is_strong_pseudoprime(367,  367), 1, "spsp(367,  367)");
  is( is_strong_pseudoprime(367, 1101), 1, "spsp(367, 1101)");
  is( is_strong_pseudoprime(49001, 921211727), 0, "spsp(49001, 921211727)");
  is( is_strong_pseudoprime(  331, 921211727), 1, "spsp(  331, 921211727)");
  is( is_strong_pseudoprime(49117, 921211727), 1, "spsp(49117, 921211727)");

  { # Verify MR base 2 for all small numbers
    my @got = map { is_strong_pseudoprime($_) } 2..4032;
    my @exp = map { ($_==2047 || $_==3277 || is_prime($_)) ? 1 : 0 } 2..4032;
    is_deeply(\@got, \@exp, "MR base 2 for 2..4032");
  }

  SKIP: {
    skip "base 2,3 without EXTENDED_TESTING",1 unless $extra;
    my @tnum = vecsample(10000, 1..1373652);
    my $fails = 0;
    for (@tnum) {
      $fails++ if !!is_strong_pseudoprime($_,2,3) != !!is_prime($_);
    }
    is($fails, 0, "is_strong_pseudoprime bases 2,3 matches is_prime");
  }

  { # strong pseudoprimes for all prime bases 2 .. pn
    my @phis = qw/2047 1373653 25326001 3215031751 2152302898747 3474749660383 341550071728321 341550071728321/;
    $#phis = 3 unless $use64;
    my @got = map { is_strong_pseudoprime($phis[$_], @sp[0 .. $_]) } 0..$#phis;
    my @exp = map { 1 } 0 .. $#phis;
    is_deeply(\@got, \@exp, "A014233: first strong pseudoprime to N prime bases");
  }

  # A 77-bit composite that passes many bases
  ok(is_strong_pseudoprime("318665857834031151167461", 2,3,5,7,11,13,17,19,23,29,31,37),"318665857834031151167461 is a pseudoprime to many bases");

  # I don't think we benefit from doing these tests.
  # '21652684502221' => [ qw/2 7 37 61 9375/ ],
  # '3825123056546413051' => [ qw/2 3 5 7 11 13 17 19 23 29 31 325 9375/ ],
  # '318665857834031151167461' => [ qw/2 3 5 7 11 13 17 19 23 29 31 37 325 9375/ ],
  # '3317044064679887385961981' => [ qw/2 3 5 7 11 13 17 19 23 29 31 37 73 325 9375/ ],
  # '6003094289670105800312596501' => [ qw/2 3 5 7 11 13 17 19 23 29 31 37 61 325 9375/ ],
  # '59276361075595573263446330101' => [ qw/2 3 5 7 11 13 17 19 23 29 31 37 325 9375/ ],
  # '564132928021909221014087501701' => [ qw/2 3 5 7 11 13 17 19 23 29 31 37 325 9375/ ],
  # '1543267864443420616877677640751301' => [ qw/2 3 5 7 11 13 17 19 23 29 31 37 61 325 9375/ ],
};


subtest 'Lucas pseudoprimes', sub {
  my @P;

  @P = qw/323 377 1159 1829 3827 5459 5777 9071 9179 10877 11419 11663 13919 14839 16109 16211 18407 18971 19043/;
  is_deeply([grep{ !is_lucas_pseudoprime($_)} @P],[],"Small Lucas");

  @P = qw/5459 5777 10877 16109 18971 22499 24569 25199 40309 58519 75077 97439 100127 113573 115639 130139/;
  is_deeply([grep{ !is_strong_lucas_pseudoprime($_)} @P],[],"Small strong Lucas");

  @P = qw/989 3239 5777 10877 27971 29681 30739 31631 39059 72389 73919 75077 100127 113573 125249 137549 137801 153931 155819/;
  is_deeply([grep{ !is_extra_strong_lucas_pseudoprime($_)} @P],[],"Small extra strong Lucas");

  @P = qw/989 3239 5777 10469 10877 27971 29681 30739 31631 39059 72389 73919 75077 100127 113573 125249 137549 137801 153931 154697 155819/;
  is_deeply([grep{ !is_almost_extra_strong_lucas_pseudoprime($_)} @P],[],"Small almost extra strong Lucas");

  @P = qw/3239 4531 5777 10877 12209 21899 31631 31831 32129 34481 36079 37949 47849 50959 51641 62479 73919 75077 97109 100127 108679 113573 116899 154697 161027/;
  is_deeply([grep{ !is_almost_extra_strong_lucas_pseudoprime($_,2)} @P],[],"Small almost extra strong Lucas (increment 2)");

  # Verify extra strong for a few small primes
  is_deeply( [grep { is_extra_strong_lucas_pseudoprime($_) } 2..100], [grep { $_ >= 2 && $_ <= 100 } @sp], "The first 100 primes are selected by is_extra_strong_lucas_pseudoprime" );

  # Verify Lucas for some small numbers
  ok(is_strong_lucas_pseudoprime(2),"is_strong_lucas_pseudoprime(2) = 1");
  my @C = (9, 16, 100, 102, 2047, 2048, 5781, 9000, 14381);
  is_deeply([grep{is_lucas_pseudoprime($_)}@C],[],"Not SLPSP: [@C]");

  if ($use64) {
    my @LPSP    = (2199055761527);
    my @SLPSP   = (4294967311,4294967357,12598021314449);
    my @ESLPSP  = (4294967311,4294967357,10099386070337);
    my @AESLPSP = (4294967311,4294967357,10071551814917);

    is_deeply([map {is_lucas_pseudoprime($_)} @LPSP], [map{1}@LPSP], "LPSP: @LPSP");
    is_deeply([map {is_strong_lucas_pseudoprime($_)} @SLPSP], [map{1}@SLPSP], "SLPSP: @SLPSP");
    is_deeply([map {is_extra_strong_lucas_pseudoprime($_)} @ESLPSP], [map{1}@ESLPSP], "ESLPSP: @ESLPSP");
    is_deeply([map {is_almost_extra_strong_lucas_pseudoprime($_)} @AESLPSP], [map{1}@AESLPSP], "AESLPSP: @AESLPSP");
  }
};


subtest 'other pseudoprimes', sub {
  my @P;

  @P = qw/1729 1905 2047 2465 3277 4033 4681 8321 12801 15841 16705 18705 25761 29341 33153 34945 41041 42799 46657 49141 52633 65281 74665 75361 80581 85489 87249 88357 90751/;
  is_deeply([grep{ !is_euler_plumb_pseudoprime($_)} @P],[],"Small Euler-Plumb");

  @P = qw/561 1105 1729 1905 2047 2465 3277 4033 4681 6601 8321 8481 10585 12801 15841 16705 18705 25761 29341 30121 33153 34945 41041 42799 46657 49141 52633 62745 65281 74665 75361 80581 85489 87249 88357 90751/;
  is_deeply([grep{ !is_euler_pseudoprime($_)} @P],[],"Small Euler base 2");

  @P = qw/121 703 1729 1891 2821 3281 7381 8401 8911 10585 12403 15457 15841 16531 18721 19345 23521 24661 28009 29341 31621 41041 44287 46657 47197 49141 50881 52633 55969 63139 63973 74593 75361 79003 82513 87913 88573 93961 97567/;
  is_deeply([grep{ !is_euler_pseudoprime($_,3)} @P],[],"Small Euler base 2");

  @P = qw/15 91 341 469 871 2257 4371 4411 5149 5185 6097 8401 8841 11581 12431 15577 15841 16471 19093 22281 25681 27613 28009 29539 31417 33001 41041 46657 48133 49141 54913 57889 79003 98301/;
  is_deeply([grep{ !is_euler_pseudoprime($_,29)} @P],[],"Small Euler base 29");

  # Following are not used for anything, but interesting

  # Euler-Jacobi (A047713)
  @P = qw/561 1105 1729 1905 2047 2465 3277 4033 4681 6601 8321 8481 10585 12801 15841 16705 18705 25761 29341 30121 33153 34945 41041 42799 46657 49141 52633 62745 65281 74665 75361 80581 85489 87249 88357 90751/;
  is_deeply([grep{ powmod(2,$_>>1,$_) != modint(kronecker(2,$_),$_) } @P],[],"Small Euler-Jacobi base 2");

  # Fibonacci (A081264, Wikipedia first definition)
  @P = qw/323 377 1891 3827 4181 5777 6601 6721 8149 10877 11663 13201 13981 15251 17119 17711 18407 19043 23407 25877 27323/;
  my @fmod5 = (0,-1,1,1,-1);
  is_deeply([grep{lucasumod(1,-1,$_+$fmod5[$_%5],$_)} @P],[],"Small Fibonacci");

  @P = qw/169 385 741 961 1121 2001 3827 4879 5719 6215 6265 6441 6479 6601 7055 7801 8119 9799 10945 11395 13067 13079 13601 15841 18241 19097 20833 20951 24727 27839 27971 29183 29953/;
  is_deeply([grep { submod(lucasumod(2,-1,$_,$_),kronecker(2,$_),$_) } @P],[],"Small Pell");
};


subtest 'Perrin pseudoprimes', sub {
  my @P = qw/271441 904631 16532714 24658561 27422714 27664033 46672291/;
  is_deeply([grep{ !is_perrin_pseudoprime($_)} @P],[],"Small Perrin");

  SKIP: {
    skip "larger Perrin",1 unless $extra;
    my @P = qw/102690901 130944133 196075949 214038533 517697641 545670533 801123451/;
    is_deeply([grep{ !is_perrin_pseudoprime($_)} @P],[],"Larger Perrin");
  }

  # Perrin restrictions
  is( is_perrin_pseudoprime(271441,0), 1, "271441 is an unrestricted Perrin pseudoprime");
  is( is_perrin_pseudoprime(271441,1), 0, "271441 is not a minimal restricted Perrin pseudoprime");
  SKIP: {
    skip "restrictions without EXTENDED_TESTING on 32-bit",4 unless $use64 || $extra;
    is( is_perrin_pseudoprime("36407440637569",1), 1, "36407440637569 is minimal restricted Perrin pseudoprime");
    is( is_perrin_pseudoprime("36407440637569",2), 0, "36407440637569 is not an Adams/Shanks Perrin pseudoprime");
    is( is_perrin_pseudoprime("364573433665",2), 1, "364573433665 is an Adams/Shanks Perrin pseudoprime");
    is( is_perrin_pseudoprime("364573433665",3), 0, "364573433665 is not a Grantham restricted Perrin pseudoprime");
  }

  # Large Perrin pseudoprime.  Very slow without GMP.
  SKIP: {
    skip "very large pseudoprime without EXTENDED_TESTING",1 unless $extra;
    skip "very large pseudoprime without GMP backend",1 unless $usegmp && $Math::Prime::Util::GMP::VERSION >= 0.40;
    my $perrinpsp = "1872702918368901354491086980308187833191468631072304770659547218657051750499825897279325406141660412842572655186363032039901203993254366727915836984799032960354882761038920216623610400227219443050113697104123375722324640843102690830473074828429679607154504449403902608511103291058038852618235905156930862492532896467422733403061010774542590301998535381232230279731082501";
    is( is_perrin_pseudoprime($perrinpsp), 1, "18727...2501 is a Perrin PRP" );
  }
};

subtest 'Catalan pseudoprimes', sub {
  my @P = (5907);

  # Without XS or GMP, these are soooo slow.
  push @P, 1194649, 12327121 if $usexs || ($extra && $usegmp && defined &Math::Prime::Util::GMP::binomial && $Math::Prime::Util::GMP::VERSION >= 0.27);

  is_deeply([grep{ !is_catalan_pseudoprime($_)} @P],[],"Catalan [@P]");
};

subtest 'Frobenius type pseudoprimes', sub {
  my @P;

  @P = qw/4181 5777 6721 10877 13201 15251 34561 51841 64079 64681 67861 68251 75077 90061 96049 97921 100127/;
  is_deeply([grep{ !is_frobenius_pseudoprime($_,1,-1)} @P],[],"Small Frobenius(1,-1)");

  @P = qw/13333 44801 486157 1615681 3125281 4219129 9006401 12589081 13404751 15576571 16719781/;
  is_deeply([grep{ !is_frobenius_pseudoprime($_,3,-5)} @P],[],"Small Frobenius(3,-5)");

  # These have no known counterexamples
  {
    my $ntests = $usexs ? 100 : 2;
    @P = map { irand() & 1 } 1 .. $ntests;   # Random odd
    @P = map { $_ % 3 ? $_ : $_-2 } @P;      # not divisible by 3
    @P = map { $_ % 5 ? $_ : $_+6 } @P;      # not divisible by 5
    my @ISP = map { [$_,!!is_prime($_)] } @P;
    is_deeply([map { [$_,!!is_frobenius_underwood_pseudoprime($_)] } @P],\@ISP,"32-bit Frobenius Underwood ($ntests random)");
    is_deeply([map { [$_,!!is_frobenius_khashin_pseudoprime($_)] } @P],\@ISP,"32-bit Frobenius Khashin ($ntests random)");
  }
  if ($use64) {
    my $ntests = 2;
    @P = map { urandomb(64) & 1 } 1 .. $ntests;
    @P = map { $_ % 3 ? $_ : $_-2 } @P;      # not divisible by 3
    @P = map { $_ % 5 ? $_ : $_+6 } @P;      # not divisible by 5
    my @ISP = map { [$_,!!is_prime($_)] } @P;
    is_deeply([map { [$_,!!is_frobenius_underwood_pseudoprime($_)] } @P],\@ISP,"64-bit Frobenius Underwood ($ntests random)");
    is_deeply([map { [$_,!!is_frobenius_khashin_pseudoprime($_)] } @P],\@ISP,"64-bit Frobenius Khashin ($ntests random)");
  }
};
