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
                         lucasumod kronecker/;

my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
my $usegmp =Math::Prime::Util::prime_get_config->{'gmp'};
my $extra = 0+(defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING});

# small primes
my @sp = qw/2 3 5 7 11 13 17 19 23 29 31 37 41 43 47 53 59 61 67 71 73 79 83 89 97/;
# strong pseudoprimes for all prime bases 2 .. pn
my @phis = qw/2047 1373653 25326001 3215031751 2152302898747 3474749660383 341550071728321 341550071728321/;
$#phis = 3 unless $use64;

# pseudoprimes from 2-100k for each prime base
#
# Using a different codebase to get reference values:
#   perl -E 'use Math::Primality ":all"; for (2 .. 100000) { print "$_ " if is_strong_pseudoprime($_,17) && !is_prime($_); } print "\n"'
#
# With large values, one of:
#
#  perl -MMath::Primality=:all -E 'my $_=$base|1; while(1) {print "$_ " if is_strong_pseudoprime($_,$base) && !is_prime($_); $_+=2; } print "\n"; BEGIN {$|=1; $base=553174392}'
#
#  perl -MMath::Primality=is_strong_pseudoprime -MMath::Prime::Util=is_prime -E 'my $_=$base|1; while(1) {print "$_ " if is_strong_pseudoprime($_,$base) && !is_prime($_); $_+=2; } print "\n"; BEGIN {$|=1; $base=553174392}'
#
# ~30x faster than Math::Primality:
#  perl -MMath::Prime::Util=:all -E 'my $_=$base|1; while(1) {print "$_ " if is_strong_pseudoprime($_,$base) && !is_prime($_); $_+=2; } print "\n"; BEGIN {$|=1; $base=553174392}'

my %pseudoprimes = (
   2 => [ qw/2047 3277 4033 4681 8321 15841 29341 42799 49141 52633 65281 74665 80581 85489 88357 90751 1194649/ ],
   3 => [ qw/121 703 1891 3281 8401 8911 10585 12403 16531 18721 19345 23521 31621 44287 47197 55969 63139 74593 79003 82513 87913 88573 97567/ ],
   5 => [ qw/781 1541 5461 5611 7813 13021 14981 15751 24211 25351 29539 38081 40501 44801 53971 79381/ ],
   7 => [ qw/25 325 703 2101 2353 4525 11041 14089 20197 29857 29891 39331 49241 58825 64681 76627 78937 79381 87673 88399 88831/ ],
  11 => [ qw/133 793 2047 4577 5041 12403 13333 14521 17711 23377 43213 43739 47611 48283 49601 50737 50997 56057 58969 68137 74089 85879 86347 87913 88831/ ],
  13 => [ qw/85 1099 5149 7107 8911 9637 13019 14491 17803 19757 20881 22177 23521 26521 35371 44173 45629 54097 56033 57205 75241 83333 85285 86347/ ],
  17 => [ qw/9 91 145 781 1111 2821 4033 4187 5365 5833 6697 7171 15805 19729 21781 22791 24211 26245 31621 33001 33227 34441 35371 38081 42127 49771 71071 74665 77293 78881 88831 96433 97921 98671/ ],
  19 => [ qw/9 49 169 343 1849 2353 2701 4033 4681 6541 6697 7957 9997 12403 13213 13747 15251 16531 18769 19729 24761 30589 31621 31861 32477 41003 49771 63139 64681 65161 66421 68257 73555 96049/ ],
  23 => [ qw/169 265 553 1271 2701 4033 4371 4681 6533 6541 7957 8321 8651 8911 9805 14981 18721 25201 31861 34133 44173 47611 47783 50737 57401 62849 82513 96049/ ],
  29 => [ qw/15 91 341 469 871 2257 4371 4411 5149 6097 8401 11581 12431 15577 16471 19093 25681 28009 29539 31417 33001 48133 49141 54913 79003/ ],
  31 => [ qw/15 49 133 481 931 6241 8911 9131 10963 11041 14191 17767 29341 56033 58969 68251 79003 83333 87061 88183/ ],
  37 => [ qw/9 451 469 589 685 817 1333 3781 8905 9271 18631 19517 20591 25327 34237 45551 46981 47587 48133 59563 61337 68101 68251 73633 79381 79501 83333 84151 96727/ ],
         61 => [ qw/217 341 1261 2701 3661 6541 6697 7613 13213 16213 22177 23653 23959 31417 50117 61777 63139 67721 76301 77421 79381 80041/ ],
         73 => [ qw/205 259 533 1441 1921 2665 3439 5257 15457 23281 24617 26797 27787 28939 34219 39481 44671 45629 64681 67069 76429 79501 93521/ ],
        325 => [ qw/341 343 697 1141 2059 2149 3097 3537 4033 4681 4941 5833 6517 7987 8911 12403 12913 15043 16021 20017 22261 23221 24649 24929 31841 35371 38503 43213 44173 47197 50041 55909 56033 58969 59089 61337 65441 68823 72641 76793 78409 85879/ ],
       9375 => [ qw/11521 14689 17893 18361 20591 28093 32809 37969 44287 60701 70801 79957 88357 88831 94249 96247 99547/ ],
      28178 => [ qw/28179 29381 30353 34441 35371 37051 38503 43387 50557 51491 57553 79003 82801 83333 87249 88507 97921 99811/ ],
      75088 => [ qw/75089 79381 81317 91001 100101 111361 114211 136927 148289 169641 176661 191407 195649/ ],
     450775 => [ qw/465991 468931 485357 505441 536851 556421 578771 585631 586249 606361 631651 638731 641683 645679/ ],
     642735 => [ qw/653251 653333 663181 676651 714653 759277 794683 805141 844097 872191 874171 894671/ ],
    9780504 => [ qw/9780505 9784915 9826489 9882457 9974791 10017517 10018081 10084177 10188481 10247357 10267951 10392241 10427209 10511201/ ],
  203659041 => [ qw/204172939 204456793 206407057 206976001 207373483 209301121 210339397 211867969 212146507 212337217 212355793 214400629 214539841 215161459/ ],
  553174392 => [ qw/553174393 553945231 554494951 554892787 555429169 557058133 557163157 557165209 558966793 559407061 560291719 561008251 563947141/ ],
 1005905886 => [ qw/1005905887 1007713171 1008793699 1010415421 1010487061 1010836369 1012732873 1015269391 1016250247 1018405741 1020182041/ ],
 1340600841 => [ qw/1345289261 1345582981 1347743101 1348964401 1350371821 1353332417 1355646961 1357500901 1361675929 1364378203 1366346521 1367104639/ ],
 1795265022 => [ qw/1795265023 1797174457 1797741901 1804469753 1807751977 1808043283 1808205701 1813675681 1816462201 1817936371 1819050257/ ],
 3046413974 => [ qw/3046413975 3048698683 3051199817 3068572849 3069705673 3070556233 3079010071 3089940811 3090723901 3109299161 3110951251 3113625601/ ],
 3613982119 => [ qw/3626488471 3630467017 3643480501 3651840727 3653628247 3654142177 3672033223 3672036061 3675774019 3687246109 3690036017 3720856369/ ],
 psp2       => [ qw/341 561 645 1105 1387 1729 1905 2047 2465 2701 2821 3277 4033 4369 4371 4681 5461 6601 7957 8321 8481 8911 10261 10585 11305 12801 13741 13747 13981 14491 15709 15841 16705 18705 18721 19951 23001 23377 25761 29341/ ],
 psp3       => [ qw/91 121 286 671 703 949 1105 1541 1729 1891 2465 2665 2701 2821 3281 3367 3751 4961 5551 6601 7381 8401 8911 10585 11011 12403 14383 15203 15457 15841 16471 16531 18721 19345 23521 24046 24661 24727 28009 29161/ ],
 plumb      => [ qw/1729 1905 2047 2465 3277 4033 4681 8321 12801 15841 16705 18705 25761 29341 33153 34945 41041 42799 46657 49141 52633 65281 74665 75361 80581 85489 87249 88357 90751/ ],
 epsp2      => [ qw/561 1105 1729 1905 2047 2465 3277 4033 4681 6601 8321 8481 10585 12801 15841 16705 18705 25761 29341 30121 33153 34945 41041 42799 46657 49141 52633 62745 65281 74665 75361 80581 85489 87249 88357 90751/ ],
 epsp3      => [ qw/121 703 1729 1891 2821 3281 7381 8401 8911 10585 12403 15457 15841 16531 18721 19345 23521 24661 28009 29341 31621 41041 44287 46657 47197 49141 50881 52633 55969 63139 63973 74593 75361 79003 82513 87913 88573 93961 97567/ ],
 epsp29     => [ qw/15 91 341 469 871 2257 4371 4411 5149 5185 6097 8401 8841 11581 12431 15577 15841 16471 19093 22281 25681 27613 28009 29539 31417 33001 41041 46657 48133 49141 54913 57889 79003 98301/ ],
 lucas      => [ qw/323 377 1159 1829 3827 5459 5777 9071 9179 10877 11419 11663 13919 14839 16109 16211 18407 18971 19043/ ],
 slucas     => [ qw/5459 5777 10877 16109 18971 22499 24569 25199 40309 58519 75077 97439 100127 113573 115639 130139/ ],
 eslucas    => [ qw/989 3239 5777 10877 27971 29681 30739 31631 39059 72389 73919 75077 100127 113573 125249 137549 137801 153931 155819/ ],
 aeslucas1  => [ qw/989 3239 5777 10469 10877 27971 29681 30739 31631 39059 72389 73919 75077 100127 113573 125249 137549 137801 153931 154697 155819/ ],
 aeslucas2  => [ qw/3239 4531 5777 10877 12209 21899 31631 31831 32129 34481 36079 37949 47849 50959 51641 62479 73919 75077 97109 100127 108679 113573 116899 154697 161027/ ],
 perrin     => [ qw/271441 904631 16532714 24658561 27422714 27664033 46672291 102690901 130944133 196075949 214038533 517697641 545670533 801123451/ ],
 catalan    => [ qw/5907 1194649 12327121/ ],
 fibonacci  => [ qw/323 377 1891 3827 4181 5777 6601 6721 8149 10877 11663 13201 13981 15251 17119 17711 18407 19043 23407 25877 27323/ ],
 pell       => [ qw/169 385 741 961 1121 2001 3827 4879 5719 6215 6265 6441 6479 6601 7055 7801 8119 9799 10945 11395 13067 13079 13601 15841 18241 19097 20833 20951 24727 27839 27971 29183 29953/ ],
 frobenius  => [ qw/4181 5777 6721 10877 13201 15251 34561 51841 64079 64681 67861 68251 75077 90061 96049 97921 100127/ ],
 frob35     => [ qw/13333 44801 486157 1615681 3125281 4219129 9006401 12589081 13404751 15576571 16719781/ ],
);

if ($use64) {
  push @{$pseudoprimes{psp3}}, 4398117272641;
  push @{$pseudoprimes{3}}, 1099558795087;
  push @{$pseudoprimes{lucas}}, 2199055761527;
  push @{$pseudoprimes{slucas}}, 4294967311,4294967357,12598021314449;
  push @{$pseudoprimes{eslucas}}, 4294967311,4294967357,10099386070337;
  push @{$pseudoprimes{aeslucas1}}, 4294967311,4294967357,10071551814917;
  push @{$pseudoprimes{aeslucas2}}, 34372519409;
}
if (!$usexs) {

  # Always cut down on the testing with PP for these mostly-novelty tests
  $#{$pseudoprimes{catalan}} = 1 unless $extra;
  $#{$pseudoprimes{perrin}} = 6 unless $extra;

  # Maybe cut even more
  if (!$usegmp || !defined &Math::Prime::Util::GMP::binomial || $Math::Prime::Util::GMP::VERSION < 0.27) {
    # Don't make Math::BigInt do large binomials
    $pseudoprimes{catalan} = [5907];
  }

}

my @small_lucas_trials = (2, 9, 16, 100, 102, 2047, 2048, 5781, 9000, 14381);

plan tests => 0 + 2
                + 4
                + scalar(keys %pseudoprimes)
                + scalar @phis
                + 1  # mr base 2    2-4k
                + 9  # mr with large bases
                + 3  # multi-base Fermat/strong pseudoprimes
                + 1  # small extra_strong
                + scalar @small_lucas_trials
                + 1  # frob-underwood
                + 2*$use64  # frob-underwood
                + 1  # frob-khashin
                + 2*$use64  # frob-khashin
                + 1*$extra
                + 6  # Perrin restrictions
                + 3  # Test implicit bases and array of bases
                + 0;

# Enforced by XS prototype.
#ok(!eval { is_strong_pseudoprime(2047); }, "MR with no base fails");
ok(!eval { is_strong_pseudoprime(2047,0); }, "MR base 0 fails");
ok(!eval { is_strong_pseudoprime(2047,1); }, "MR base 1 fails");

is( is_strong_pseudoprime(0, 2), 0, "MR with 0 shortcut composite");
is( is_strong_pseudoprime(1, 2), 0, "MR with 0 shortcut composite");
is( is_strong_pseudoprime(2, 2), 1, "MR with 2 shortcut prime");
is( is_strong_pseudoprime(3, 2), 1, "MR with 3 shortcut prime");

# Check that for each test, small pseudoprimes pass the given test.
for my $base (sort keys %pseudoprimes) {
  my @c = @{$pseudoprimes{$base}};
  my @fails;
  my $text;
  if      ($base =~ /^psp(\d+)/) {
    my $pbase = $1;
    @fails = grep { !is_pseudoprime($_,$pbase) } @c;
    $text = "pseudoprimes base $pbase (i.e. Fermat)";
  } elsif ($base =~ /^epsp(\d+)/) {
    my $pbase = $1;
    @fails = grep { !is_euler_pseudoprime($_,$pbase) } @c;
    $text = "Euler pseudoprimes base $pbase";
  } elsif ($base =~ /^aeslucas(\d+)/) {
    my $inc = $1;
    @fails = grep { !is_almost_extra_strong_lucas_pseudoprime($_,$inc) } @c;
    $text = "almost extra strong Lucas pseudoprimes (inc $inc)";
  } elsif ($base eq 'eslucas') {
    @fails = grep { !is_extra_strong_lucas_pseudoprime($_) } @c;
    $text = "extra strong Lucas pseudoprimes";
  } elsif ($base eq 'slucas') {
    @fails = grep { !is_strong_lucas_pseudoprime($_) } @c;
    $text = "strong Lucas pseudoprimes";
  } elsif ($base eq 'lucas') {
    @fails = grep { !is_lucas_pseudoprime($_) } @c;
    $text = "Lucas pseudoprimes";
  } elsif ($base eq 'plumb') {
    @fails = grep { !is_euler_plumb_pseudoprime($_) } @c;
    $text = "Euler-Plumb pseudoprimes";
  } elsif ($base eq 'perrin') {
    @fails = grep { !is_perrin_pseudoprime($_) } @c;
    $text = "Unrestricted Perrin pseudoprimes";
  } elsif ($base eq 'catalan') {
    @fails = grep { !is_catalan_pseudoprime($_) } @c;
    $text = "Catalan pseudoprimes";
  } elsif ($base eq 'frobenius') {
    @fails = grep { !is_frobenius_pseudoprime($_,1,-1) } @c;
    $text = "Frobenius(1,-1) pseudoprimes";
  } elsif ($base eq 'frob35') {
    @fails = grep { !is_frobenius_pseudoprime($_,3,-5) } @c;
    $text = "Frobenius(3,-5) pseudoprimes";
  } elsif ($base eq 'fibonacci') {
    @fails = grep {
      my $t = (($_%5)==2||($_%5)==3) ? $_+1 : $_-1;
      my $is_fib = !lucasumod(1, -1, $t, $_);
      !$is_fib;
    } @c;
    $text = "Fibonacci pseudoprimes";
  } elsif ($base eq 'pell') {
    if ($] < 5.008) {  # Work around a fault in ancient Perl
      @fails = grep { "" . ((lucasumod(2,-1,$_,$_) - kronecker(2,$_)) % $_) } @c;
    } else {
      @fails = grep { ((lucasumod(2,-1,$_,$_) - kronecker(2,$_)) % $_) } @c;
    }
    $text = "Pell pseudoprimes";
  } else {
    @fails = grep { !is_strong_pseudoprime($_,$base) } @c;
    $text = "strong pseudoprimes base $base (i.e. Miller-Rabin)";
  }

  is_deeply(\@fails, [], "Small $text");
}

# Check that phi_n makes passes MR with all prime bases < pn
for my $phi (1 .. scalar @phis) {
  #next if ($phi > 4) && (!$use64);
  ok( is_strong_pseudoprime($phis[$phi-1], @sp[0 .. $phi-1]), "phi_$phi passes MR with first $phi primes");
}

# Verify MR base 2 for all small numbers
{
  my $mr2fail = 0;
  for (2 .. 4032) {
    next if $_ == 2047 || $_ == 3277;
    if (is_prime($_)) {
      if (!is_strong_pseudoprime($_,2)) { $mr2fail = $_; last; }
    } else {
      if (is_strong_pseudoprime($_,2))  { $mr2fail = $_; last; }
    }
  }
  is($mr2fail, 0, "MR base 2 matches is_prime for 2-4032 (excl 2047,3277)");
}

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

# Some examples of Fermat and strong pseudoprimes
is(is_pseudoprime(143168581, 2, 3, 5, 7, 11), 1, "143168581 is a Fermat pseudoprime to bases 2,3,5,7,11");
is(is_strong_pseudoprime(3215031751, 2, 3, 5, 7), 1, "3215031751 is a strong pseudoprime to bases 2,3,5,7");
is(is_strong_pseudoprime("2152302898747", 2, 3, 5, 7, 11), 1, "2152302898747 is a strong pseudoprime to bases 2,3,5,7,11");

# Verify extra strong for a few small primes
is_deeply( [grep { is_extra_strong_lucas_pseudoprime($_) } 2..100], [grep { $_ >= 2 && $_ <= 100 } @sp], "The first 100 primes are selected by is_extra_strong_lucas_pseudoprime" );

# Verify Lucas for some small numbers
for my $n (@small_lucas_trials) {
  next if $n == 5459 || $n == 5777 || $n == 10877 || $n == 16109 || $n == 18971;
  if (is_prime($n)) {
    # Technically if it is a prime it isn't a pseudoprime.
    ok(is_strong_lucas_pseudoprime($n), "$n is a prime and a strong Lucas-Selfridge pseudoprime");
  } else {
    ok(!is_strong_lucas_pseudoprime($n), "$n is not a prime and not a strong Lucas-Selfridge pseudoprime");
  }
}

# Verify MR base 2-3 for many small numbers (up to phi2)
if ($extra) {
  my $mr2fail = 0;
  foreach my $i (1 .. 50000) {
    my $n = int(rand(1373652)) + 1;
    my $isp23 = !!is_strong_pseudoprime($n,2,3);
    my $prime = !!is_prime($n);
    if ($isp23 != $prime) {
      $mr2fail = $n;
      last;
    }
  }
  is($mr2fail, 0, "is_strong_pseudoprime bases 2,3 matches is_prime");
}

{
  my $fufail = 0;
  my $ntests = ($usexs) ? 100 : 2;
  foreach my $i (1 .. $ntests) {
    my $n = 2*int(rand(1000000000)) + 1;
    my $ispfu = !!is_frobenius_underwood_pseudoprime($n);
    my $prime = !!is_prime($n);
    if ($ispfu != $prime) {
      $fufail = $n;
      last;
    }
  }
  is($fufail, 0, "is_frobenius_underwood_pseudoprime matches is_prime");
  if ($use64) {
    is( is_frobenius_underwood_pseudoprime("2727480595375747"), 1, "Frobenius Underwood with 52-bit prime" );
    is( is_frobenius_underwood_pseudoprime(10099386070337), 0, "Frobenius Underwood with 44-bit Lucas pseudoprime" );
  }
}

{
  my $fufail = 0;
  my $ntests = ($usexs) ? 100 : 2;
  foreach my $i (1 .. $ntests) {
    my $n = 2*int(rand(1000000000)) + 1;
    my $ispfu = !!is_frobenius_khashin_pseudoprime($n);
    my $prime = !!is_prime($n);
    if ($ispfu != $prime) {
      $fufail = $n;
      last;
    }
  }
  is($fufail, 0, "is_frobenius_khashin_pseudoprime matches is_prime");
  if ($use64) {
    is( is_frobenius_khashin_pseudoprime("2727480595375747"), 1, "Frobenius Khashin with 52-bit prime" );
    is( is_frobenius_khashin_pseudoprime(10099386070337), 0, "Frobenius Khashin with 44-bit Lucas pseudoprime" );
  }
}

# Perrin restrictions
is( is_perrin_pseudoprime("40814059160177",0), 1, "40814059160177 is an unrestricted Perrin pseudoprime");
is( is_perrin_pseudoprime("40814059160177",1), 0, "40814059160177 is not a minimal restricted Perrin pseudoprime");
is( is_perrin_pseudoprime("36407440637569",1), 1, "36407440637569 is minimal restricted Perrin pseudoprime");
is( is_perrin_pseudoprime("36407440637569",2), 0, "36407440637569 is not an Adams/Shanks Perrin pseudoprime");
is( is_perrin_pseudoprime("364573433665",2), 1, "364573433665 is an Adams/Shanks Perrin pseudoprime");
is( is_perrin_pseudoprime("364573433665",3), 0, "364573433665 is not a Grantham restricted Perrin pseudoprime");

# Implicit base 2
{
  my(@ns) = (91,133,145,697,703,1065,1541,1891,2047,3277,3781,6697,12403);


  my(@carm) = (1729, 2821, 6601, 8911);

  is_deeply( [map { is_pseudoprime($_, 2) } @ns],
             [map { is_pseudoprime($_)    } @ns],
             "implicit base 2" );

  my(@bases235) = (2,3,5);
  is_deeply( [map { is_pseudoprime($_, 2, 3, 5)   } @carm],
             [map { is_pseudoprime($_, @bases235) } @carm],
             "is_pseudoprime can take array of bases" );

  my(@basesnull) = ();
  is_deeply( [map { is_pseudoprime($_, 2)          } @ns],
             [map { is_pseudoprime($_, @basesnull) } @ns],
             "empty array of bases is implicit base 2" );
}
