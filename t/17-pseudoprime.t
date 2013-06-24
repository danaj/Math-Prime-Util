#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_prime
                         is_pseudoprime
                         is_strong_pseudoprime
                         is_lucas_pseudoprime
                         is_strong_lucas_pseudoprime
                         is_extra_strong_lucas_pseudoprime
                         is_frobenius_underwood_pseudoprime
                         lucas_sequence/;

my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};

# small primes
my @sp = qw/2 3 5 7 11 13 17 19 23 29 31 37/;
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
 lucas      => [ qw/323 377 1159 1829 3827 5459 5777 9071 9179 10877 11419 11663 13919 14839 16109 16211 18407 18971 19043/ ],
 slucas     => [ qw/5459 5777 10877 16109 18971 22499 24569 25199 40309 58519 75077 97439 100127 113573 115639 130139/ ],
 eslucas    => [ qw/989 3239 5777 10877 27971 29681 30739 31631 39059 72389 73919 75077 100127 113573 125249 137549 137801/ ],
);
my $num_pseudoprimes = 0;
foreach my $ppref (values %pseudoprimes) {
  $num_pseudoprimes += scalar @$ppref;
}
my @small_lucas_trials = (2, 9, 16, 100, 102, 2047, 2048, 5781, 9000, 14381);

my %lucas_sequences = (
  "323 1 1 324" => [0,2,1],
  "323 4 1 324" => [170,308,1],
  "323 4 5 324" => [194,156,115],
  "323 3 1 324" => [0,2,1],
  "323 3 1  81" => [0,287,1],
  "323 5 -1 81" => [153,195,322],
  "49001 25 117 24501" => [20933,18744,19141],
  "18971 10001 -1 4743" => [5866,14421,18970],
);

plan tests => 0 + 3
                + 4
                + $num_pseudoprimes
                + scalar @phis
                + 1  # mr base 2    2-4k
                + 9  # mr with large bases
                + scalar @small_lucas_trials
                + scalar(keys %lucas_sequences)
                + 1  # frob-underwood
                + 1*$extra;

ok(!eval { is_strong_pseudoprime(2047); }, "MR with no base fails");
ok(!eval { is_strong_pseudoprime(2047,0); }, "MR base 0 fails");
ok(!eval { is_strong_pseudoprime(2047,1); }, "MR base 1 fails");

is( is_strong_pseudoprime(0, 2), 0, "MR with 0 shortcut composite");
is( is_strong_pseudoprime(1, 2), 0, "MR with 0 shortcut composite");
is( is_strong_pseudoprime(2, 2), 1, "MR with 2 shortcut prime");
is( is_strong_pseudoprime(3, 2), 1, "MR with 3 shortcut prime");


# Check that each strong pseudoprime base b makes it through MR with that base
while (my($base, $ppref) = each (%pseudoprimes)) {
  foreach my $p (@$ppref) {
    if      ($base eq 'eslucas') {
      ok(is_extra_strong_lucas_pseudoprime($p), "$p is an extra strong Lucas pseudoprime");
    } elsif ($base eq 'slucas') {
      ok(is_strong_lucas_pseudoprime($p), "$p is a strong Lucas-Selfridge pseudoprime");
    } elsif ($base eq 'lucas') {
      ok(is_lucas_pseudoprime($p), "$p is a Lucas-Selfridge pseudoprime");
    } elsif ($base =~ /^psp(\d+)/) {
      my $base = $1;
      ok(is_pseudoprime($p, $base), "$p is a pseudoprime to base $base");
    } else {
      ok(is_strong_pseudoprime($p, $base), "Pseudoprime (base $base) $p passes MR");
    }
  }
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

# Lucas sequences, used for quite a few tests
sub lucas_sequence_to_native {
  map { (ref($_) eq 'Math::BigInt') ? int($_->bstr) : $_ } lucas_sequence(@_);
}
while (my($params, $expect) = each (%lucas_sequences)) {
  is_deeply( [lucas_sequence_to_native(split(' ', $params))], $expect, "Lucas sequence $params" );
}

{
  my $fufail = 0;
  foreach my $i (1 .. 5000) {
    my $n = int(rand(1000000000)) + 1;
    my $ispfu = !!is_frobenius_underwood_pseudoprime($n);
    my $prime = !!is_prime($n);
    if ($ispfu != $prime) {
      $fufail = $n;
      last;
    }
  }
  is($fufail, 0, "is_frobenius_underwood_pseudoprime matches is_prime");
}
