#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/next_prime prev_prime/;

my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;

plan tests => 2 + 3*2 + 6 + 1 + 2 + 148 + 148 + 1;

my @small_primes = qw/
2 3 5 7 11 13 17 19 23 29 31 37 41 43 47 53 59 61 67 71
73 79 83 89 97 101 103 107 109 113 127 131 137 139 149 151 157 163 167 173
179 181 191 193 197 199 211 223 227 229 233 239 241 251 257 263 269 271 277 281
283 293 307 311 313 317 331 337 347 349 353 359 367 373 379 383 389 397 401 409
419 421 431 433 439 443 449 457 461 463 467 479 487 491 499 503 509 521 523 541
547 557 563 569 571 577 587 593 599 601 607 613 617 619 631 641 643 647 653 659
661 673 677 683 691 701 709 719 727 733 739 743 751 757 761 769 773 787 797 809
811 821 823 827 829 839 853 857 859 863 877 881 883 887 907 911 919 929 937 941
947 953 967 971 977 983 991 997 1009 1013 1019 1021 1031 1033 1039 1049 1051 1061 1063 1069
1087 1091 1093 1097 1103 1109 1117 1123 1129 1151 1153 1163 1171 1181 1187 1193 1201 1213 1217 1223
1229 1231 1237 1249 1259 1277 1279 1283 1289 1291 1297 1301 1303 1307 1319 1321 1327 1361 1367 1373
1381 1399 1409 1423 1427 1429 1433 1439 1447 1451 1453 1459 1471 1481 1483 1487 1489 1493 1499 1511
1523 1531 1543 1549 1553 1559 1567 1571 1579 1583 1597 1601 1607 1609 1613 1619 1621 1627 1637 1657
1663 1667 1669 1693 1697 1699 1709 1721 1723 1733 1741 1747 1753 1759 1777 1783 1787 1789 1801 1811
1823 1831 1847 1861 1867 1871 1873 1877 1879 1889 1901 1907 1913 1931 1933 1949 1951 1973 1979 1987
1993 1997 1999 2003 2011 2017 2027 2029 2039 2053 2063 2069 2081 2083 2087 2089 2099 2111 2113 2129
2131 2137 2141 2143 2153 2161 2179 2203 2207 2213 2221 2237 2239 2243 2251 2267 2269 2273 2281 2287
2293 2297 2309 2311 2333 2339 2341 2347 2351 2357 2371 2377 2381 2383 2389 2393 2399 2411 2417 2423
2437 2441 2447 2459 2467 2473 2477 2503 2521 2531 2539 2543 2549 2551 2557 2579 2591 2593 2609 2617
2621 2633 2647 2657 2659 2663 2671 2677 2683 2687 2689 2693 2699 2707 2711 2713 2719 2729 2731 2741
2749 2753 2767 2777 2789 2791 2797 2801 2803 2819 2833 2837 2843 2851 2857 2861 2879 2887 2897 2903
2909 2917 2927 2939 2953 2957 2963 2969 2971 2999 3001 3011 3019 3023 3037 3041 3049 3061 3067 3079
3083 3089 3109 3119 3121 3137 3163 3167 3169 3181 3187 3191 3203 3209 3217 3221 3229 3251 3253 3257
3259 3271 3299 3301 3307 3313 3319 3323 3329 3331 3343 3347 3359 3361 3371 3373 3389 3391 3407 3413
3433 3449 3457 3461 3463 3467 3469 3491 3499 3511 3517 3527 3529 3533 3539 3541 3547 3557 3559 3571
/;

{
  # insert primes before and after
  unshift @small_primes, undef;
  push @small_primes, 3581;
  # Now test next_prime and prev_prime for all numbers 0 to 3572
  my $prev_index = 0;
  my $next_index = 1;
  # We'll collect them here
  my(@got_next, @got_prev, @exp_next, @exp_prev);
  foreach my $n (0 .. 3572) {
    $next_index++ if $n >= $small_primes[$next_index];
    $prev_index++ if $n > $small_primes[$prev_index+1];
    push @got_next, next_prime($n);
    push @got_prev, prev_prime($n);
    push @exp_next, $small_primes[$next_index];
    push @exp_prev, $small_primes[$prev_index];
  }
  is_deeply( \@got_next, \@exp_next, "next_prime 0 .. 3572" );
  is_deeply( \@got_prev, \@exp_prev, "prev_prime 0 .. 3572" );
}

my %primegaps = (
 19609 => 52,
 360653 => 96,
 2010733 => 148,
);

while (my($base, $range) = each (%primegaps)) {
  is( next_prime($base), $base+$range, "next prime of $base is $base+$range" );
  is( prev_prime($base+$range), $base, "prev prime of $base+$range is $base" );
}

is( next_prime(19608), 19609, "next prime of 19608 is 19609" );
is( next_prime(19610), 19661, "next prime of 19610 is 19661" );
is( next_prime(19660), 19661, "next prime of 19660 is 19661" );
is( prev_prime(19662), 19661, "prev prime of 19662 is 19661" );
is( prev_prime(19660), 19609, "prev prime of 19660 is 19609" );
is( prev_prime(19610), 19609, "prev prime of 19610 is 19609" );

is( next_prime(10019), 10037, "next prime of 10019 is 10037" );

is( prev_prime(2), undef, "Previous prime of 2 returns undef" );
if ($use64) {
  # With 5.8.8 and earlier, this can cause problems due to Perl getting lost
  # when the return value is turned into a Math::BigInt.  Fixed in 5.8.9.
  is( "".next_prime(18446744073709551611), "18446744073709551629", "Next prime of ~0-4 returns bigint next prime" );
} else {
  is( "".next_prime(4294967291), "4294967311", "Next prime of ~0-4 returns bigint next prime" );
}

# Turns out the testing of prev/next from 0-3572 still misses some cases.
foreach my $n (2010733 .. 2010880) {
  is(next_prime($n), 2010881, "next_prime($n) == 2010881");
}
foreach my $n (2010734 .. 2010881) {
  is(prev_prime($n), 2010733, "prev_prime($n) == 2010733");
}
# Similar test case to 2010870, where m=0 and next_prime is at m=1
is(next_prime(1234567890), 1234567891, "next_prime(1234567890) == 1234567891)");
