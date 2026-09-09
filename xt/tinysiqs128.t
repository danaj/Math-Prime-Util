#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw(modint mulint powint prev_prime subint);

plan skip_all => 'tinysiqs128 requires XS with uint128_t support'
  unless Math::Prime::Util->can('_XS_tinysiqs128') &&
         Math::Prime::Util::_XS_has_uint128();

sub check_zero {
  my ($n, $name) = @_;
  is(Math::Prime::Util::_XS_tinysiqs128($n), 0, $name);
}

sub check_factor {
  my ($n, $name, @expected) = @_;
  my $f = Math::Prime::Util::_XS_tinysiqs128($n);
  my $expected = !@expected || grep { $f == $_ } @expected;
  ok($f > 1 && $f < $n && modint($n, $f) == 0 && $expected, $name)
    or diag("n=$n factor=$f");
}

check_zero($_, "no proper factor for $_") for qw(0 1 2 3 5 7 97 65537);
check_zero('4294967291', '32-bit prime returns zero');
check_zero(prev_prime(powint(2, 33)), '33-bit prime returns zero');
check_zero(prev_prime(powint(2, 128)), '128-bit prime returns zero');

check_factor(4, 'small even composite', 2);
check_factor(9, 'small square', 3);
check_factor(25, 'small square caught by the front end', 5);
check_factor(77, 'small odd composite', 7, 11);

for my $p (7, 11, 13, 17, 19, 23, 29, 31, 37) {
  check_factor($p * 1009, "mod-30 trial wheel reaches factor $p", $p);
}

for my $bits (33, 36, 37, 41, 42, 49, 50, 52, 64, 65, 80, 81, 95,
              96, 106, 107, 116, 117, 120, 124, 126, 127, 128) {
  my $pbits = int($bits / 2);
  my $qbits = $bits - $pbits;
  my $p = prev_prime(subint(powint(2, $pbits), 101 + 17*$bits));
  my $q = prev_prime(subint(powint(2, $qbits), 1001 + 31*$bits));
  check_factor(mulint($p, $q), "balanced semiprime near $bits bits", $p, $q);
}

# The historical SIQS-path labels for 8509504187 and 36346174480237 record
# where the inputs came from; compact SQUFOF may now intercept them.
check_factor('8509504187', '33-bit polynomial regression',
             '65423', '130069');
check_factor('991947540576721',
             '50-bit compact splitter to SIQS fallback regression');
check_factor('1983486546265867',
             '51-bit compact splitter to SIQS fallback regression');
check_factor('36346174480237', 'q=2 to q=1 recovery regression');
check_factor('85070591730229016725614958824927593363',
             '126-bit extra-relations regression');
check_factor('170141183460466514018999784278050528421',
             '127-bit dependency-margin regression');
check_factor('340282366920923092713875188997681850397',
             '128-bit dependency-margin regression');
check_factor('340282366920938463463374607431768211455',
             'full uint128 input is accepted', 3);
check_factor('18446744073709551616', 'decimal input at 2^64', 2);

my $fifth_base = prev_prime(powint(2, 25));
check_factor(powint($fifth_base, 5), 'large odd perfect power', $fifth_base);

for my $sizes ([17, 111], [32, 72], [33, 95], [40, 80]) {
  my ($pbits, $qbits) = @$sizes;
  my $p = prev_prime(subint(powint(2, $pbits), 12345));
  my $q = prev_prime(subint(powint(2, $qbits), 54321));
  check_factor(mulint($p, $q), "unbalanced $pbits x $qbits semiprime", $p, $q);
}

my $cube_base = prev_prime(powint(2, 42));
check_factor(powint($cube_base, 3), 'large odd perfect cube', $cube_base);

my $square_base = prev_prime(powint(2, 64));
check_factor(powint($square_base, 2), '128-bit prime square', $square_base);

my $repeated_p = prev_prime(powint(2, 30));
my $repeated_q = prev_prime(powint(2, 60));
check_factor(mulint(powint($repeated_p, 2), $repeated_q),
             'mixed repeated factors');

my @three_primes = map { prev_prime(powint(2, $_)) } (38, 41, 44);
check_factor(mulint(mulint($three_primes[0], $three_primes[1]),
                    $three_primes[2]),
             'three distinct large factors');

require Math::BigInt;
my $big_input = Math::BigInt->new('8509504187');
check_factor($big_input, 'Math::BigInt input is parsed exactly',
             '65423', '130069');
check_factor('0000000000008509504187', 'leading-zero decimal input',
             '65423', '130069');

for my $bad ('-1', '340282366920938463463374607431768211456',
             '1.5', 'not-an-integer') {
  my $ok = eval { Math::Prime::Util::_XS_tinysiqs128($bad); 1 };
  ok(!$ok && $@ =~ /must fit uint128_t/,
     "invalid input '$bad' is rejected");
}

for (1 .. 3) {
  check_factor('340282366920923092713875188997681850397',
               "repeated 128-bit call $_");
  check_zero('18446744073709551557', "interleaved prime call $_");
}

done_testing();
