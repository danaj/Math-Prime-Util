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

for my $bits (33, 36, 37, 40, 41, 42, 48, 49, 50, 52, 56, 64, 65, 80,
              81, 95, 96, 106, 107, 116, 117, 120, 124, 126, 127, 128) {
  my $pbits = int($bits / 2);
  my $qbits = $bits - $pbits;
  my $p = prev_prime(subint(powint(2, $pbits), 101 + 17*$bits));
  my $q = prev_prime(subint(powint(2, $qbits), 1001 + 31*$bits));
  check_factor(mulint($p, $q), "balanced semiprime near $bits bits", $p, $q);
}

# Retain the independent balanced corpus from the former tinyqs128 test.
for my $bits (65, 72, 80, 88, 96, 104, 112, 120, 124, 126, 127, 128) {
  my $pbits = int($bits / 2);
  my $qbits = $bits - $pbits;
  my $p = prev_prime(subint(powint(2, $pbits), 1009 + 17*$bits));
  my $q = prev_prime(subint(powint(2, $qbits), 100003 + 31*$bits));
  check_factor(mulint($p, $q), "retained balanced case near $bits bits",
               $p, $q);
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
             '126-bit extra-relations regression',
             '9223372036854217439', '9223372036854727117');
check_factor('170141183460466514018999784278050528421',
             '127-bit dependency-margin regression',
             '9223372036854750727', '18446744073709307123');
check_factor('340282366920923092713875188997681850397',
             '128-bit dependency-margin regression',
             '18446744073708783169', '18446744073709486813');
check_factor('340282366920938463463374607431768211455',
             'full uint128 input is accepted', 3);
check_factor('18446744073709551616', 'decimal input at 2^64', 2);
check_factor('409927641983158062491994823',
             'retained 89-bit unbalanced split regression',
             '33316981933', '12303864822075331');

# These difficult high-end inputs formerly exercised tinyqs128's adaptive
# relation gathering.  Their implementation-specific labels are gone, but
# they remain useful independent health cases for the SIQS tail.
my %retained_high_end = (
  127 => [qw(
    154527564022878220982748729324103304743
    161179732548189293760223240490395881353
    127883814751042649647666397127788550683
    125318416433446384843272724113296422871
    140223444857892303368277270633688022093
    162557045538178298486647651339443005663
    122923458607306217422801229530464182507
    115765532451100685824914918408188032091
  )],
  128 => [qw(
    267212883610666843931009408162450623523
    249916340136674620760793264938393468693
    236865378567171197852729155721764556063
    288119424046590105263765021583065156483
    311864632355300169938435240169024757037
    210556783077736902179313303657825846647
  )],
  126 => [qw(
    83656988431739306362112439976867878317
    64754360628803519037675047708191697231
    83595435226116927807160860895179770021
    84813749665103582857694079360041304941
  )],
  125 => [qw(
    39931171521012620279401119962518358699
    24123049789078476044084420737829877593
    37378683740769755305321586855523751313
  )],
  124 => [qw(
    20319590604088821620875280075918123633
    17529763909310581579258330658921405041
  )],
);
for my $bits (sort { $a <=> $b } keys %retained_high_end) {
  my $case = 0;
  check_factor($_, "retained $bits-bit high-end case " . ++$case)
    for @{$retained_high_end{$bits}};
}

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
