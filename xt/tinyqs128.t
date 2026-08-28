#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw(mulint powint prev_prime subint modint);

plan skip_all => 'tinyqs128 requires XS with uint128_t support'
  unless Math::Prime::Util->can('_XS_tinyqs128') &&
         Math::Prime::Util::_XS_has_uint128();

sub check_factor {
  my ($n, $name, @expected) = @_;
  my $f = Math::Prime::Util::_XS_tinyqs128($n);
  my $expected = !@expected || grep { $f == $_ } @expected;
  ok($f > 1 && $f < $n && modint($n, $f) == 0 && $expected, $name)
    or diag("n=$n factor=$f");
}

for my $bits (33, 40, 48, 56, 64) {
  my $pbits = int($bits / 2);
  my $qbits = $bits - $pbits;
  my $p = prev_prime(subint(powint(2, $pbits), 101 + 17*$bits));
  my $q = prev_prime(subint(powint(2, $qbits), 1001 + 31*$bits));
  check_factor(mulint($p, $q), "balanced semiprime near $bits bits", $p, $q);
}

check_factor('8509504187', '33-bit polynomial-pool regression',
             '65423', '130069');

for my $bits (65, 72, 80, 88, 96, 104, 112, 120, 124, 126, 127, 128) {
  my $pbits = int($bits / 2);
  my $qbits = $bits - $pbits;
  my $p = prev_prime(subint(powint(2, $pbits), 1009 + 17*$bits));
  my $q = prev_prime(subint(powint(2, $qbits), 100003 + 31*$bits));
  check_factor(mulint($p, $q), "balanced semiprime near $bits bits", $p, $q);
}

for my $sizes ([32, 72], [40, 80]) {
  my ($pbits, $qbits) = @$sizes;
  my $p = prev_prime(subint(powint(2, $pbits), 12345));
  my $q = prev_prime(subint(powint(2, $qbits), 54321));
  check_factor(mulint($p, $q), "unbalanced $pbits x $qbits semiprime", $p, $q);
}

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
check_factor('409927641983158062491994823',
             'input-derived random seed regression',
             '33316981933', '12303864822075331');

my %adaptive_relations = (
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
for my $bits (sort { $a <=> $b } keys %adaptive_relations) {
  my $case = 0;
  check_factor($_, "$bits-bit adaptive relation case " . ++$case)
    for @{$adaptive_relations{$bits}};
}

my $p = prev_prime(subint(powint(2, 61), 12345));
check_factor(mulint($p, $p), 'perfect square shortcut', $p);
check_factor(mulint(3, $p), 'small factor 3 shortcut', 3);
check_factor(mulint(5, $p), 'small factor 5 shortcut', 5);
check_factor(mulint(7, $p), 'factor-base divisor shortcut', 7);

done_testing();
