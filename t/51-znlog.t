#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/znlog/;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
my $usegmp= Math::Prime::Util::prime_get_config->{'gmp'};
my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
$use64 = 0 if $use64 && 18446744073709550592 == ~0;

my @znlogs = (
 [ [5,2,1019], 10],
 [ [2,4,17], undef],
 [ [7,3,8], undef],
 [ [7,17,36], undef],       # No solution (Pari #1463)
 [ [1,8,9], 0],
 [ [3,3,8], 1],
 [ [10,2,101], 25],
 [ [2,55,101], 73],         # 2 = 55^73 mod 101
 [ [5,2,401], 48],          # 5 = 2^48 mod 401  (Pari #1285)
 [ [228,2,383], 110],
 [ [3061666278, 499998, 3332205179], 22],
 [ [5678,5,10007], 8620],   # 5678 = 5^8620 mod 10007
 [ [7531,6,8101], 6689],    # 7531 = 6^6689 mod 8101
 # Some odd cases.  Pari pre-2.6 and post 2.6 have issues with them.
 [ [0,30,100], 2],          # 0 = 30^2 mod 100
 [ [1,1,101], 0],           # 1 = 1^0 mod 101
 [ [8,2,102], 3],           # 8 = 2^3 mod 102
 [ [18,18,102], 1],         # 18 = 18^1 mod 102
);
if ($usexs || $extra) {
  push @znlogs, [[5675,5,10000019], 2003974];  # 5675 = 5^2003974 mod 10000019
  push @znlogs, [[18478760,5,314138927], 34034873];
  push @znlogs, [[553521,459996,557057], 15471];
  push @znlogs, [[7443282,4,13524947], 6762454];
}
if ($usexs && $use64) {
  # Nice case for PH
  push @znlogs, [[32712908945642193,5,71245073933756341], 5945146967010377];
}

plan tests => scalar(@znlogs);

###### znlog
foreach my $arg (@znlogs) {
  my($aref, $exp) = @$arg;
  my ($a, $g, $p) = @$aref;
  my $k = znlog($a,$g,$p);
  is( $k, $exp, "znlog($a,$g,$p) = " . ((defined $exp) ? $exp : "<undef>") );
}
