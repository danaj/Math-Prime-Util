#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_congruent_number/;
my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
#my $usegmp = Math::Prime::Util::prime_get_config->{'gmp'};
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};

my @cn200 = (5,6,7,13,14,15,20,21,22,23,24,28,29,30,31,34,37,38,39,41,45,46,47,52,53,54,55,56,60,61,62,63,65,69,70,71,77,78,79,80,84,85,86,87,88,92,93,94,95,96,101,102,103,109,110,111,112,116,117,118,119,120,124,125,126,127,133,134,135,136,137,138,141,142,143,145,148,149,150,151,152,154,156,157,158,159,161,164,165,166,167,173,174,175,180,181,182,183,184,188,189,190,191,194,197,198,199);

my @cn1e6 = (1,5,6,7,9);

plan tests => 0
            + 2
            + 6
            ;

# This covers all simple special cases in the code
is_deeply([grep { is_congruent_number($_) } 1..200], \@cn200, "congruent numbers to 200");

is_deeply([grep { is_congruent_number(1000000+$_) } 1..10], \@cn1e6, "congruent numbers 10^6 + (1..10)");

SKIP: {
  skip "PP doesn't have all the NC families",6 unless $usexs;
  is_deeply([map { is_congruent_number($_) } (1419,7611,840873)],
            [1,1,0],
            "Selected values");

  my @iskra = (qw/2451 7923 8643/);
  is_deeply([map { is_congruent_number($_) } @iskra],
            [map { 0 } @iskra],
            "Non-congruent for Iskra 1996");

  my @rheinholz = (qw/2211 5379 10923/);
  is_deeply([map { is_congruent_number($_) } @rheinholz],
            [map { 0 } @rheinholz],
            "Non-congruent for Reinholz, Spearman, and Yang 2013");

  my @cheng2019 = (qw/1947 3363 5907/);
  is_deeply([map { is_congruent_number($_) } @cheng2019],
            [map { 0 } @cheng2019],
            "Non-congruent for Cheng and Guo 2019");

  my @das = (17*3*409*19,17*3*859*3697,19*409*3697*859,17*3*409*19*3697*859, 5*7*29*79,5*7*821*151,29*79*821*151,5*7*29*79*821*151);
  @das = grep { $_ <= ~0 } @das;  # Only test native integers
  is_deeply([map { is_congruent_number($_) } @das],
            [map { 0 } @das],
            "Non-congruent examples from Das and Saikia 2020 section 5");

  my @cheng2018 = (qw/595 714 795 1290 29715 36330 83130 106530 110922 218139 466570 505155 1529745 2679170 5322345 5676555 12063594 24649707/);
  is_deeply([map { is_congruent_number($_) } @cheng2018],
            [map { 0 } @cheng2018],
            "Non-congruent for Cheng and Guo 2018");
}
