#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_pillai/;
my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
my $usegmp = Math::Prime::Util::prime_get_config->{'gmp'};
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};

my @pillai = (23,29,59,61,67,71,79,83,109,137,139,149,193,227,233,239,251,257,269,271,277,293,307,311,317,359,379,383,389,397,401,419,431,437,449,461,463,467,479,499,503,521,551,557,563,569,571,577,593,599,601,607,613,619,631,641,647,661,673,683,691,709,719,727,733,739,787,797,809,811,823,829,853,857,881,883,887,907,919,947,953,967,983,991);
# This is horribly slow in pure Perl.  Run fewer tests.
$#pillai = 10 unless $usexs;   #if (!$usexs) { $#pillai = ($usegmp) ? 50 : 10; }

plan tests => 0
            + 2
            ;

is(is_pillai(1059511), 16, "1059511 is a Pillai prime");

is_deeply( [grep { is_pillai($_) } 0 .. $pillai[-1]], \@pillai, "is_pillai from -10 to 1000" );
