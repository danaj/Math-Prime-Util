#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/omega_prime_count prime_get_config/;

my $config = prime_get_config();
plan skip_all => "large native omega_prime_count test requires 64-bit XS"
  unless $config->{'xs'} && $config->{'maxbits'} > 32;

diag "checking final-prime-power regression; this may take about 35 seconds";
is(omega_prime_count(7, "27111750936270"), "400287212536",
   "omega_prime_count includes final prime powers beyond its prime cache");

done_testing();
