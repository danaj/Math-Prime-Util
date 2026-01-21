#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_prime logint verify_prime
                         random_maurer_prime_with_cert
                         random_shawe_taylor_prime_with_cert
                         random_proven_prime_with_cert
                        /;

my $usegmp= Math::Prime::Util::prime_get_config->{'gmp'};
my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
$use64 = 0 if $use64 && 18446744073709550592 == ~0;

my $do_st = 1;
$do_st = 0 unless eval { require Digest::SHA;
                         my $version = $Digest::SHA::VERSION;
                         $version =~ s/[^\d.]//g;
                         $version >= 4.00; };

plan tests => 3*2;

my $bits = $usegmp          ?  80
         : ~0 > 4294967295  ?  67
                            :  35;

{
  my($n,$cert) = random_maurer_prime_with_cert($bits);
  ok( is_prime($n) && logint($n,2) == $bits-1, "Random Maurer prime returns a $bits-bits prime" );
  ok( verify_prime($cert), "   with a valid certificate" );
}

SKIP: {
  skip "random Shawe-Taylor prime generation requires Digest::SHA",2 unless $do_st;
  my($n,$cert) = random_shawe_taylor_prime_with_cert($bits);
  ok( is_prime($n) && logint($n,2) == $bits-1, "Random Shawe-Taylor prime returns a $bits-bits prime" );
  ok( verify_prime($cert), "   with a valid certificate" );
}

{
  my($n,$cert) = random_proven_prime_with_cert($bits);
  ok( is_prime($n) && logint($n,2) == $bits-1, "Random proven prime returns a $bits-bits prime" );
  ok( verify_prime($cert), "   with a valid certificate" );
}
