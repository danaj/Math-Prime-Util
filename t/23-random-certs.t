#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_prime verify_prime
                         random_maurer_prime_with_cert
                         random_shawe_taylor_prime_with_cert
                         random_proven_prime_with_cert
                        /;

use Math::BigInt try=>"GMP,Pari";

plan tests => 3*2;

{
  my($n,$cert) = random_maurer_prime_with_cert(80);
  ok( is_prime($n), "Random Maurer prime returns a prime" );
  ok( verify_prime($cert), "   with a valid certificate" );
}

{
  my($n,$cert) = random_shawe_taylor_prime_with_cert(80);
  ok( is_prime($n), "Random Shawe-Taylor prime returns a prime" );
  ok( verify_prime($cert), "   with a valid certificate" );
}

{
  my($n,$cert) = random_proven_prime_with_cert(80);
  ok( is_prime($n), "Random proven prime returns a prime" );
  ok( verify_prime($cert), "   with a valid certificate" );
}
