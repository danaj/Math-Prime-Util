#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/csrand random_bytes/;

my $use64 = (~0 > 4294967295);
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $maxbits = $use64 ? 64 : 32;

plan tests => 6;

########

my $plaintext = 'Now I’m not gonna make a lot of extravagant claims for this little machine. Sure, it’ll change your whole life for the better, but that’s all.';
my $key = "A lovely day for a ride";
my $nonce = "20170412";

csrand( pack("A32A8",$key, $nonce) );
my $ciphertext = $plaintext ^ random_bytes(length($plaintext));

if (unpack("L",$ciphertext) == 351607655) {
  isnt( $ciphertext, $plaintext, "Ciphertext is probably ChaCha/20 expected result" );
} else {
  isnt( $ciphertext, $plaintext, "We at least vaguely changed the text" );
}

my $key2 = "The city needs a car like a fish needs a bicycle.";
csrand( pack("A32A8",$key2, $nonce) );
my $ciphertext2 = $plaintext ^ random_bytes(length($plaintext));
isnt( $ciphertext2, $plaintext, "We at least vaguely changed the text" );

if (unpack("L",$ciphertext2) == 3391833874) {
  isnt( $ciphertext2, $ciphertext, "Different key makes different ChaCha/20 result" );
} else {
  isnt( $ciphertext2, $ciphertext, "Different key produces different data" );
}

csrand( pack("A32A8",$key, $nonce) );
my $ciphertext3 = $plaintext ^ random_bytes(length($plaintext));
is( $ciphertext3, $ciphertext, "We can reproduce the cipher" );

csrand( pack("A32A8",$key, $nonce) );
my $decodetext = $ciphertext ^ random_bytes(length($ciphertext));
is( $decodetext, $plaintext, "We can decode using the same key." );


csrand( pack("A32A8",$key, "Berlin") );
my $ciphertext4 = $plaintext ^ random_bytes(length($plaintext));
isnt( $ciphertext4, $ciphertext, "Different nonce produces different data" );
