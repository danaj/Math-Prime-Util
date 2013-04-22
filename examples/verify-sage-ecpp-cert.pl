#!/usr/bin/env perl
use warnings;
use strict;
use Math::BigInt try=>"GMP,Pari";
use Math::Prime::Util qw/:all/;
use Data::Dump qw/dump/;

# Takes the output of one of the Sage functions:
#     goldwasser_kilian(n)
#     atkin_morain(n)
#     ecpp(n)
# and run it through MPU's verifier.
#
# Example:
#   perl verify-sage-ecpp-cert.pl <<EOP
#   [100000000000000000039, 0, 100000000000000000038, 99999999982060129348, 964019216509,
#   (19568017359532265472 : 20079688430965599468 : 1), 964019216509, 0, 858481324253L, 964020818668, 241005204667,
#   (799886189130 : 601302851255 : 1), 241005204667, 0, 11167880189, 241006096857, 626011,
#   (140197349256 : 38879868868 : 1)]
#   EOP

my @input;
while (<>) {
  chomp;
  push @input, split(/\s*,\s*/, $_);
}
my $N;
my ($n, $a, $b, $m, $q, $Px, $Py);
my @cert;
while (@input) {
  $_ = shift @input;
  $_ =~ s/L\s*$//;
  if (!defined $N) {
    die "Need brackets around certificate\n" unless s/^\s*\[//;
    $N = $_;
    $n = $_;
    @cert = ($N, "AGKM");
  }
  elsif (!defined $n) { $n = $_; }
  elsif (!defined $a) { $a = $_; }
  elsif (!defined $b) { $b = $_; }
  elsif (!defined $m) { $m = $_; }
  elsif (!defined $q) { $q = $_; }
  elsif (!defined $Px) {
    die "Can't parse point" unless /\(\s*(\d+)\s*:\s*(\d+)\s*:\s*(\d+)\s*\)/;
    $Px = $1;
    $Py = $2;
    die "Bad input\n"
        unless defined $n && defined $a && defined $b && defined $m
            && defined $q && defined $Px && defined $Py;
    push @cert, [$n, $a, $b, $m, $q, [$Px,$Py]];
    undef $n;
    undef $a;
    undef $b;
    undef $m;
    undef $q;
    undef $Px;
    undef $Py;
  }
}

print dump(\@cert), "\n";
print verify_prime(@cert) ? "SUCCESS\n" : "FAILURE\n";
