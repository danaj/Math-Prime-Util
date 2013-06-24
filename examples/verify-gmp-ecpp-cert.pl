#!/usr/bin/env perl
use warnings;
use strict;
use Math::BigInt try=>"GMP,Pari";
use Math::Prime::Util qw/:all/;
use Data::Dump qw/dumpf/;
my $bifilter = sub { my($ctx, $n) = @_;
                     return {dump=>"$n"} if ref($n) eq "Math::BigInt";
                     undef; };

# Takes the output of GMP-ECPP, creates a certificate in the format used
# by MPU, and runs it through the verifier.
#
# Example:
#
# perl -MMath::Prime::Util -E 'say random_ndigit_prime(60)' | \
# gmp-ecpp -q | \
# perl examples/verify-gmp-eccp-cert.pl


my $early_check = 0;

my $N;
my ($n, $a, $b, $m, $q, $Px, $Py);
my @cert;

while (<>) {
  if (/^N\[(\d+)\]\s*=\s*(\d+)/) {
    $n = $2;
    if ($1 == 0) {
      if (defined $N) {
        # I guess we're done with the last one...
        print verify_prime(@cert) ? "SUCCESS\n" : "FAILURE\n";
      }
      #die "Bad input" if defined $N;
      $N = $n;
      @cert = ($n, "AGKM");
    }
  }
  elsif (/^a\s*=\s*(\d+)/)       { $a = $1; }
  elsif (/^b\s*=\s*(\d+)/)       { $b = $1; }
  elsif (/^m\s*=\s*(\d+)/)       { $m = $1; }
  elsif (/^q\s*=\s*(\d+)/)       { $q = $1; }
  elsif (/^P\s*=\s*\(\s*(\d+)\s*,\s*(\d+)\s*\)/) {
    $Px = $1;
    $Py = $2;
    die "Bad input\n"
        unless defined $N && defined $a && defined $b && defined $m
            && defined $q && defined $Px && defined $Py;

    # If for a given q value, is_prime returns 2, that indicates it can
    # produce an n-1 primality proof very quickly, so we could stop now.
    if ($early_check) {
      my $bq = Math::BigInt->new("$q");
      if (is_prime($bq) == 2) {
        push @cert, [$n, $a, $b, $m, [prime_certificate($bq)], [$Px,$Py]];
        last;
      }
    }
    push @cert, [$n, $a, $b, $m, $q, [$Px,$Py]];
  }
  else {
    undef $N if /^proven prime/;
  }
}

print dumpf(\@cert, $bifilter), "\n";
print verify_prime(@cert) ? "SUCCESS\n" : "FAILURE\n";
