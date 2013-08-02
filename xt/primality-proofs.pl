#!/usr/bin/perl
use warnings;
use strict;
use Math::Prime::Util ':all';
use Math::BigInt lib=>"GMP";
#use Crypt::Primes 'maurer';
use Math::Pari qw/isprime/;
use Crypt::Random 'makerandom';
use Data::Dump::Filtered 'dump_filtered';
my $bifilter = sub { my($ctx, $n) = @_;
                     return {dump=>"$n"} if ref($n) eq "Math::BigInt";
                     undef; };

$|++;
my $num = 71;
my $size = 300;
my $prime_method = 'pari';   # mpu, pari, or cpmaurer

my @ns;
print "Generate ";
die "invalid size, must be > 4" unless $size > 4;
foreach my $i (1..$num) {
  my $bits = int(rand($size-4)) + 4;
  my $n;

  # How do we get random primes?
  #   MPU is the fastest, but it's our own code with identical primality tests.
  #   Pari + Crypt::Random works pretty well if you have them.
  #   Crypt::Primes::maurer will sometimes output composites (!!!).
  if      ($prime_method eq 'cpmaurer') {
    $n = Crypt::Primes::maurer(Size=>$bits); 
  } elsif ($prime_method eq 'pari') {
    do { $n = makerandom(Size=>$bits,Strength=>0); } while !isprime($n);
  } elsif ($prime_method eq 'mpu') {
    $n = random_nbit_prime($bits);
  } else {
    die "Unknown random prime generation method\n";
  }

  push @ns, Math::BigInt->new("$n");
  # print a number corresponding to hundreds of bits
  print int(3.322*length("$n")/100);
}
print "\n";

my @certs;
print "Prove    ";
foreach my $n (@ns) {
  my ($isp,$cert) = is_provable_prime_with_cert($n);
  die "$n is reported as $isp\n" unless $isp == 2;
  push @certs, [$n, $cert];
  print proof_mark($cert);
}
print "\n";

print "Verify   ";
foreach my $certn (@certs) {
  my $v = verify_prime($certn->[1]);
  print proof_mark($certn->[1]);
  next if $v;
  print "\n\n$certn->[0] didn't verify!\n\n";
  print dump_filtered($certn->[1], $bifilter);
  die;
}
print "\n";

sub proof_mark {
  my $cert = shift;
  my $type;
  if (ref($cert) eq 'ARRAY') {
    $type = (scalar @$cert == 1) ? "bpsw" : $cert->[1];
    if ($type =~ /n-1/i) {
      $type = ($cert->[2]->[0] eq 'B') ? 'BLS7' : 'BLS5';
    }
  } else {
    ($type) = $cert =~ /Type (\S+)/;
    $type = 'ECPP' if $cert =~ /Type\s+ECPP/;
  }
  if (!defined $type) { die "\nNo cert:\n\n", dump_filtered($cert, $bifilter); }
  if    ($type =~ /bls5/i)      { return '5'; }
  elsif ($type =~ /bls7/i)      { return '7'; }
  if    ($type =~ /bls3/i)      { return '-'; }
  elsif ($type =~ /bls15/i)     { return '+'; }
  elsif ($type =~ /bpsw/i)      { return '.'; }
  elsif ($type =~ /ecpp|agkm/i) { return 'E'; }
  warn "type: $type\n";
  return '?';
}
