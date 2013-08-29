#!/usr/bin/env perl
use warnings;
use strict;
use Math::Prime::Util ':all';
use Math::BigInt lib=>"GMP";

$|++;

# The number of tests performed.  71 makes a nice display for 80 columns.
my $num = 71;
# Select random primes with sizes randomly between 4 and this number of bits.
my $size = 300;
# Which selection method?
#    mpu is 2x faster than pari, but it's our code
#    pari works pretty well, and is 2x faster than Crypt::Primes
#    cpmaurer is slow and can produce composites
my $prime_method = 'pari';   # mpu, pari, or cpmaurer

my @ns;
print "Generate ";
die "invalid size, must be > 4" unless $size > 4;
foreach my $i (1..$num) {
  my $bits = int(rand($size-4)) + 4;
  my $n;

  if      ($prime_method eq 'cpmaurer') {
    require Crypt::Primes;
    $n = Crypt::Primes::maurer(Size=>$bits); 
  } elsif ($prime_method eq 'pari') {
    require Math::Pari;
    require Crypt::Random;
    do { $n = Crypt::Random::makerandom(Size=>$bits,Strength=>0); }
       while !Math::Pari::isprime($n);
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
prime_set_config(verbose=>1);
foreach my $certn (@certs) {
  my $v = verify_prime($certn->[1]);
  print proof_mark($certn->[1]);
  next if $v;
  print "\n\n$certn->[0] didn't verify!\n\n";
  {
    my $c = $certn->[1];
    $c =~ s/^/  /smg;
    print $c;
  }
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
    return 'E' if $cert =~ /Type\s+ECPP/;
    ($type) = $cert =~ /Type (\S+)/;
  }
  if (!defined $type) { die "\nNo type:\n\n$cert"; }
  if    ($type =~ /bls5/i)      { return '5'; }
  elsif ($type =~ /bls7/i)      { return '7'; }
  if    ($type =~ /bls3/i)      { return '-'; }
  elsif ($type =~ /bls15/i)     { return '+'; }
  elsif ($type =~ /bpsw|small/i){ return '.'; }
  elsif ($type =~ /ecpp|agkm/i) { return 'E'; }
  warn "type: $type\n";
  return '?';
}
