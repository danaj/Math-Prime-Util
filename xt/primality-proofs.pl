#!/usr/bin/env perl
use warnings;
use strict;
use Math::Prime::Util ':all';
use Math::BigInt lib=>"GMP";

if (!Math::Prime::Util::prime_get_config->{gmp}) {
  print "\nYou should install Math::Prime::Util::GMP.\n\n";
}
$|++;

print "random prime proofs: 50, 100, 200, 300, 400 +/- 50 digits\n";
test_proofs(  4, 100, 71, 'mpu');  print "\n";
test_proofs( 50, 150, 71, 'mpu');  print "\n";
test_proofs(150, 250, 71, 'mpu');  print "\n";
test_proofs(250, 350, 71, 'mpu');  print "\n";
test_proofs(350, 450, 71, 'mpu');  print "\n";

# size: random primes with bit sizes randomly between 4 and this number
# num:  this many tests performed.  71 makes a nice 80-column display
# method: how to generate random primes:
#     Ideally we would use some independent code.  Time for one thousand
#     random primes from rand(4-300) or rand(4-600) bits:
#            300bits  600bits  which
#               2sec     6sec  mpu (with mpu::gmp installed)
#              31sec   124sec  pari
#              97sec   254sec  cpmaurer
#     We don't seem to have any practical choice other than MPU's
#     random_nbit_prime as the other random prime code is just so slow.
sub test_proofs {
  my($minsize, $size, $num, $prime_method) = @_;

  if      ($prime_method eq 'cpmaurer') {
    require Crypt::Primes;
  } elsif ($prime_method eq 'pari') {
    require Math::Pari;
    require Crypt::Random;
  } elsif ($prime_method eq 'mpu') {
    # nothing
  } else {
    die "Unknown random prime generation method\n";
  }

  my @ns;
  print "Generate ";
  $minsize = 4 if $minsize < 4;
  $minsize = $size if $minsize > $size;
  die "invalid size, must be > 4" unless $size > 4;
  foreach my $i (1..$num) {
    my $bits = int(rand($size-$minsize)) + $minsize;
    my $n;
    if      ($prime_method eq 'cpmaurer') {
      $n = Crypt::Primes::maurer(Size=>$bits);
    } elsif ($prime_method eq 'pari') {
      # This is ~4x faster, has awful distribution.  Still much slower than MPU.
      # $n = Math::Pari::nextprime( ...makerandom... );
      do { $n = Crypt::Random::makerandom(Size=>$bits,Strength=>0); }
         while !Math::Pari::isprime($n);
    } else {
      $n = random_nbit_prime($bits);
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
  prime_set_config(verbose=>0);
  print "\n";
}

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
