package Math::Prime::Util::RandomPrimesGMP;
use strict;
use warnings;
use Carp qw/carp croak confess/;
use Math::Prime::Util::GMP;
use Math::Prime::Util qw/ verify_prime
                          is_provable_prime is_provable_prime_with_cert /;

BEGIN {
  $Math::Prime::Util::RandomPrimesGMP::AUTHORITY = 'cpan:DANAJ';
  $Math::Prime::Util::RandomPrimesGMP::VERSION = '0.61';
}

BEGIN {
  do { require Math::BigInt;  Math::BigInt->import(try=>"GMP,Pari"); }
    unless defined $Math::BigInt::VERSION;

  if (!Math::Prime::Util::GMP::is_csprng_well_seeded()) {
    require Bytes::Random::Secure::Tiny;
    my $brs = Bytes::Random::Secure::Tiny->new(NonBlocking=>1);
    Math::Prime::Util::GMP::seed_csprng(256, $brs->bytes(256));
  }
}

################################################################################

sub random_prime {
  my($low,$high) = @_;
  Math::Prime::Util::_reftyped($_[1],
    Math::Prime::Util::GMP::random_prime($low,$high)
  );
}

sub random_maurer_prime {
  my($bits) = @_;
  Math::Prime::Util::_reftyped($_[0],
    Math::Prime::Util::GMP::random_maurer_prime($bits)
  );
}

sub random_shawe_taylor_prime {
  my($bits) = @_;
  Math::Prime::Util::_reftyped($_[0],
    Math::Prime::Util::GMP::random_shawe_taylor_prime($bits)
  );
}

sub random_maurer_prime_with_cert {
  my($bits) = @_;
  my($n,$cert) = Math::Prime::Util::GMP::random_maurer_prime_with_cert($bits);
  $n = Math::Prime::Util::_reftyped($_[0], $n);
  ($n,$cert);
}

sub random_shawe_taylor_prime_with_cert {
  my($bits) = @_;
  my($n,$cert) = Math::Prime::Util::GMP::random_shawe_taylor_prime_with_cert($bits);
  $n = Math::Prime::Util::_reftyped($_[0], $n);
  ($n,$cert);
}

sub random_proven_prime_with_cert {
  my $k = shift;
  my($n, $isp, $cert);
  if ($Math::Prime::Util::GMP::VERSION >= 0.43) {
    ($n,$cert) = random_maurer_prime_with_cert($k);
  } else {
    $n = random_nbit_prime($k);
    ($isp, $cert) = is_provable_prime_with_cert($n);
    croak "${k}-bit prime could not be proven" if $isp != 2;
  }
  return ($n, $cert);
}

sub miller_rabin_random {
  my($n, $k, $seed) = @_;
  if (defined $seed) {
    Math::Prime::Util::GMP::miller_rabin_random($n, $k, $seed);
  } else {
    Math::Prime::Util::GMP::miller_rabin_random($n, $k);
  }
}


1;

__END__


# ABSTRACT:  Generate random primes using MPU::GMP

=pod

=encoding utf8

=head1 NAME

Math::Prime::Util::RandomPrimesGMP - Generate random primes using MPU::GMP


=head1 VERSION

Version 0.61


=head1 SYNOPSIS

=head1 DESCRIPTION

Routines to generate random primes.


=head1 RANDOM PRIME FUNCTIONS

=head2 random_prime

Generate a random prime between C<low> and C<high>.

=head2 random_ndigit_prime

Generate a random prime with C<n> digits.  C<n> must be at least 1.

=head2 random_nbit_prime

Generate a random prime with C<n> bits.  C<n> must be at least 2.

=head2 random_strong_prime

Generate a random strong prime with C<n> bits.  C<n> must be at least 128.

=head2 random_maurer_prime

Generate a random proven prime with C<n> bits using Maurer's algorithm.
C<n> must be at least 2.

=head2 random_shawe_taylor_prime

Generate a random proven prime with C<n> bits using Shawe-Taylor's algorithm
from FIPS 186-4.
C<n> must be at least 2.

=head2 random_maurer_prime_with_cert

As L</random_maurer_prime> but also returns a certificate string.

=head2 random_shawe_taylor_prime_with_cert

As L</random_shawe_taylor_prime> but also returns a certificate string.

=head2 random_proven_prime

Generate or construct a random provable prime of C<n> bits.  C<n> must
be at least 2.

=head2 random_proven_prime_with_cert

Generate or construct a random provable prime of C<n> bits.  C<n> must
be at least 2.  Returns a list of two items: the prime and the certificate.


=head1 RANDOM PRIMALITY FUNCTIONS

=head2 miller_rabin_random

Given a number C<n> and a number of tests to perform C<k>, this does C<k>
Miller-Rabin tests on C<n> using randomly selected bases.  The return value
is 1 if all bases are a witness to to C<n>, or 0 if any of them fail.

=head1 SEE ALSO

L<Math::Prime::Util>

=head1 AUTHORS

Dana Jacobsen E<lt>dana@acm.orgE<gt>


=head1 COPYRIGHT

Copyright 2017 by Dana Jacobsen E<lt>dana@acm.orgE<gt>

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut
