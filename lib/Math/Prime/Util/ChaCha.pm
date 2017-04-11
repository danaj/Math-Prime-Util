package Math::Prime::Util::ChaCha;
use strict;
use warnings;
use Carp qw/carp croak confess/;

BEGIN {
  $Math::Prime::Util::ChaCha::AUTHORITY = 'cpan:DANAJ';
  $Math::Prime::Util::ChaCha::VERSION = '0.61';
}

###############################################################################
# Begin ChaCha core, reference RFC 7539
# with change to make blockcount/nonce be 64/64 from 32/96
# Dana Jacobsen, 9 Apr 2017

BEGIN {
  use constant ROUNDS => 20;
}

# Note: 32-bit will break "normal" code.
# We 'use integer' to fix one problem, then do extra masking and
# a final signed->unsigned conversion to fix the extra problems our
# fix gave us.  We really want 'use uinteger'.

sub _quarterround {
  my($a,$b,$c,$d) = @_;
  use integer;
  $a=($a+$b)&0xFFFFFFFF; $d^=$a; $d=(($d<<16)&0xFFFFFFFF)|(($d>>16)& 0xFFFF);
  $c=($c+$d)&0xFFFFFFFF; $b^=$c; $b=(($b<<12)&0xFFFFFFFF)|(($b>>20)& 0xFFF);
  $a=($a+$b)&0xFFFFFFFF; $d^=$a; $d=(($d<< 8)&0xFFFFFFFF)|(($d>>24)& 0xFF);
  $c=($c+$d)&0xFFFFFFFF; $b^=$c; $b=(($b<< 7)&0xFFFFFFFF)|(($b>>25)& 0x7F);
  unpack("L*",pack("L*",$a,$b,$c,$d));
}

sub _test_qr {
  return unless ROUNDS == 20;
  my($a,$b,$c,$d);

  ($a,$b,$c,$d) = _quarterround(0x11111111,0x01020304,0x9b8d6f43,0x01234567);
  #($a,$b,$c,$d) = unpack("L*",pack("L*",$a,$b,$c,$d));
  #printf "  %08x  %08x  %08x  %08x\n", $a,$b,$c,$d;
  die "QR test 2.1.1 fail 1" unless $a == 0xea2a92f4 && $b == 0xcb1cf8ce && $c == 0x4581472e && $d == 0x5881c4bb;

  ($a,$b,$c,$d) = _quarterround(0x516461b1,0x2a5f714c, 0x53372767, 0x3d631689);
  #printf "  %08x  %08x  %08x  %08x\n", $a,$b,$c,$d;
  die "QR test 2.2.1 fail 2" unless $a == 0xbdb886dc && $b == 0xcfacafd2 && $c == 0xe46bea80 && $d == 0xccc07c79;
}
_test_qr;

#  State is:
#       cccccccc  cccccccc  cccccccc  cccccccc
#       kkkkkkkk  kkkkkkkk  kkkkkkkk  kkkkkkkk
#       kkkkkkkk  kkkkkkkk  kkkkkkkk  kkkkkkkk
#       bbbbbbbb  nnnnnnnn  nnnnnnnn  nnnnnnnn
#
#     c=constant k=key b=blockcount n=nonce

sub _core {
  my($j) = @_;
  use integer;
  #die "Invalid ChaCha state" unless scalar(@$j) == 16;
  my @x = @$j;
  for (1 .. ROUNDS/2) {

    # Unrolling is ugly but makes a big performance diff.
    # Generated with unroll-chacha.pl

    #@x[ 0, 4, 8,12] = _quarterround(@x[ 0, 4, 8,12]);
    #@x[ 1, 5, 9,13] = _quarterround(@x[ 1, 5, 9,13]);
    #@x[ 2, 6,10,14] = _quarterround(@x[ 2, 6,10,14]);
    #@x[ 3, 7,11,15] = _quarterround(@x[ 3, 7,11,15]);
    #@x[ 0, 5,10,15] = _quarterround(@x[ 0, 5,10,15]);
    #@x[ 1, 6,11,12] = _quarterround(@x[ 1, 6,11,12]);
    #@x[ 2, 7, 8,13] = _quarterround(@x[ 2, 7, 8,13]);
    #@x[ 3, 4, 9,14] = _quarterround(@x[ 3, 4, 9,14]);

    $x[ 0]=($x[ 0]+$x[ 4])&0xFFFFFFFF; $x[12]^=$x[ 0]; $x[12]=(($x[12]<<16)&0xFFFFFFFF)|(($x[12]>>16)& 0xFFFF);
    $x[ 8]=($x[ 8]+$x[12])&0xFFFFFFFF; $x[ 4]^=$x[ 8]; $x[ 4]=(($x[ 4]<<12)&0xFFFFFFFF)|(($x[ 4]>>20)& 0xFFF);
    $x[ 0]=($x[ 0]+$x[ 4])&0xFFFFFFFF; $x[12]^=$x[ 0]; $x[12]=(($x[12]<< 8)&0xFFFFFFFF)|(($x[12]>>24)& 0xFF);
    $x[ 8]=($x[ 8]+$x[12])&0xFFFFFFFF; $x[ 4]^=$x[ 8]; $x[ 4]=(($x[ 4]<< 7)&0xFFFFFFFF)|(($x[ 4]>>25)& 0x7F);

    $x[ 1]=($x[ 1]+$x[ 5])&0xFFFFFFFF; $x[13]^=$x[ 1]; $x[13]=(($x[13]<<16)&0xFFFFFFFF)|(($x[13]>>16)& 0xFFFF);
    $x[ 9]=($x[ 9]+$x[13])&0xFFFFFFFF; $x[ 5]^=$x[ 9]; $x[ 5]=(($x[ 5]<<12)&0xFFFFFFFF)|(($x[ 5]>>20)& 0xFFF);
    $x[ 1]=($x[ 1]+$x[ 5])&0xFFFFFFFF; $x[13]^=$x[ 1]; $x[13]=(($x[13]<< 8)&0xFFFFFFFF)|(($x[13]>>24)& 0xFF);
    $x[ 9]=($x[ 9]+$x[13])&0xFFFFFFFF; $x[ 5]^=$x[ 9]; $x[ 5]=(($x[ 5]<< 7)&0xFFFFFFFF)|(($x[ 5]>>25)& 0x7F);

    $x[ 2]=($x[ 2]+$x[ 6])&0xFFFFFFFF; $x[14]^=$x[ 2]; $x[14]=(($x[14]<<16)&0xFFFFFFFF)|(($x[14]>>16)& 0xFFFF);
    $x[10]=($x[10]+$x[14])&0xFFFFFFFF; $x[ 6]^=$x[10]; $x[ 6]=(($x[ 6]<<12)&0xFFFFFFFF)|(($x[ 6]>>20)& 0xFFF);
    $x[ 2]=($x[ 2]+$x[ 6])&0xFFFFFFFF; $x[14]^=$x[ 2]; $x[14]=(($x[14]<< 8)&0xFFFFFFFF)|(($x[14]>>24)& 0xFF);
    $x[10]=($x[10]+$x[14])&0xFFFFFFFF; $x[ 6]^=$x[10]; $x[ 6]=(($x[ 6]<< 7)&0xFFFFFFFF)|(($x[ 6]>>25)& 0x7F);

    $x[ 3]=($x[ 3]+$x[ 7])&0xFFFFFFFF; $x[15]^=$x[ 3]; $x[15]=(($x[15]<<16)&0xFFFFFFFF)|(($x[15]>>16)& 0xFFFF);
    $x[11]=($x[11]+$x[15])&0xFFFFFFFF; $x[ 7]^=$x[11]; $x[ 7]=(($x[ 7]<<12)&0xFFFFFFFF)|(($x[ 7]>>20)& 0xFFF);
    $x[ 3]=($x[ 3]+$x[ 7])&0xFFFFFFFF; $x[15]^=$x[ 3]; $x[15]=(($x[15]<< 8)&0xFFFFFFFF)|(($x[15]>>24)& 0xFF);
    $x[11]=($x[11]+$x[15])&0xFFFFFFFF; $x[ 7]^=$x[11]; $x[ 7]=(($x[ 7]<< 7)&0xFFFFFFFF)|(($x[ 7]>>25)& 0x7F);

    $x[ 0]=($x[ 0]+$x[ 5])&0xFFFFFFFF; $x[15]^=$x[ 0]; $x[15]=(($x[15]<<16)&0xFFFFFFFF)|(($x[15]>>16)& 0xFFFF);
    $x[10]=($x[10]+$x[15])&0xFFFFFFFF; $x[ 5]^=$x[10]; $x[ 5]=(($x[ 5]<<12)&0xFFFFFFFF)|(($x[ 5]>>20)& 0xFFF);
    $x[ 0]=($x[ 0]+$x[ 5])&0xFFFFFFFF; $x[15]^=$x[ 0]; $x[15]=(($x[15]<< 8)&0xFFFFFFFF)|(($x[15]>>24)& 0xFF);
    $x[10]=($x[10]+$x[15])&0xFFFFFFFF; $x[ 5]^=$x[10]; $x[ 5]=(($x[ 5]<< 7)&0xFFFFFFFF)|(($x[ 5]>>25)& 0x7F);

    $x[ 1]=($x[ 1]+$x[ 6])&0xFFFFFFFF; $x[12]^=$x[ 1]; $x[12]=(($x[12]<<16)&0xFFFFFFFF)|(($x[12]>>16)& 0xFFFF);
    $x[11]=($x[11]+$x[12])&0xFFFFFFFF; $x[ 6]^=$x[11]; $x[ 6]=(($x[ 6]<<12)&0xFFFFFFFF)|(($x[ 6]>>20)& 0xFFF);
    $x[ 1]=($x[ 1]+$x[ 6])&0xFFFFFFFF; $x[12]^=$x[ 1]; $x[12]=(($x[12]<< 8)&0xFFFFFFFF)|(($x[12]>>24)& 0xFF);
    $x[11]=($x[11]+$x[12])&0xFFFFFFFF; $x[ 6]^=$x[11]; $x[ 6]=(($x[ 6]<< 7)&0xFFFFFFFF)|(($x[ 6]>>25)& 0x7F);

    $x[ 2]=($x[ 2]+$x[ 7])&0xFFFFFFFF; $x[13]^=$x[ 2]; $x[13]=(($x[13]<<16)&0xFFFFFFFF)|(($x[13]>>16)& 0xFFFF);
    $x[ 8]=($x[ 8]+$x[13])&0xFFFFFFFF; $x[ 7]^=$x[ 8]; $x[ 7]=(($x[ 7]<<12)&0xFFFFFFFF)|(($x[ 7]>>20)& 0xFFF);
    $x[ 2]=($x[ 2]+$x[ 7])&0xFFFFFFFF; $x[13]^=$x[ 2]; $x[13]=(($x[13]<< 8)&0xFFFFFFFF)|(($x[13]>>24)& 0xFF);
    $x[ 8]=($x[ 8]+$x[13])&0xFFFFFFFF; $x[ 7]^=$x[ 8]; $x[ 7]=(($x[ 7]<< 7)&0xFFFFFFFF)|(($x[ 7]>>25)& 0x7F);

    $x[ 3]=($x[ 3]+$x[ 4])&0xFFFFFFFF; $x[14]^=$x[ 3]; $x[14]=(($x[14]<<16)&0xFFFFFFFF)|(($x[14]>>16)& 0xFFFF);
    $x[ 9]=($x[ 9]+$x[14])&0xFFFFFFFF; $x[ 4]^=$x[ 9]; $x[ 4]=(($x[ 4]<<12)&0xFFFFFFFF)|(($x[ 4]>>20)& 0xFFF);
    $x[ 3]=($x[ 3]+$x[ 4])&0xFFFFFFFF; $x[14]^=$x[ 3]; $x[14]=(($x[14]<< 8)&0xFFFFFFFF)|(($x[14]>>24)& 0xFF);
    $x[ 9]=($x[ 9]+$x[14])&0xFFFFFFFF; $x[ 4]^=$x[ 9]; $x[ 4]=(($x[ 4]<< 7)&0xFFFFFFFF)|(($x[ 4]>>25)& 0x7F);

  }
  pack("L16", map { $x[$_] + $j->[$_] } 0..15);
}
sub _test_core {
  return unless ROUNDS == 20;
  my $init_state = '617078653320646e79622d326b20657403020100070605040b0a09080f0e0d0c13121110171615141b1a19181f1e1d1c00000001090000004a00000000000000';
  my @state = map { hex("0x$_") } unpack "(a8)*", $init_state;
  my $instr = join("",map { sprintf("%08x",$_) } @state);
  die "Block function fail test 2.3.2 input" unless $instr eq '617078653320646e79622d326b20657403020100070605040b0a09080f0e0d0c13121110171615141b1a19181f1e1d1c00000001090000004a00000000000000';
  my @out = unpack("L16", _core(\@state));
  my $outstr = join("",map { sprintf("%08x",$_) } @out);
  #printf "  %08x  %08x  %08x  %08x\n  %08x  %08x  %08x  %08x\n  %08x  %08x  %08x  %08x\n  %08x  %08x  %08x  %08x\n", @state;
  die "Block function fail test 2.3.2 output" unless $outstr eq 'e4e7f11015593bd11fdd0f50c47120a3c7f4d1c70368c0339aaa22044e6cd4c3466482d209aa9f0705d7c214a2028bd9d19c12b5b94e16dee883d0cb4e3c50a2';
}
_test_core();

sub _keystream {
  my($nbytes, $rstate) = @_;
  croak "Keystream invalid state" unless scalar(@$rstate) == 16;
  my $stream = '';
  while ($nbytes > 0) {
    $stream .= _core($rstate);
    $nbytes -= 64;
    if (++$rstate->[12] > 4294967295) {
      $rstate->[12] = 0;  $rstate->[13]++;
    }
  }
  return $stream;
}
sub _test_keystream {
  return unless ROUNDS == 20;
  my $init_state = '617078653320646e79622d326b20657403020100070605040b0a09080f0e0d0c13121110171615141b1a19181f1e1d1c00000001000000004a00000000000000';
  my @state = map { hex("0x$_") } unpack "(a8)*", $init_state;
  my $instr = join("",map { sprintf("%08x",$_) } @state);
  croak "Block function fail test 2.4.2 input" unless $instr eq '617078653320646e79622d326b20657403020100070605040b0a09080f0e0d0c13121110171615141b1a19181f1e1d1c00000001000000004a00000000000000';
  my $keystream = _keystream(114, \@state);
  # Verify new state
  my $outstr = join("",map { sprintf("%08x",$_) } @state);
  croak "Block function fail test 2.4.2 output" unless $outstr eq '617078653320646e79622d326b20657403020100070605040b0a09080f0e0d0c13121110171615141b1a19181f1e1d1c00000003000000004a00000000000000';
  # Rather tediaus way of doing it, but very explicit
  my $ksstr = join("",map { sprintf("%02x",$_) } map { ord } split(//,$keystream));
  croak "Block function fail test 2.4.2 keystream" unless substr($ksstr,0,2*114) eq '224f51f3401bd9e12fde276fb8631ded8c131f823d2c06e27e4fcaec9ef3cf788a3b0aa372600a92b57974cded2b9334794cba40c63e34cdea212c4cf07d41b769a6749f3f630f4122cafe28ec4dc47e26d4346d70b98c73f3e9c53ac40c5945398b6eda1a832c89c167eacd901d7e2bf363';
}
_test_keystream();

# End ChaCha core
###############################################################################

# Simple PRNG used to fill small seeds
sub _prng_next {
  my($s) = @_;
  my $oldstate = $s->[0];
  $s->[0] = ($s->[0] * 747796405 + $s->[1]) & 0xFFFFFFFF;
  my $word = ((($oldstate >> (($oldstate >> 28) + 4)) ^ $oldstate) * 277803737) & 0xFFFFFFFF;
  ($word >> 22) ^ $word;
}
sub _prng_new {
  my($a,$b,$c,$d) = @_;
  my @s = (0, (($b << 1) | 1) & 0xFFFFFFFF);
  _prng_next(\@s);
  $s[0] = ($s[0] + $a) & 0xFFFFFFFF;
  _prng_next(\@s);
  \@s;
}
###############################################################################

my $_goodseed;
my $_state;     # Not accessible outside this file by standard means
my $_stream;
my $_have;
my $_sptr;

sub _is_csprng_well_seeded { $_goodseed }

sub seed_csprng {
  my($seed) = @_;
  $_goodseed = length($seed) >= 16;
  my @seed = unpack("L*",substr($seed,0,40));
  # If not enough data, fill rest using simple RNG
  if ($#seed < 9) {
    my $rng = _prng_new(map { $_ <= $#seed ? $seed[$_] : 0 } 0..3);
    push @seed, _prng_next($rng) while $#seed < 9;
  }
  croak "Seed count failure" unless $#seed == 9;
  $_state = [0x61707865, 0x3320646e, 0x79622d32, 0x6b206574, @seed[0..7], 0, 0, @seed[8..9]];
  $_stream = '';
  $_have = 0;
  $_sptr = 0;
}
sub srand {
  my $seed = shift;
  $seed = CORE::rand unless defined $seed;
  my $str = (~0 == 4294967295)
          ? pack("L",$seed)
          : pack("L2", $seed, $seed >> 32);
  seed_csprng($str);
  $seed;
}
sub irand {
  if ($_have < 4) {
    my $keystream = _keystream(64, $_state);
    $_stream = substr($_stream,$_sptr,$_have) . $keystream;
    $_sptr = 0;
    $_have = length($_stream);
  }
  $_have -= 4;
  $_sptr += 4;
  return unpack("L",substr($_stream, $_sptr-4, 4));
}
sub irand64 {
  return irand() if ~0 == 4294967295;
  if ($_have < 8) {
    my $keystream = _keystream(64, $_state);
    $_stream = substr($_stream,$_sptr,$_have) . $keystream;
    $_sptr = 0;
    $_have = length($_stream);
  }
  $_have -= 8;
  $_sptr += 8;
  return unpack("Q", substr($_stream, $_sptr-8, 8));
}
sub random_bytes {
  my($bytes) = @_;
  $bytes = (defined $bytes) ? int abs $bytes : 0;

  if ($_have < $bytes) {
    return _keystream($bytes, $_state) if $_have == 0;
    my($s,$h)=($_sptr,$_have);
    ($_sptr,$_have)=($s+$h,0);
    return substr($_stream,$s,$h) . _keystream($bytes-$h, $_state);
  }
  $_have -= $bytes;
  $_sptr += $bytes;
  return substr($_stream, $_sptr-$bytes, $bytes);
}

1;

__END__


# ABSTRACT:  Pure Perl ChaCha20 CSPRNG

=pod

=encoding utf8

=head1 NAME

Math::Prime::Util::ChaCha - Pure Perl ChaCha20 CSPRNG


=head1 VERSION

Version 0.61


=head1 SYNOPSIS

=head1 DESCRIPTION

A pure Perl implementation of ChaCha20 with a CSPRNG interface.

=head1 FUNCTIONS

=head2 seed_csprng

Takes a binary string as input and seeds the internal CSPRNG.

=head2 srand

A method for sieving the CSPRNG with a small value.  This will not be secure
but can be useful for simulations and emulating the system C<srand>.

With no argument, chooses a random number, seeds and returns the number.
With a single integer argument, seeds and returns the number.

=head2 irand

Returns a random 32-bit integer.

=head2 irand64

Returns a random 64-bit integer.

=head2 random_bytes

Takes an unsigned number C<n> as input and returns that many random bytes
as a single binary string.

=head2

=head1 AUTHORS

Dana Jacobsen E<lt>dana@acm.orgE<gt>

=head1 ACKNOWLEDGEMENTS

Daniel J. Bernstein wrote the ChaCha family of stream ciphers in 2008 as
an update to the popular Salsa20 cipher from 2005.

RFC7539: "ChaCha20 and Poly1305 for IETF Protocols" was used to create both
the C and Perl implementations.  Test vectors from that document are used
here as well.

=head1 COPYRIGHT

Copyright 2017 by Dana Jacobsen E<lt>dana@acm.orgE<gt>

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut
