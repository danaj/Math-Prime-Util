package Math::Prime::Util::ISAAC;
use strict;
use warnings;
use Carp qw/carp croak confess/;

BEGIN {
  $Math::Prime::Util::ISAAC::AUTHORITY = 'cpan:DANAJ';
  $Math::Prime::Util::ISAAC::VERSION = '0.62';
}

###############################################################################
#   ISAAC, adapted from Bytes::Random::Secure::Tiny and Math::Random::ISAAC
###############################################################################

sub _isaac {
    my($randcnt, $r, $aa, $bb, $cc, $mm) = @_;
    use integer;
    $cc = ($cc+1) & 0xFFFFFFFF;
    $bb = ($bb+$cc) & 0xFFFFFFFF;
    my ($x, $y, $i4);
    for my $i (0 .. 255) {
        $x = $mm->[$i];
        $i4 = $i % 4;
        if    ($i4 == 0) { $aa = $aa ^ ($aa << 13); }
        elsif ($i4 == 1) { $aa = $aa ^ (($aa >>  6) & 0x03ffffff); }
        elsif ($i4 == 2) { $aa = $aa ^ ($aa <<  2); }
        elsif ($i4 == 3) { $aa = $aa ^ (($aa >> 16) & 0x0000ffff); }
        $aa &= 0xFFFFFFFF;

        $aa             = ($mm->[($i+128)   & 0xFF] + $aa)       & 0xFFFFFFFF;;
        $mm->[$i] =  $y = ($mm->[($x >> 2)  & 0xFF] + $aa + $bb) & 0xFFFFFFFF;
        $r->[$i]  = $bb = ($mm->[($y >> 10) & 0xFF] + $x)        & 0xFFFFFFFF;
    }
    if (2147483647+1 < 0) {  # Bulk convert them to unsigned ints
      @$r = unpack("L*",pack("L*",@$r));
    }
    $randcnt = 0;  # We've got 256 new 32-bit values to use
    return ($randcnt, $r, $aa, $bb, $cc, $mm);
}

sub _randinit {
    use integer;
    my($randcnt, $r, $aa, $bb, $cc, $mm) = @_;
    $aa = $bb = $cc = 0;
    my ($c, $d, $e, $f, $g, $h, $j, $k) = (0x9e3779b9)x8; # The golden ratio.
    for (1..4) {
        $c ^= $d << 11;                     $f += $c;       $d += $e;
        $d ^= 0x3fffffff & ($e >> 2);       $g += $d;       $e += $f;
        $e ^= $f << 8;                      $h += $e;       $f += $g;
        $f ^= 0x0000ffff & ($g >> 16);      $j += $f;       $g += $h;
        $g ^= $h << 10;                     $k += $g;       $h += $j;
        $h ^= 0x0fffffff & ($j >> 4);       $c += $h;       $j += $k;
        $j ^= $k << 8;                      $d += $j;       $k += $c;
        $k ^= 0x007fffff & ($c >> 9);       $e += $k;       $c += $d;
    }
    for (my $i = 0; $i < 256; $i += 8) {
        $c += $r->[$i  ];   $d += $r->[$i+1];
        $e += $r->[$i+2];   $f += $r->[$i+3];
        $g += $r->[$i+4];   $h += $r->[$i+5];
        $j += $r->[$i+6];   $k += $r->[$i+7];
        $c ^= $d << 11;                     $f += $c;       $d += $e;
        $d ^= 0x3fffffff & ($e >> 2);       $g += $d;       $e += $f;
        $e ^= $f << 8;                      $h += $e;       $f += $g;
        $f ^= 0x0000ffff & ($g >> 16);      $j += $f;       $g += $h;
        $g ^= $h << 10;                     $k += $g;       $h += $j;
        $h ^= 0x0fffffff & ($j >> 4);       $c += $h;       $j += $k;
        $j ^= $k << 8;                      $d += $j;       $k += $c;
        $k ^= 0x007fffff & ($c >> 9);       $e += $k;       $c += $d;
        @$mm[$i..$i+7] = ($c,$d,$e,$f,$g,$h,$j,$k);
    }
    for (my $i = 0; $i < 256; $i += 8) {
        $c += $mm->[$i  ];  $d += $mm->[$i+1];
        $e += $mm->[$i+2];  $f += $mm->[$i+3];
        $g += $mm->[$i+4];  $h += $mm->[$i+5];
        $j += $mm->[$i+6];  $k += $mm->[$i+7];
        $c ^= $d << 11;                     $f += $c;       $d += $e;
        $d ^= 0x3fffffff & ($e >> 2);       $g += $d;       $e += $f;
        $e ^= $f << 8;                      $h += $e;       $f += $g;
        $f ^= 0x0000ffff & ($g >> 16);      $j += $f;       $g += $h;
        $g ^= $h << 10;                     $k += $g;       $h += $j;
        $h ^= 0x0fffffff & ($j >> 4);       $c += $h;       $j += $k;
        $j ^= $k << 8;                      $d += $j;       $k += $c;
        $k ^= 0x007fffff & ($c >> 9);       $e += $k;       $c += $d;
        @$mm[$i..$i+7] = ($c,$d,$e,$f,$g,$h,$j,$k);
    }
    my @ctx = _isaac($randcnt, $r, $aa, $bb, $cc, $mm);
    $ctx[0] = 256;  # Force running isaac again first use
    return @ctx;
}

###############################################################################

# Simple PRNG used to fill small seeds
sub _prng_next {
  my($s) = @_;
  my $word;
  my $oldstate = $s->[0];
  if (~0 > 0xFFFFFFFF) {
    $s->[0] = ($s->[0] * 747796405 + $s->[1]) & 0xFFFFFFFF;
    $word = ((($oldstate >> (($oldstate >> 28) + 4)) ^ $oldstate) * 277803737) & 0xFFFFFFFF;
  } else {
    { use integer; $s->[0] = unpack("L",pack("L", $s->[0] * 747796405 + $s->[1] )); }
    $word = (($oldstate >> (($oldstate >> 28) + 4)) ^ $oldstate) & 0xFFFFFFFF;
    { use integer; $word = unpack("L",pack("L", $word * 277803737)); }
  }
  ($word >> 22) ^ $word;
}
sub _prng_new {
  my($a,$b,$c,$d) = @_;
  my @s = (0, (($b & 0x7FFFFFFF) << 1) | 1);
  _prng_next(\@s);
  if (~0 > 0xFFFFFFFF) {
    $s[0] = ($s[0] + $a) & 0xFFFFFFFF;
  } else {
    use integer; $s[0] = unpack("L",pack("L",$s[0] + $a));
  }
  _prng_next(\@s);
  $s[0] = ($s[0] ^ $c) & 0xFFFFFFFF;
  _prng_next(\@s);
  $s[0] = ($s[0] ^ $d) & 0xFFFFFFFF;
  _prng_next(\@s);
  \@s;
}

my $_goodseed = 0;
my @_CTX = (0,[],0,0,0,[]); # (randcnt, randrsl[256], aa, bb, cc, mm[256])

sub _is_csprng_well_seeded { $_goodseed }

sub csrand {
  my($seed) = @_;
  $_goodseed = length($seed) >= 16;
  my @mm = (0) x 256;
  while (length($seed) % 4) { $seed .= pack("C",0); }  # zero pad end word
  my @r = unpack("L*",substr($seed,0,1024));
  # If not enough data, fill rest using simple RNG
  if ($#r < 255) {
    my $rng = _prng_new(map { $_ <= $#r ? $r[$_] : 0 } 0..3);
    push @r, _prng_next($rng) while $#r < 255;
  }
  @_CTX = _randinit(0, \@r, 0, 0, 0, \@mm);
  1;
}

sub srand {
  my $seed = shift;
  $seed = CORE::rand unless defined $seed;
  if ($seed <= 4294967295) { csrand(pack("V",$seed)); }
  else                     { csrand(pack("V2",$seed,$seed>>32)); }
  $seed;
}
sub irand {
  @_CTX = _isaac(@_CTX) if $_CTX[0] > 255;
  return $_CTX[1]->[$_CTX[0]++];
}

sub irand64 {
  return irand() if ~0 == 4294967295;
  my($a,$b);
  if ($_CTX[0] < 254) {
    $a = $_CTX[1]->[$_CTX[0]++];
    $b = $_CTX[1]->[$_CTX[0]++];
  } else {
    ($a,$b) = (irand(), irand());
  }
  ($a << 32) | $b;
}

sub random_bytes {
  my($bytes) = @_;
  $bytes = (defined $bytes) ? int abs $bytes : 0;
  my $str = '';

  while ($bytes >= 4) {
    @_CTX = _isaac(@_CTX) if $_CTX[0] > 255;
    my ($rcnt, $r) = @_CTX;
    my($copywords, $havewords) = ($bytes >> 2, 256 - $rcnt);
    $copywords = $havewords if $havewords < $copywords;
    $str .= pack("L*", @$r[$rcnt .. $rcnt+$copywords-1]);
    $_CTX[0] = $rcnt + $copywords;
    $bytes -= 4*$copywords;
  }
  if ($bytes > 0) {
    my $rem = pack("L", irand());
    $str .= substr($rem, 0, $bytes);
  }
  return $str;
}

1;

__END__


# ABSTRACT:  Pure Perl ISAAC CSPRNG

=pod

=encoding utf8

=head1 NAME

Math::Prime::Util::ISAAC - Pure Perl ISAAC CSPRNG


=head1 VERSION

Version 0.62


=head1 SYNOPSIS

=head1 DESCRIPTION

A pure Perl implementation of ISAAC with a CSPRNG interface.

=head1 FUNCTIONS

=head2 csrand

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

Bob Jenkins wrote ISAAC in 1996, which is a seriously fast CSPRNG.

John Allen did the port to Perl in 2000.

Jonathan Yu released L<Math::Random::ISAAC> in 2009 and has maintained it since.

David Oswald trimmed the code substantially for L<Bytes::Random::Secure::Tiny>.
Our code is based on that code.

=head1 COPYRIGHT

Copyright 2017 by Dana Jacobsen E<lt>dana@acm.orgE<gt>

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut
