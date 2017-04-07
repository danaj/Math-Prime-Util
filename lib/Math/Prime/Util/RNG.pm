package Math::Prime::Util::RNG;
use strict;
use warnings;
use Carp qw/carp croak confess/;

BEGIN {
  $Math::Prime::Util::RNG::AUTHORITY = 'cpan:DANAJ';
  $Math::Prime::Util::RNG::VERSION = '0.61';
}

######################################################################################
# ISAAC
######################################################################################

my $_goodseed = 0;
my @_CTX = (0,[],0,0,0,[]); # (randcnt, randrsl[256], aa, bb, cc, mm[256])

sub _isaac {
    use integer;
    my($randcnt, $r, $aa, $bb, $cc, $mm) = @_;
    $cc = ($cc+1) & 0xFFFFFFFF;
    $bb = ($bb+$cc) & 0xFFFFFFFF;
    my ($x, $y); # temporary storage
    for my $i (0 .. 255) {
        my $x = $mm->[$i];
        my $i4 = $i % 4;
        if    ($i4 == 0) { $aa = $aa ^ ($aa << 13); }
        elsif ($i4 == 1) { $aa = $aa ^ ($aa >>  6); }
        elsif ($i4 == 2) { $aa = $aa ^ ($aa <<  2); }
        elsif ($i4 == 3) { $aa = $aa ^ ($aa >> 16); }
        $aa &= 0xFFFFFFFF;

        $aa             = ($mm->[($i+128)   & 0xFF] + $aa)       & 0xFFFFFFFF;;
        $mm->[$i] =  $y = ($mm->[($x >> 2)  & 0xFF] + $aa + $bb) & 0xFFFFFFFF;
        $r->[$i]  = $bb = ($mm->[($y >> 10) & 0xFF] + $x)        & 0xFFFFFFFF;
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

sub _isaac_seed {
  my($seed) = @_;
  $_goodseed = length($seed) >= 16;
  my @mm = (0) x 256;
  my @r = (0) x 256;
  # Replicate seed and put in randrsl
  if (length($seed) > 0) {
    $seed .= $seed while length($seed) < 1024;
    @r = unpack("L256",$seed);
  }
  @_CTX = _randinit(0, \@r, 0, 0, 0, \@mm);
  1;
}

sub _is_csprng_well_seeded { $_goodseed }

sub seed_csprng {
  my($seed) = @_;
  _isaac_seed($seed);
}
sub srand {
  my $seed = shift;
  $seed = CORE::rand unless defined $seed;
  _isaac_seed(pack("L2", ($seed >> 32) & 0xFFFFFFFF, $seed & 0xFFFFFFFF));
  $seed;
}
sub irand {
  @_CTX = _isaac(@_CTX) if $_CTX[0] > 255;
  return $_CTX[1]->[$_CTX[0]++];
}
sub random_bytes {
  my($bytes) = @_;
  $bytes = (defined $bytes) ? int abs $bytes : 0;
  my $str = '';
  #while ($bytes >= 4) {
  #  $str .= pack("L", irand());
  #  $bytes -= 4;
  #}
  while ($bytes >= 4) {
    @_CTX = _isaac(@_CTX) if $_CTX[0] > 255;
    my ($rcnt, $r) = @_CTX;
    my $remwords = 256 - $rcnt;
    my $copywords = ($remwords > $bytes*4) ? $bytes*4 : $remwords;
    $str .= pack("L*", map { $r->[$rcnt++] } 1 .. $copywords);
    $_CTX[0] = $rcnt;
    $bytes -= 4*$copywords;
  }
  if ($bytes > 0) {
    my $rem = pack("L", irand());
    $str .= substr($rem, 0, $bytes);
  }
  return $str;
}

# TODO: Ensure this is right with 32-bit and big endian
sub irand64 { irand() | (irand() << 32); };
# TODO: 32-bit and check correctness
sub drand {
  my $m = shift;
  my $d = (irand64() / (~0 + 1.0));
  $d *= $m if $m;
  $d;
}
#sub drand { irand() / 4294967296.0 }

1;

__END__


# ABSTRACT:  Pure Perl CSPRNG

=pod

=encoding utf8

=head1 NAME

Math::Prime::Util::RNG - Pure Perl CSPRNG


=head1 VERSION

Version 0.61


=head1 SYNOPSIS

=head1 DESCRIPTION

An alternative pure Perl implementation of a CSPRNG.

=head1 AUTHORS

Dana Jacobsen E<lt>dana@acm.orgE<gt>


=head1 COPYRIGHT

Copyright 2017 by Dana Jacobsen E<lt>dana@acm.orgE<gt>

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut
