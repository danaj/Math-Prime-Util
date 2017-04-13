#!/usr/bin/env perl
use warnings;
use strict;
use v5.16;
use ntheory;

sub outqr {
  my($bits,$a,$b,$c,$d) = @_;

  ($a,$b,$c,$d) = map { length($_) == 1 ? "$_ " : $_ } ($a,$b,$c,$d);

  my $qr32 = <<'EOT';
  {use integer;$a+=$b;} $d^=$a; $d=($d<<16)|($d>>16);
  {use integer;$c+=$d;} $b^=$c; $b=($b<<12)|($b>>20);
  {use integer;$a+=$b;} $d^=$a; $d=($d<< 8)|($d>>24);
  {use integer;$c+=$d;} $b^=$c; $b=($b<< 7)|($b>>25);
EOT
  my $qr64 = <<'EOT';
  $a=($a+$b)&0xFFFFFFFF; $d^=$a; $d=(($d<<16)|($d>>16))&0xFFFFFFFF;
  $c=($c+$d)&0xFFFFFFFF; $b^=$c; $b=(($b<<12)|($b>>20))&0xFFFFFFFF;
  $a=($a+$b)&0xFFFFFFFF; $d^=$a; $d=(($d<< 8)|($d>>24))&0xFFFFFFFF;
  $c=($c+$d)&0xFFFFFFFF; $b^=$c; $b=(($b<< 7)|($b>>25))&0xFFFFFFFF;
EOT

  my $qr = ($bits == 32) ? $qr32 : $qr64;

  $qr =~ s/\$a/\$x$a/g;
  $qr =~ s/\$b/\$x$b/g;
  $qr =~ s/\$c/\$x$c/g;
  $qr =~ s/\$d/\$x$d/g;
  $qr =~ s/^/      /mg;
  $qr =~ s/\n$//;

  say $qr;
}

say "      if (BITS == 64) {";
say "        use integer;";
  outqr(64,0,4,8,12);
  outqr(64,1,5,9,13);
  outqr(64,2,6,10,14);
  outqr(64,3,7,11,15);
  outqr(64,0,5,10,15);
  outqr(64,1,6,11,12);
  outqr(64,2,7,8,13);
  outqr(64,3,4,9,14);
say "      } else { # 32-bit";
  outqr(32,0,4,8,12);
  outqr(32,1,5,9,13);
  outqr(32,2,6,10,14);
  outqr(32,3,7,11,15);
  outqr(32,0,5,10,15);
  outqr(32,1,6,11,12);
  outqr(32,2,7,8,13);
  outqr(32,3,4,9,14);
say "      }";
