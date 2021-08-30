#!/usr/bin/env perl
use warnings;
use strict;
use ntheory ":all";
use Math::Prime::Util::PP;


show(7,3);
#show(13,4);
#show(13,13);

sub show {
  my($D,$d) = @_;
  printf "Name         %2d/%2d    %2d/%2d    %2d/%2d    %2d/%2d\n",
    $D,$d, $D,-$d,  -$D,$d,  -$D,-$d;
  printf "---------    -----    -----    -----    -----\n";
  showd($D,$d,"Perl naive", \&ndivrem);
  showd($D,$d,"Perl floor", \&nfdivrem);
  print "\n";
  showd($D,$d,"Perl ui",    \&uidivrem);
  showd($D,$d,"tdivrem",    \&tdivrem);
  showd($D,$d,"PP tdivrem", \&Math::Prime::Util::PP::tdivrem);
  print "\n";
  showd($D,$d,"fdivrem",    \&fdivrem);
  showd($D,$d,"div/mod",    \&div_and_mod);
  showd($D,$d,"PP fdivrem", \&Math::Prime::Util::PP::fdivrem);
  print "\n";
  showd($D,$d,"divrem",     \&divrem);
  showd($D,$d,"PP divrem",  \&Math::Prime::Util::PP::divrem);
}

sub showd {
  my($D,$d,$name,$func) = @_;
  #printf "%-8s  %2d/%2d = %2d %2d     %2d/%2d = %2d %2d\n", "$name:", -$D,$d, $func->(-$D,$d), $D,-$d, $func->($D,-$d);
  printf "%-11s  %2d %2d    %2d %2d    %2d %2d    %2d %2d\n",
    $name, $func->($D,$d), $func->($D,-$d), $func->(-$D,$d), $func->(-$D,-$d);
}

sub ndivrem {
  my($D,$d) = @_;
  ( int($D/$d), $D % $d );
}
sub uidivrem {
  my($D,$d) = @_;
  use integer;
  ( int($D/$d), $D % $d );
}
sub nfdivrem {
  my($D,$d) = @_;
  use POSIX;
  ( POSIX::floor($D/$d), $D % $d );
}
sub div_and_mod {
  my($D,$d) = @_;
  ( divint($D,$d), modint($D,$d) );
}
sub ivmod {
  my($a,$n) = @_;
  die "wrong usage: n must be positive" unless $n >= 0;
  return $a % $n if $a >= 0;
  my $amodn = -$a % $n;
  return ($amodn == 0) ? 0 : $n-$amodn;
}

