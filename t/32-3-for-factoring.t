#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/fordivisors forfactored forsquarefree forsquarefreeint
                         vecsum sqrtint/;

my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;

subtest 'fordivisors' => sub {
  my $a = 0;
  fordivisors { $a += $_ + $_*$_ } 54321;
  is($a, 3287796520, "fordivisors: d|54321: a+=d+d^2");

  my @A027750 = (1,1,2,1,3,1,2,4,1,5,1,2,3,6,1,7,1,2,4,8,1,3,9,1,2,5,10,1,11,1,2,3,4,6,12,1,13,1,2,7,14,1,3,5,15,1,2,4,8,16,1,17,1,2,3,6,9,18,1,19,1,2,4,5,10,20,1,3,7,21,1,2,11,22,1,23,1,2,3,4,6,8,12,24,1,5,25,1,2,13,26,1,3,9,27,1,2,4,7,14,28,1,29,1,2,3,5,6,10,15,30,1,31,1,2,4,8,16,32,1,3,11,33,1,2,17,34,1,5,7,35,1,2,3,4,6,9,12,18,36,1,37,1,2,19,38,1,3,13,39,1,2,4,5,8,10,20,40,1,41,1,2,3,6,7,14,21,42,1,43,1,2,4,11,22,44,1,3,5,9,15,45,1,2,23,46,1,47,1,2,3,4,6,8,12,16,24,48,1,7,49,1,2,5,10,25,50);
  my @a;
  for my $n (1..50) {  fordivisors { push @a, $_ } $n;  }
  is_deeply(\@a, \@A027750, "A027750 using fordivisors");
};

subtest 'forfactored and squarefree iterators' => sub {
  my $s;
  $s=0; forfactored      { $s += 1+$_ } 0,0; is($s, 0, "forfactored {} 0,0");
  $s=0; forsquarefree    { $s += 1+$_ } 0,0; is($s, 0, "forsquarefree {} 0,0");
  $s=0; forsquarefreeint { $s += 1+$_ } 0,0; is($s, 0, "forsquarefreeint {} 0,0");
  $s=0; forfactored      { $s += 1+$_ } 0,1; is($s, 2, "forfactored {} 0,1");
  $s=0; forsquarefree    { $s += 1+$_ } 0,1; is($s, 2, "forsquarefree {} 0,1");
  $s=0; forsquarefreeint { $s += 1+$_ } 0,1; is($s, 2, "forsquarefreeint {} 0,1");

  $s=0; forfactored { $s += $_ } 1; is($s, 1, "forfactored {} 1");
  $s=0; forfactored { $s += vecsum($_,@_) } 100; is($s, 7330, "forfactored {} 100");
  $s=0; forsquarefree { $s += vecsum($_,@_) } 100; is($s, 4763, "forsquarefree {} 100");
  $s=0; forfactored { $s += vecsum($_,@_) } 1e8,1e8+10; is($s, 1208835222, "forfactored {} 10^8,10^8+10");
  is( a053462(6), 607926, "A053462 using forsquarefree");

  $s = 0; forsquarefree    { $s += $_ } 7193953,7195732; is($s, 7813597636, "forsquarefree {} 7193953,7195732");
  $s = 0; forsquarefreeint { $s += $_ } 7193953,7195732; is($s, 7813597636, "forsquarefreeint {} 7193953,7195732");
};

subtest 'bigint ranges' => sub {
  my($E,%d);
  if ($use64) {
    $E = 66;
    $d{sf} = [26,29,27,[7,13,19,223,2683,16981,4200451],29,[3,31,379,"2093425718354419"]];
    $d{sfint} = [26,29,27,29];
    $d{factored} = [29,30,29,[qw/3 31 379 2093425718354419/],30,[qw/2 17 1129 1922236656459079/]];
  } else {
    $E = 36;
    $d{sf} = [305,308,306,[2,79,3617,120247]];
    $d{sfint} = [636,639,637];
    $d{factored} = [170,171,170,[2,3,3,19,89,2257687],171,[17,149,1033,26263]];
  }

  { my @r;  my($a1,$a2,$arg1,$arg2,@res) = split_d($E,$d{sf});
    forsquarefree { push @r,"$_",[map{"$_"}@_] } $arg1, $arg2;
    is_deeply(\@r, \@res, "forsquarefree {} 2^$E+$a1, 2^$E+$a2");
  }
  { my @r;  my($a1,$a2,$arg1,$arg2,@res) = split_d($E,$d{sfint});
    forsquarefreeint { push @r,"$_" } $arg1, $arg2;
    is_deeply(\@r, \@res, "forsquarefreeint {} 2^$E+$a1, 2^$E+$a2");
  }
  { my @r;  my($a1,$a2,$arg1,$arg2,@res) = split_d($E,$d{factored});
    forfactored { push @r,"$_",[map{"$_"}@_] } $arg1, $arg2;
    is_deeply(\@r, \@res, "forfactored {} 2^$E+$a1, 2^$E+$a2");
  }

  my @r;
  fordivisors { push @r,"$_" } "73786976294838225404";
  is_deeply(\@r, [qw/1 2 4 137 274 548 134647766961383623 269295533922767246 538591067845534492 18446744073709556351 36893488147419112702 73786976294838225404/], "fordivisors {} 2^66+18940");
};

sub a053462 {
  my($s,$n)=(0,10**$_[0]-1);
  forsquarefree { $s += int($n / ($_*$_)) * ((scalar(@_) & 1)?-1:1); } sqrtint($n);
  $s;
}

sub split_d {
  my($E,$arr) = @_;
  my $a1 = $arr->[0];
  my $a2 = $arr->[1];
  return ($a1,$a2,map { ref($_) ? $_ : 2**$E + $_    } @$arr) if $E < 48;
  return ($a1,$a2,map { ref($_) ? $_ : plus_2_66($_) } @$arr) if $E == 66;
  die "unsupported test exponent $E";
}

sub plus_2_66 {
  my $add = shift;
  my $final = $add + 6464;
  die "too large" if $final > 99999;
  return "737869762948382" . sprintf("%05d",$final);
}

done_testing();
