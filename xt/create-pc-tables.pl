#!/usr/bin/env perl
use strict;
use warnings;
use ntheory ":all";
use v5.20;

my $s = 0;
if(0) {
  $s+= make_table(1,    15,   0,   9);
  $s+= make_table(2,    30,   9,  39);
  $s+= make_table(3,    60,  39,  63);
  $s+= make_table(4,   300,  60, 120);
  $s+= make_table(5, 30000, 120,3000);
}
if(0) {
  $s+= make_table(1,    15,   0,   9);
  $s+= make_table(2,    30,   9,  39);
  $s+= make_table(3,    60,  39,  63);
  $s+= make_table(4,    60,  63,  90);
  $s+= make_table(5, 15000,  90,3000);
}
if(0) {
  $s+= make_table(1,     5,   0,   1.5);
  $s+= make_table(2,    15,   1.5, 12);
  $s+= make_table(3,    30,  12,  39);
  $s+= make_table(4,    30,  39,  66);
  $s+= make_table(5,    60,  66,  90);
  $s+= make_table(6, 30000,  90,3000);
}
if(1) {
  #                      k    M       M
  $s+= make_table(0,     3,   0,      0.30);
  $s+= make_table(1,     6,   0.30,   3.0 );
  $s+= make_table(2,    15,   3.0,   15   );
  $s+= make_table(3,    30,  15,     42   );
  $s+= make_table(4,    30,  42,     69   );
  $s+= make_table(5,    60,  69,     90   );
  $s+= make_table(6, 30000,  90,   3000   );
}
say "/* $s bytes */";


sub make_table {
  my($name, $stepk, $start, $stop) = @_;
  my $step = 1000 * $stepk;
  $start *= 1_000_000;
  $stop  *= 1_000_000;

  die "start must be less than stop" unless $start < $stop;
  die "start must be divisible by step" unless ($start % $step) == 0;
  die "stop  must be divisible by step" unless ($stop  % $step) == 0;
  my $s = $start / $step;
  my $pc = prime_count($start);
  my $nsteps = ($stop - $start) / $step;

  if ($start == 0) {
    $s = 0;
    $pc = prime_count(5);
  }

  my @c;
  {
    my($npc,$spc) = ($pc);
    @c = map { ($spc,$npc) = ($npc, prime_count(($s+$_)*$step)); $npc-$spc; } 1 .. $nsteps;
  }
  my $min = vecmin(@c);
  @c = map { $_-$min } @c;
  my $max = vecmax(@c);

  say "#define NSTEP_STEP_$name   $step";
  say "#define NSTEP_START_$name  $start";
  say "#define NSTEP_COUNT_$name  $pc";
  say "#define NSTEP_BASE_$name   $min";
  my $type = ($max <= 255) ? "char" : ($max <= 65535) ? "short" : "int";
  say "static const unsigned $type step_counts_${name}[] =";
  say "{",join(",",@c),"};";
  say "#define NSTEP_NUM_$name  (sizeof(step_counts_$name)/sizeof(step_counts_${name}[0]))";
  say "";
  return scalar(@c) * (($max <= 255) ? 1 : ($max <= 65535) ? 2 : 4);
}
