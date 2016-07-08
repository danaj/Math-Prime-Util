#!/usr/bin/env perl
use warnings;
use strict;
use ntheory ":all";
use Math::GMPz;

# https://oeis.org/A104217/b104217.txt
my %mods;
open(my $pfile, '<', 'b104217.txt') or die "Cannot open b104217.txt\n";
while (<$pfile>) {
  next unless /^(\d+)\s+(\d+)/;
  $mods{$1} = $2;
}
close($pfile) or die "Error on close\n";


my @maskdata;
my @struct;
my $offset = 0;
for my $mod (sort {$a<=>$b} keys %mods) {
  last if $offset > 65535;
  my $period = $mods{$mod};
  next if $mod < 2 || $period > 65535;
  next unless is_prime($mod) || (is_power($mod,2) && is_prime(sqrtint($mod)));

  # Find the zeros
  my @P = (3,0,2);
  my @zeros;
  for (0 .. $period-1) {
    push @zeros, $_ if ($P[0] % $mod) == 0;
    @P = ($P[1], $P[2], ($P[0]+$P[1]) % $mod);
  }
  my $nzeros = scalar(@zeros);

  my $pwords = int(($period+31)/32);
  next unless $pwords < 5000;
  my @nums = (0) x $pwords;
  for (@zeros) {
    $nums[int($_/32)] |= 1 << ($_ % 32);
  }

  my $bytesperzero = $pwords*4 / $nzeros;
  my $expect = (1/$mod) * $nzeros;
  next unless $expect > 0.003;
  next unless $bytesperzero < 100;
  #print "mod $mod  nzeros $nzeros  bpz $bytesperzero  exp $expect\n";

  push @struct, "  {$mod, $period, $offset}";
  push @maskdata, @nums;
  $offset += scalar(@nums);
}
print "#define NPERRINDIV ", scalar(@struct), "\n";
print "/* ", 4*scalar(@maskdata), " mask bytes */\n";
print "static const uint32_t _perrinmask[] = {",
      join(",", map { ($_ > 2147483647) ? "${_}U" : $_ } @maskdata),
      "};\n";
print "static _perrin _perrindata[NPERRINDIV] = {\n",
      join(",\n", @struct),
      "\n};\n";
