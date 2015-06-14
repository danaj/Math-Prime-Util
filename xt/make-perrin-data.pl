#!/usr/bin/env perl
use warnings;
use strict;
use ntheory ":all";
use Math::GMPz;

# https://oeis.org/A104217/b104217.txt
my %mods;
open(my $pfile, '<', '/home/dana/Downloads/b104217.txt') or die "Cannot open b104217.txt\n";
while (<$pfile>) {
  next unless /^(\d+)\s+(\d+)/;
  $mods{$1} = $2;
}
close($pfile) or die "Error on close\n";


my @maskdata;
my @struct;
my $offset = 0;
for my $mod (sort {$a<=>$b} keys %mods) {
  next unless $mod >= 2;
  next unless is_prime($mod) || (is_power($mod,2) && is_prime(sqrtint($mod)));
  my $period = $mods{$mod};
  my $pwords = int(($period+31)/32);
  next if $period > 65535 || $offset > 65535;
  #last if $mod > 59;
  next unless $mod * $pwords < 3000;

  my @P = map { Math::GMPz->new($_) } (3,0,2);
  my @nums = (0) x $pwords;
  for (0 .. $period-1) {
    my $v = $P[0] % $mod;
    if ($v == 0) {
      $nums[int($_/32)] |= 1 << ($_ % 32);
    }
    @P=($P[1],$P[2],$P[0]+$P[1]);
  }

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
