#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long qw(GetOptions);
use Time::HiRes qw(time);

# Compare natural APIs for a cyclic workload on one 2048-bit positive value.
# Each round performs 56 operations and restores the original value.  Parsing,
# mask construction, module loading, and correctness checks are not timed.
#
# Math::BigInt backends are measured in separate processes because a process
# can load only one backend.  Optional modules are reported as skipped.
#
#   perl -Mblib bench/bitops-2048.pl
#   perl -Mblib bench/bitops-2048.pl --seconds 2

my $WIDTH = 2048;
my $POSITION_COUNT = 8;
my $OPS_PER_ROUND = 56;
my $seconds = 0.75;
my $worker = '';

GetOptions(
  'seconds=f' => \$seconds,
  'worker=s'  => \$worker,
) or die "Usage: $0 [--seconds N]\n";
die "seconds must be positive\n" unless $seconds > 0;

sub bits_to_hex {
  my ($bits) = @_;
  my $hex = '';
  for (my $high = $WIDTH; $high > 0; $high -= 4) {
    my $low = $high - 4;
    my $v = 0;
    for my $j (0 .. 3) {
      $v |= 1 << $j if $bits->[$low+$j];
    }
    $hex .= sprintf('%x', $v);
  }
  $hex =~ s/^0+//;
  return length($hex) ? $hex : '0';
}

sub bits_to_vec {
  my ($bits) = @_;
  my $v = "\0" x (($WIDTH + 7) >> 3);
  for my $k (0 .. $WIDTH-1) {
    vec($v, $k, 1) = 1 if $bits->[$k];
  }
  return $v;
}

sub pick_positions {
  my ($bits, $wanted, $offset) = @_;
  my @positions;
  my $k = $offset;
  while (@positions < $POSITION_COUNT) {
    push @positions, $k if $bits->[$k] == $wanted;
    $k = ($k + 251) % ($WIDTH-1);
  }
  return \@positions;
}

sub make_dataset {
  my (@initial, @and, @or, @xor, @andnot);

  for my $k (0 .. $WIDTH-1) {
    $initial[$k] = (($k * 73 + int($k / 8) + 19) % 29) < 14 ? 1 : 0;
    $and[$k]     = (($k * 37 + int($k / 5) +  7) % 31) < 23 ? 1 : 0;
    $or[$k]      = (($k * 43 + int($k / 9) + 11) % 37) < 15 ? 1 : 0;
    $xor[$k]     = (($k * 53 + int($k / 7) +  3) % 41) < 19 ? 1 : 0;
    $andnot[$k]  = (($k * 61 + int($k / 6) + 17) % 43) < 17 ? 1 : 0;
  }

  # Every intermediate remains a 2048-bit positive value.
  $initial[$WIDTH-1] = 1;
  $and[$WIDTH-1] = 1;
  $xor[$WIDTH-1] = 0;
  $andnot[$WIDTH-1] = 0;

  my @and_restore    = map { $initial[$_] && !$and[$_] ? 1 : 0 }
                           0 .. $WIDTH-1;
  my @or_keep        = map { $initial[$_] || !$or[$_] ? 1 : 0 }
                           0 .. $WIDTH-1;
  my @andnot_keep    = map { !$andnot[$_] ? 1 : 0 } 0 .. $WIDTH-1;
  my @andnot_restore = map { $initial[$_] && $andnot[$_] ? 1 : 0 }
                           0 .. $WIDTH-1;

  my $setpos   = pick_positions(\@initial, 0, 17);
  my $clearpos = pick_positions(\@initial, 1, 89);
  my @flippos  = map { (137 + 337 * $_) % ($WIDTH-1) }
                     0 .. $POSITION_COUNT-1;

  my %bits = (
    initial        => \@initial,
    and            => \@and,
    and_restore    => \@and_restore,
    or             => \@or,
    or_keep        => \@or_keep,
    xor            => \@xor,
    andnot         => \@andnot,
    andnot_keep    => \@andnot_keep,
    andnot_restore => \@andnot_restore,
  );
  my %hex = map { $_ => bits_to_hex($bits{$_}) } keys %bits;

  my (%single_hex, %clear_hex);
  my %seen;
  for my $k (@$setpos, @$clearpos, @flippos) {
    next if $seen{$k}++;
    my @single = (0) x $WIDTH;
    my @clear  = (1) x $WIDTH;
    $single[$k] = 1;
    $clear[$k] = 0;
    $single_hex{$k} = bits_to_hex(\@single);
    $clear_hex{$k} = bits_to_hex(\@clear);
  }

  return {
    bits       => \%bits,
    hex        => \%hex,
    single_hex => \%single_hex,
    clear_hex  => \%clear_hex,
    setpos     => $setpos,
    clearpos   => $clearpos,
    flippos    => \@flippos,
  };
}

sub load_module {
  my ($module) = @_;
  (my $file = $module) =~ s{::}{/}g;
  $file .= '.pm';
  return eval { require $file; 1 };
}

sub object_masks {
  my ($data, $make) = @_;
  my %dense = map { $_ => $make->($data->{hex}{$_}) }
                    keys %{$data->{hex}};
  my %single = map { $_ => $make->($data->{single_hex}{$_}) }
                     keys %{$data->{single_hex}};
  my %clear = map { $_ => $make->($data->{clear_hex}{$_}) }
                    keys %{$data->{clear_hex}};
  return (\%dense, \%single, \%clear);
}

sub build_mpu {
  my ($data, $class, $module, $label) = @_;

  if ($class eq 'Math::BigInt') {
    return (undef, $@) unless eval {
      require Math::BigInt;
      Math::BigInt->import(only => 'Calc');
      1;
    };
  } else {
    return (undef, "$module is not installed") unless load_module($module);
  }
  return (undef, 'Math::Prime::Util is not available')
    unless load_module('Math::Prime::Util');
  Math::Prime::Util::prime_set_config(bigint => $class);

  my $make = sub { Math::Prime::Util::toint('0x' . $_[0]) };
  my ($d, $single, $clear) = object_masks($data, $make);
  my $x = $make->($data->{hex}{initial});
  my $expected = "$x";
  my @setpos = @{$data->{setpos}};
  my @clearpos = @{$data->{clearpos}};
  my @flippos = @{$data->{flippos}};

  my $run = sub {
    for my $k (@setpos)   { $x = Math::Prime::Util::bitset($x, $k) }
    for my $k (@setpos)   { $x = Math::Prime::Util::bitclear($x, $k) }
    for my $k (@clearpos) { $x = Math::Prime::Util::bitclear($x, $k) }
    for my $k (@clearpos) { $x = Math::Prime::Util::bitset($x, $k) }
    for my $k (@flippos)  { $x = Math::Prime::Util::bitflip($x, $k) }
    for my $k (@flippos)  { $x = Math::Prime::Util::bitflip($x, $k) }
    $x = Math::Prime::Util::bitand($x, $d->{and});
    $x = Math::Prime::Util::bitor($x, $d->{and_restore});
    $x = Math::Prime::Util::bitor($x, $d->{or});
    $x = Math::Prime::Util::bitand($x, $d->{or_keep});
    $x = Math::Prime::Util::bitxor($x, $d->{xor});
    $x = Math::Prime::Util::bitxor($x, $d->{xor});
    $x = Math::Prime::Util::bitandnot($x, $d->{andnot});
    $x = Math::Prime::Util::bitor($x, $d->{andnot_restore});
  };
  my $check = sub { "$x" eq $expected };
  return ({label => $label, model => 'immutable', run => $run,
           check => $check, detail => 'output ' . (ref($x) || 'scalar')}, '');
}

sub build_bigint {
  my ($data, $backend) = @_;
  return (undef, "Math::BigInt::$backend is not installed") unless eval {
    require Math::BigInt;
    Math::BigInt->import(only => $backend);
    1;
  };
  my $lib = Math::BigInt->config()->{lib};
  my $make = sub { Math::BigInt->from_hex('0x' . $_[0]) };
  my ($d, $single, $clear) = object_masks($data, $make);
  my $x = $make->($data->{hex}{initial});
  my $expected = "$x";
  my @setpos = @{$data->{setpos}};
  my @clearpos = @{$data->{clearpos}};
  my @flippos = @{$data->{flippos}};

  my $run = sub {
    for my $k (@setpos)   { $x->bior($single->{$k}) }
    for my $k (@setpos)   { $x->band($clear->{$k}) }
    for my $k (@clearpos) { $x->band($clear->{$k}) }
    for my $k (@clearpos) { $x->bior($single->{$k}) }
    for my $k (@flippos)  { $x->bxor($single->{$k}) }
    for my $k (@flippos)  { $x->bxor($single->{$k}) }
    $x->band($d->{and});
    $x->bior($d->{and_restore});
    $x->bior($d->{or});
    $x->band($d->{or_keep});
    $x->bxor($d->{xor});
    $x->bxor($d->{xor});
    $x->band($d->{andnot_keep});
    $x->bior($d->{andnot_restore});
  };
  my $check = sub { "$x" eq $expected };
  return ({label => "Math::BigInt/$backend", model => 'mutable', run => $run,
           check => $check, detail => $lib}, '');
}

sub build_math_gmp {
  my ($data) = @_;
  return (undef, 'Math::GMP is not installed') unless load_module('Math::GMP');
  my $make = sub { Math::GMP->new('0x' . $_[0]) };
  my ($d, $single, $clear) = object_masks($data, $make);
  my $x = $make->($data->{hex}{initial});
  my $expected = "$x";
  my @setpos = @{$data->{setpos}};
  my @clearpos = @{$data->{clearpos}};
  my @flippos = @{$data->{flippos}};

  my $run = sub {
    for my $k (@setpos)   { $x = $x->bior($single->{$k}, 0) }
    for my $k (@setpos)   { $x = $x->band($clear->{$k}, 0) }
    for my $k (@clearpos) { $x = $x->band($clear->{$k}, 0) }
    for my $k (@clearpos) { $x = $x->bior($single->{$k}, 0) }
    for my $k (@flippos)  { $x = $x->bxor($single->{$k}, 0) }
    for my $k (@flippos)  { $x = $x->bxor($single->{$k}, 0) }
    $x = $x->band($d->{and}, 0);
    $x = $x->bior($d->{and_restore}, 0);
    $x = $x->bior($d->{or}, 0);
    $x = $x->band($d->{or_keep}, 0);
    $x = $x->bxor($d->{xor}, 0);
    $x = $x->bxor($d->{xor}, 0);
    $x = $x->band($d->{andnot_keep}, 0);
    $x = $x->bior($d->{andnot_restore}, 0);
  };
  my $check = sub { "$x" eq $expected };
  return ({label => 'Math::GMP', model => 'immutable', run => $run,
           check => $check, detail => "version $Math::GMP::VERSION"}, '');
}

sub build_math_gmpz {
  my ($data) = @_;
  return (undef, 'Math::GMPz is not installed') unless load_module('Math::GMPz');
  my $make = sub { Math::GMPz->new('0x' . $_[0]) };
  my ($d, $single, $clear) = object_masks($data, $make);
  my $x = $make->($data->{hex}{initial});
  my $expected = "$x";
  my @setpos = @{$data->{setpos}};
  my @clearpos = @{$data->{clearpos}};
  my @flippos = @{$data->{flippos}};

  my $run = sub {
    for my $k (@setpos)   { Math::GMPz::Rmpz_setbit($x, $k) }
    for my $k (@setpos)   { Math::GMPz::Rmpz_clrbit($x, $k) }
    for my $k (@clearpos) { Math::GMPz::Rmpz_clrbit($x, $k) }
    for my $k (@clearpos) { Math::GMPz::Rmpz_setbit($x, $k) }
    for my $k (@flippos)  { Math::GMPz::Rmpz_combit($x, $k) }
    for my $k (@flippos)  { Math::GMPz::Rmpz_combit($x, $k) }
    Math::GMPz::Rmpz_and($x, $x, $d->{and});
    Math::GMPz::Rmpz_ior($x, $x, $d->{and_restore});
    Math::GMPz::Rmpz_ior($x, $x, $d->{or});
    Math::GMPz::Rmpz_and($x, $x, $d->{or_keep});
    Math::GMPz::Rmpz_xor($x, $x, $d->{xor});
    Math::GMPz::Rmpz_xor($x, $x, $d->{xor});
    Math::GMPz::Rmpz_and($x, $x, $d->{andnot_keep});
    Math::GMPz::Rmpz_ior($x, $x, $d->{andnot_restore});
  };
  my $check = sub { "$x" eq $expected };
  return ({label => 'Math::GMPz', model => 'mutable', run => $run,
           check => $check, detail => "version $Math::GMPz::VERSION"}, '');
}

sub build_mpu_gmp {
  my ($data) = @_;
  return (undef, 'Math::Prime::Util::GMP is not installed')
    unless load_module('Math::Prime::Util::GMP');
  my $make = sub {
    my @digits = map { hex($_) } split //, $_[0];
    return Math::Prime::Util::GMP::fromdigits(\@digits, 16);
  };
  my ($d, $single, $clear) = object_masks($data, $make);
  my $x = $make->($data->{hex}{initial});
  my $expected = "$x";
  my @setpos = @{$data->{setpos}};
  my @clearpos = @{$data->{clearpos}};
  my @flippos = @{$data->{flippos}};

  my $run = sub {
    for my $k (@setpos)   { $x = Math::Prime::Util::GMP::setbit($x, $k) }
    for my $k (@setpos)   { $x = Math::Prime::Util::GMP::clrbit($x, $k) }
    for my $k (@clearpos) { $x = Math::Prime::Util::GMP::clrbit($x, $k) }
    for my $k (@clearpos) { $x = Math::Prime::Util::GMP::setbit($x, $k) }
    for my $k (@flippos)  { $x = Math::Prime::Util::GMP::notbit($x, $k) }
    for my $k (@flippos)  { $x = Math::Prime::Util::GMP::notbit($x, $k) }
    $x = Math::Prime::Util::GMP::bitand($x, $d->{and});
    $x = Math::Prime::Util::GMP::bitor($x, $d->{and_restore});
    $x = Math::Prime::Util::GMP::bitor($x, $d->{or});
    $x = Math::Prime::Util::GMP::bitand($x, $d->{or_keep});
    $x = Math::Prime::Util::GMP::bitxor($x, $d->{xor});
    $x = Math::Prime::Util::GMP::bitxor($x, $d->{xor});
    $x = Math::Prime::Util::GMP::bitand($x, $d->{andnot_keep});
    $x = Math::Prime::Util::GMP::bitor($x, $d->{andnot_restore});
  };
  my $check = sub { "$x" eq $expected };
  return ({label => 'MPU::GMP', model => 'immutable', run => $run,
           check => $check,
           detail => "version $Math::Prime::Util::GMP::VERSION"}, '');
}

sub build_vec {
  my ($data) = @_;
  my %d = map { $_ => bits_to_vec($data->{bits}{$_}) }
                keys %{$data->{bits}};
  my $x = $d{initial};
  my $expected = $x;
  my @setpos = @{$data->{setpos}};
  my @clearpos = @{$data->{clearpos}};
  my @flippos = @{$data->{flippos}};

  my $run = sub {
    for my $k (@setpos)   { vec($x, $k, 1) = 1 }
    for my $k (@setpos)   { vec($x, $k, 1) = 0 }
    for my $k (@clearpos) { vec($x, $k, 1) = 0 }
    for my $k (@clearpos) { vec($x, $k, 1) = 1 }
    for my $k (@flippos)  { vec($x, $k, 1) ^= 1 }
    for my $k (@flippos)  { vec($x, $k, 1) ^= 1 }
    $x &= $d{and};
    $x |= $d{and_restore};
    $x |= $d{or};
    $x &= $d{or_keep};
    $x ^= $d{xor};
    $x ^= $d{xor};
    $x &= $d{andnot_keep};
    $x |= $d{andnot_restore};
  };
  my $check = sub { $x eq $expected };
  return ({label => 'Perl vec()', model => 'mutable', run => $run,
           check => $check, detail => 'packed 256-byte string'}, '');
}

sub build_bit_vector {
  my ($data) = @_;
  return (undef, 'Bit::Vector is not installed') unless load_module('Bit::Vector');
  my $make = sub {
    my ($bits) = @_;
    my $v = Bit::Vector->new($WIDTH);
    for my $k (0 .. $WIDTH-1) {
      $v->Bit_On($k) if $bits->[$k];
    }
    return $v;
  };
  my %d = map { $_ => $make->($data->{bits}{$_}) }
                keys %{$data->{bits}};
  my $x = $d{initial}->Clone();
  my $expected = $x->to_Hex();
  my @setpos = @{$data->{setpos}};
  my @clearpos = @{$data->{clearpos}};
  my @flippos = @{$data->{flippos}};

  my $run = sub {
    for my $k (@setpos)   { $x->Bit_On($k) }
    for my $k (@setpos)   { $x->Bit_Off($k) }
    for my $k (@clearpos) { $x->Bit_Off($k) }
    for my $k (@clearpos) { $x->Bit_On($k) }
    for my $k (@flippos)  { $x->bit_flip($k) }
    for my $k (@flippos)  { $x->bit_flip($k) }
    $x->And($x, $d{and});
    $x->Or($x, $d{and_restore});
    $x->Or($x, $d{or});
    $x->And($x, $d{or_keep});
    $x->Xor($x, $d{xor});
    $x->Xor($x, $d{xor});
    $x->AndNot($x, $d{andnot});
    $x->Or($x, $d{andnot_restore});
  };
  my $check = sub { $x->to_Hex() eq $expected };
  return ({label => 'Bit::Vector', model => 'mutable', run => $run,
           check => $check, detail => "version $Bit::Vector::VERSION"}, '');
}

sub build_worker {
  my ($id, $data) = @_;
  return build_mpu($data, 'Math::BigInt', 'Math::BigInt',
                   'MPU (Math::BigInt/Calc)') if $id eq 'mpu_calc';
  return build_mpu($data, 'Math::GMPz', 'Math::GMPz',
                   'MPU (Math::GMPz)') if $id eq 'mpu_gmpz';
  return build_mpu($data, 'Math::GMP', 'Math::GMP',
                   'MPU (Math::GMP)') if $id eq 'mpu_math_gmp';
  return build_bigint($data, 'Calc') if $id eq 'bigint_calc';
  return build_bigint($data, 'GMP') if $id eq 'bigint_gmp';
  return build_math_gmp($data) if $id eq 'math_gmp';
  return build_math_gmpz($data) if $id eq 'math_gmpz';
  return build_mpu_gmp($data) if $id eq 'mpu_gmp';
  return build_vec($data) if $id eq 'vec';
  return build_bit_vector($data) if $id eq 'bit_vector';
  return (undef, "unknown worker '$id'");
}

sub measure_rate {
  my ($code, $minimum) = @_;
  my $batch = 1;

  while (1) {
    my $start = time;
    my $left = $batch;
    $code->() while $left--;
    my $elapsed = time - $start;
    last if $elapsed >= 0.02;
    my $scale = $elapsed > 0 ? int(0.02 / $elapsed) + 1 : 10;
    $scale = 10 if $scale > 10;
    $batch *= $scale;
  }

  my $calls = 0;
  my $start = time;
  my $elapsed;
  do {
    my $left = $batch;
    $code->() while $left--;
    $calls += $batch;
    $elapsed = time - $start;
  } while ($elapsed < $minimum);
  return $calls / $elapsed;
}

sub worker_main {
  my ($id) = @_;
  my $data = make_dataset();
  my ($bench, $reason) = build_worker($id, $data);
  if (!$bench) {
    $reason =~ s/[\r\n\t]+/ /g;
    print "MPU_BITOPS_SKIP\t$id\t$reason\n";
    return;
  }

  $bench->{run}->();
  die "$bench->{label} did not restore its input\n" unless $bench->{check}->();
  my $rate = measure_rate($bench->{run}, $seconds);
  die "$bench->{label} did not restore its input\n" unless $bench->{check}->();
  print join("\t", 'MPU_BITOPS_RESULT', $id, $bench->{label},
             $bench->{model}, sprintf('%.12g', $rate), $bench->{detail}), "\n";
}

if (length $worker) {
  worker_main($worker);
  exit 0;
}

my @workers = qw(
  mpu_calc mpu_gmpz mpu_math_gmp bigint_calc bigint_gmp
  math_gmp math_gmpz mpu_gmp vec bit_vector
);
my (@results, @skips);
for my $id (@workers) {
  my @inc = map { "-I$_" } grep { !ref($_) && length($_) } @INC;
  open(my $fh, '-|', $^X, @inc, $0, "--worker=$id", "--seconds=$seconds")
    or die "Cannot start $id worker: $!\n";
  while (my $line = <$fh>) {
    chomp $line;
    if ($line =~ /^MPU_BITOPS_RESULT\t/) {
      my (undef, $wid, $label, $model, $rate, $detail) = split /\t/, $line, 6;
      push @results, {id => $wid, label => $label, model => $model,
                      rate => 0 + $rate, detail => $detail};
    } elsif ($line =~ /^MPU_BITOPS_SKIP\t/) {
      my (undef, $wid, $reason) = split /\t/, $line, 3;
      push @skips, "$wid: $reason";
    } else {
      print "$id: $line\n";
    }
  }
  close($fh) or die "$id worker failed\n";
}

die "No benchmark implementation is available\n" unless @results;
@results = sort { $a->{rate} <=> $b->{rate} } @results;
my $slowest = $results[0]{rate};

printf "2048-bit cyclic bit-operation benchmark\n";
printf "%d operations/round, %.2fs per implementation\n\n",
       $OPS_PER_ROUND, $seconds;
printf "%-27s  %-9s  %12s  %10s  %9s\n",
       'implementation', 'model', 'rounds/s', 'ns/op', 'relative';
for my $r (@results) {
  printf "%-27s  %-9s  %12.1f  %10.1f  %8.2fx\n",
         $r->{label}, $r->{model}, $r->{rate},
         1.0e9 / ($r->{rate} * $OPS_PER_ROUND), $r->{rate} / $slowest;
}

if (@skips) {
  print "\nSkipped:\n";
  print "  $_\n" for @skips;
}

print "\nSet, clear, flip, and AND-NOT use prebuilt masks where an\n";
print "implementation has no corresponding direct operation.  Details:\n";
for my $r (@results) {
  print "  $r->{label}: $r->{detail}\n";
}
