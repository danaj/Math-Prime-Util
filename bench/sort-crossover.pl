#!perl
use strict;
use warnings;
use Benchmark qw(cmpthese);
use Config;
use Time::HiRes qw(time);
use Math::Prime::Util qw(
  csrand irand32 irand64 subint vecsorti
);

# Every array is sorted exactly once.  Working copies are prepared outside
# the timed region, so the timings contain neither input generation nor array
# copying.  Useful controls:
#
#   MPU_SORT_TYPES=sint64,uint64
#   MPU_SORT_PATTERNS=random,organ_pipe
#   MPU_SORT_VALUES=1000000
#   MPU_SORT_SAMPLES=7
#   perl -Mblib bench/sort-crossover.pl 384 399 400 432

my @sizes = @ARGV ? @ARGV : (256, 384, 399, 400, 512, 640, 768, 896, 1024);
my @types = split /,/, ($ENV{MPU_SORT_TYPES} ||
                        ($Config{uvsize} >= 8
                          ? 'sint32,uint32,sint64,uint64'
                          : 'sint32,uint32'));
my @patterns = split /,/, ($ENV{MPU_SORT_PATTERNS} ||
                           'random,ascending,descending,organ_pipe,sawtooth');
my $target_values = $ENV{MPU_SORT_VALUES} || 500_000;
my $samples = $ENV{MPU_SORT_SAMPLES} || 5;

my %valid_type = map { $_ => 1 } qw(sint32 uint32 sint64 uint64);
my %valid_pattern = map { $_ => 1 }
                    qw(random ascending descending organ_pipe sawtooth);

die "Sizes must be integers of at least 2\n"
  if grep { !/^\d+\z/ || $_ < 2 } @sizes;
die "Unknown input type\n" if grep { !$valid_type{$_} } @types;
die "Unknown input pattern\n" if grep { !$valid_pattern{$_} } @patterns;
die "MPU_SORT_VALUES and MPU_SORT_SAMPLES must be positive\n"
  if $target_values < 1 || $samples < 1;
die "The 64-bit categories require a 64-bit Perl\n"
  if $Config{uvsize} < 8 && grep { /64\z/ } @types;

my $have_djb = eval { require Sort::DJB; 1 };
my $have_sort_key = eval { require Sort::Key; 1 };
my $have_radix = eval { require Sort::Key::Radix; 1 };
my $have_packed = eval { require Sort::Packed; 1 };

my %type_info = (
  sint32 => {
    djb    => 'sort_int32',
    djb_label => 'int32',
    skey   => 'isort',
    radix  => 'isort',
    packed => ['l', 'Packed_int32'],
  },
  uint32 => {
    djb    => 'sort_uint32',
    djb_label => 'uint32',
    skey   => 'usort',
    radix  => 'usort',
    packed => ['L', 'Packed_u32'],
  },
  sint64 => {
    djb    => 'sort_int64',
    djb_label => 'int64',
    skey   => 'isort',
    radix  => 'isort',
    packed => ['q', 'Packed_i64'],
  },
  uint64 => {
    djb    => 'sort_uint64',
    djb_label => 'uint64',
    skey   => 'usort',
    radix  => 'usort',
    packed => ['Q', 'Packed_u64'],
  },
);

sub random_value {
  my ($type) = @_;
  return irand32() if $type eq 'uint32';
  return subint(irand32(),  '2147483648') if $type eq 'sint32';
  return irand64() if $type eq 'uint64';
  return subint(irand64(), '9223372036854775808');
}

sub numeric_sort {
  return sort { $a <=> $b } @_;
}

sub make_input {
  my ($type, $pattern, $size) = @_;
  my @v = map { random_value($type) } 1 .. $size;

  return \@v if $pattern eq 'random';

  if ($pattern eq 'ascending') {
    @v = numeric_sort(@v);
    return \@v;
  }
  if ($pattern eq 'descending') {
    @v = reverse numeric_sort(@v);
    return \@v;
  }
  if ($pattern eq 'organ_pipe') {
    my $nhalf = int(($size + 1) / 2);
    my @half = numeric_sort(@v[0 .. $nhalf-1]);
    my @tail = reverse @half[0 .. $size-$nhalf-1];
    return [@half, @tail];
  }

  # A duplicate-heavy sawtooth with 32 distinct full-width values.
  my $nteeth = $size < 32 ? $size : 32;
  my @teeth = numeric_sort(@v[0 .. $nteeth-1]);
  return [map { $teeth[$_ % $nteeth] } 0 .. $size-1];
}

sub same_array {
  my ($x, $y) = @_;
  return 0 unless @$x == @$y;
  for my $i (0 .. $#$x) {
    return 0 unless "$x->[$i]" eq "$y->[$i]";
  }
  return 1;
}

sub median {
  my @v = sort { $a <=> $b } @_;
  my $n = @v;
  return $v[int($n/2)] if $n & 1;
  return ($v[$n/2-1] + $v[$n/2]) / 2;
}

sub time_sorter {
  my ($code, $source) = @_;
  my @elapsed;

  {
    my $warm = [@{$source->[0]}];
    $code->($warm);
  }

  for (1 .. $samples) {
    my @work = map { [@$_] } @$source;
    my $start = time;
    $code->($_) for @work;
    push @elapsed, time - $start;
  }

  my $elapsed = median(@elapsed);
  return bless [$elapsed, $elapsed, 0, 0, 0, scalar(@$source)], 'Benchmark';
}

sub packed_supported {
  my ($format) = @_;
  return 0 unless $have_packed;
  return eval {
    my $buf = pack("$format*", 2, 1);
    Sort::Packed::sort_packed($format, $buf);
    my @v = unpack("$format*", $buf);
    @v == 2 && $v[0] == 1 && $v[1] == 2;
  };
}

sub sorters_for {
  my ($type) = @_;
  my $info = $type_info{$type};
  my @sorters;

  if ($have_djb) {
    my $name = "DJB_$info->{djb_label}";
    my $code = Sort::DJB->can($info->{djb});
    push @sorters, [$name, sub { $code->($_[0]) }] if $code;
  }

  push @sorters, ['Perl_sort', sub {
    my ($v) = @_;
    return [sort { $a <=> $b } @$v];
  }];

  if ($have_sort_key) {
    my $method = $info->{skey};
    my $code = Sort::Key->can($method);
    push @sorters, ["SKey_$method", sub {
      my ($v) = @_;
      return [$code->(@$v)];
    }] if $code;
  }

  if ($have_radix) {
    my $method = $info->{radix};
    my $code = Sort::Key::Radix->can($method);
    push @sorters, ["Radix_$method", sub {
      my ($v) = @_;
      return [$code->(@$v)];
    }] if $code;
  }

  my ($format, $name) = @{$info->{packed}};
  if (packed_supported($format)) {
    push @sorters, [$name, sub {
      my ($v) = @_;
      my $buf = pack("$format*", @$v);
      Sort::Packed::sort_packed($format, $buf);
      return [unpack("$format*", $buf)];
    }];
  }

  push @sorters, ['MPU_vecsorti', sub { vecsorti($_[0]) }];
  return @sorters;
}

print "=" x 72, "\n";
print "Math::Prime::Util integer sort crossover benchmark\n";
print "=" x 72, "\n\n";
printf "Perl                 : %s\n", $^V;
printf "UV size              : %d bits\n", 8*$Config{uvsize};
printf "Sort::DJB            : %s\n",
       $have_djb ? "v$Sort::DJB::VERSION" : 'not installed';
printf "Sort::Key            : %s\n",
       $have_sort_key ? "v$Sort::Key::VERSION" : 'not installed';
printf "Sort::Key::Radix     : %s\n",
       $have_radix ? "v$Sort::Key::Radix::VERSION" : 'not installed';
printf "Sort::Packed         : %s\n",
       $have_packed ? "v$Sort::Packed::VERSION" : 'not installed';
printf "Timing samples       : %d (median reported)\n", $samples;
printf "Values per sample    : approximately %d\n", $target_values;

for my $type (@types) {
  my @sorters = sorters_for($type);

  for my $pattern (@patterns) {
    print "\n", "-" x 72, "\n";
    print "BENCHMARK: $type / $pattern\n";
    print "-" x 72, "\n\n";

    for my $size (@sizes) {
      my $narrays = int($target_values / $size);
      $narrays = 16 if $narrays < 16;

      # Per-case seeding keeps filtered, complete, old, and new runs directly
      # comparable.
      csrand("bench/sort-crossover.pl:$type:$pattern:$size");
      my @source = map { make_input($type, $pattern, $size) }
                   1 .. $narrays;
      my @expected = numeric_sort(@{$source[0]});
      my %results;

      printf "  n = %4d elements (%d independent arrays):\n",
             $size, $narrays;

      for my $sorter (@sorters) {
        my ($name, $code) = @$sorter;
        my $check = [@{$source[0]}];
        $check = $code->($check);
        die "$name failed for $type/$pattern/$size\n"
          unless same_array($check, \@expected);
        $results{$name} = time_sorter($code, \@source);
      }
      cmpthese(\%results);
      print "\n";
    }
  }
}
