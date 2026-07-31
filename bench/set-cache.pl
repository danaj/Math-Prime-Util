#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long qw(GetOptions);
use Time::HiRes qw(time);
use Math::Prime::Util qw(prime_get_config setcontains setcontainsany
                         setinsert setremove);

die "This benchmark requires the XS implementation\n"
  unless prime_get_config()->{xs};

my $seconds = 0.5;
my $max_queries = 400;
my $domain = 'nonnegative';
GetOptions(
  'seconds=f' => \$seconds,
  'queries=i' => \$max_queries,
  'domain=s'  => \$domain,
) or die "Usage: $0 [--seconds N] [--queries N] "
       . "[--domain nonnegative|negative|mixed] [set sizes ...]\n";

die "seconds must be positive\n" if $seconds <= 0;
die "queries must be at least 6\n" if $max_queries < 6;
die "domain must be nonnegative, negative, or mixed\n"
  unless $domain =~ /\A(?:nonnegative|negative|mixed)\z/;

my @sizes = @ARGV ? @ARGV : (30, 500, 8_000, 2_000_000);
for my $size (@sizes) {
  die "set sizes must be integers of at least 21\n"
    if $size !~ /^\d+$/ || $size < 21;
}

my $sink = 0;

sub value_at {
  my ($size, $index) = @_;
  return 2 * $index if $domain eq 'nonnegative';
  return 2 * $index - 2 * $size if $domain eq 'negative';
  return 2 * $index - $size;
}

sub make_set {
  my ($size) = @_;
  my @set;
  for (my $i = 0; $i < $size; $i++) {
    $set[$i] = value_at($size, $i);
  }
  return \@set;
}

sub make_queries {
  my ($size, $count, $clustered) = @_;
  my (@hits, @misses);
  my $start = int(($size - $count) / 2);

  for (my $i = 0; $i < $count; $i++) {
    my $index = $clustered
              ? $start + $i
              : int(($i + 1) * $size / ($count + 1));
    my $value = value_at($size, $index);
    push @hits,   $value;
    push @misses, $value + 1;
  }
  return (\@hits, \@misses);
}

sub measure_rate {
  my ($code, $min_time) = @_;
  my $batch = 1;

  # Keep clock and loop overhead small compared with the operation itself.
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
  } while ($elapsed < $min_time);

  return $calls / $elapsed;
}

sub check_unchanged {
  my ($set, $size, $name) = @_;
  die "$name changed the source set\n"
    if @$set != $size || $set->[0] != value_at($size, 0)
                     || $set->[-1] != value_at($size, $size-1);
}

printf "set cache benchmark: %s, %.2fs per workload, at most %d queries\n",
       $domain, $seconds, $max_queries;
printf "%10s  %7s  %-25s  %12s  %10s  %10s\n",
       qw(set_size queries workload calls/s us/call ns/query);

for my $size (@sizes) {
  my $count = int($size * 0.4);
  $count = $max_queries if $count > $max_queries;
  $count = 6 if $count < 6;

  my $set = make_set($size);
  my ($spread_hits, $spread_misses) = make_queries($size, $count, 0);
  my ($cluster_hits, $cluster_misses) = make_queries($size, $count, 1);

  my @workloads = (
    [ 'contains spread hits',    sub { $sink += setcontains($set, $spread_hits) } ],
    [ 'contains cluster hits',   sub { $sink += setcontains($set, $cluster_hits) } ],
    [ 'containsany spread miss', sub { $sink += setcontainsany($set, $spread_misses) } ],
    [ 'containsany cluster miss',sub { $sink += setcontainsany($set, $cluster_misses) } ],
    [ 'insert existing spread',  sub { $sink += setinsert($set, $spread_hits) } ],
    [ 'insert existing cluster', sub { $sink += setinsert($set, $cluster_hits) } ],
    [ 'remove absent spread',    sub { $sink += setremove($set, $spread_misses) } ],
    [ 'remove absent cluster',   sub { $sink += setremove($set, $cluster_misses) } ],
  );

  for my $workload (@workloads) {
    my ($name, $code) = @$workload;
    $code->();
    check_unchanged($set, $size, $name);
    my $rate = measure_rate($code, $seconds);
    check_unchanged($set, $size, $name);
    printf "%10d  %7d  %-25s  %12.1f  %10.3f  %10.1f\n",
           $size, $count, $name, $rate, 1_000_000 / $rate,
           1_000_000_000 / ($rate * $count);
  }
  print "\n";
}

# Make the benchmarked return values observable.
print "sink: $sink\n" if $ENV{MPU_BENCH_SHOW_SINK};
