#!/usr/bin/env perl
use strict;
use warnings;

$| = 1;

use Math::BigInt;
use Math::Prime::Util qw/
  addint subint negint powint urandomb csrand vecequal
  toset setunion setintersect setminus setdelta
  sumset setbinop setcontains setcontainsany
  setinsert setremove setinvert
  set_is_disjoint set_is_equal
  set_is_subset set_is_proper_subset
  set_is_superset set_is_proper_superset
  set_is_proper_intersection
  is_sidon_set is_sumfree_set
/;
use Math::Prime::Util::PP ();

my @names = qw/
  toset setunion setintersect setminus setdelta
  sumset setbinop setcontains setcontainsany
  setinsert setremove setinvert
  set_is_disjoint set_is_equal
  set_is_subset set_is_proper_subset
  set_is_superset set_is_proper_superset
  set_is_proper_intersection
  is_sidon_set is_sumfree_set
/;
my %M = map { $_ => Math::Prime::Util->can($_) } @names;
my %P = map { $_ => Math::Prime::Util::PP->can($_) } @names;

my $seed = @ARGV ? shift : 4;
die "usage: $0 [seed]\n" if @ARGV;
csrand($seed);

my @widths = (3, 8, 16, 31, 32, 33, 53, 63, 64, 65,
              80, 127, 128, 129, 192);

my ($ncases, $nchecks) = (0, 0);

sub brief {
  my ($values) = @_;
  my @v = map { "$_" } @$values;
  splice(@v, 6, @v-10, '...') if @v > 10;
  return '[' . join(',', @v) . ']';
}

sub fail_case {
  my ($label, $ra, $rb, $detail) = @_;
  die "$label: $detail\nA=" . brief($ra) . "\nB=" . brief($rb) . "\n";
}

sub same_values {
  my ($ra, $rb) = @_;
  return @$ra == @$rb && vecequal($ra, $rb);
}

sub check_values {
  my ($label, $ra, $rb, $got, $expected) = @_;
  $nchecks++;
  fail_case($label, $ra, $rb,
            'got ' . brief($got) . ', expected ' . brief($expected))
    if !same_values($got, $expected);
}

sub check_scalar {
  my ($label, $ra, $rb, $got, $expected) = @_;
  $nchecks++;
  fail_case($label, $ra, $rb, "got $got, expected $expected")
    if !defined($got) || "$got" ne "$expected";
}

sub rotate {
  my ($values, $n) = @_;
  return [@$values] if @$values < 2;
  $n %= @$values;
  return [@$values] if $n == 0;
  return [@$values[$n .. $#$values], @$values[0 .. $n-1]];
}

sub random_values {
  my ($bits, $length, $signed) = @_;
  my @v = map { urandomb($bits) } 1 .. $length;
  if ($signed) {
    my $middle = powint(2, $bits-1);
    @v = map { subint($_, $middle) } @v;
  }
  return \@v;
}

sub series {
  my ($bits, $length, $offset, $signed) = @_;
  my $origin = powint(2, $bits-1);
  $origin = negint($origin) if $signed;
  return [map { addint($origin, $offset+$_) } 0 .. $length-1];
}

sub decorate {
  my ($values, $style) = @_;
  my @v;
  for my $i (0 .. $#$values) {
    my $kind = $style == 3 ? $i % 4 : $style;
    my $value = $values->[$i];
    if ($kind == 0) {
      push @v, $value;
    } elsif ($kind == 1) {
      push @v, "$value";
    } elsif ($kind == 2) {
      push @v, Math::BigInt->new("$value");
    } else {
      push @v, Math::BigInt->new("$value")->bstr;
    }
  }
  return \@v;
}

sub make_case {
  my ($bits, $profile, $na, $nb, $serial) = @_;
  my $signmode = $serial % 4;
  my $asigned = $signmode == 1 || $signmode == 2;
  my $bsigned = $signmode == 1 || $signmode == 3;
  my ($ra, $rb);

  if ($profile eq 'empty') {
    ($ra, $rb) = ([], []);
  } elsif ($profile eq 'random') {
    $ra = random_values($bits, $na, $asigned);
    $rb = random_values($bits, $nb, $bsigned);
  } elsif ($profile eq 'equal') {
    $ra = series($bits, $na, -int($na/2), $asigned);
    $rb = rotate($ra, $serial);
    push @$ra, $ra->[0], $ra->[-1];
    push @$rb, $rb->[0];
  } elsif ($profile eq 'b_subset_a') {
    $ra = series($bits, $na, -int($na/2), $asigned);
    my $step = int($na/$nb);
    $rb = [map { $ra->[$_*$step] } 0 .. $nb-1];
    $ra = rotate($ra, $serial);
  } elsif ($profile eq 'a_subset_b') {
    $rb = series($bits, $nb, -int($nb/2), $bsigned);
    my $step = int($nb/$na);
    $ra = [map { $rb->[$_*$step] } 0 .. $na-1];
    $rb = rotate($rb, $serial);
  } elsif ($profile eq 'disjoint') {
    my $nmax = $na > $nb ? $na : $nb;
    my $all = series($bits, 2*$nmax, 0, $asigned);
    $ra = [map { $all->[2*$_] } 0 .. $na-1];
    $rb = [map { $all->[2*$_+1] } 0 .. $nb-1];
  } elsif ($profile eq 'overlap') {
    $ra = series($bits, $na, 0, $asigned);
    $rb = series($bits, $nb, int($na/2), $asigned);
  } elsif ($profile eq 'late_miss') {
    $ra = series($bits, $na, 0, $asigned);
    $rb = [@$ra];
    $rb->[-1] = addint($rb->[-1], 2);
  } elsif ($profile eq 'duplicates') {
    my $poolsize = $na > $nb ? $na : $nb;
    $poolsize = 13 if $poolsize > 13;
    $poolsize = 1 if $poolsize < 1;
    my $pool = series($bits, $poolsize, -int($poolsize/2), $asigned);
    $ra = [map { $pool->[($_*7+$serial) % $poolsize] } 0 .. $na-1];
    $rb = [map { $pool->[($_*5+$serial+3) % $poolsize] } 0 .. $nb-1];
  } elsif ($profile eq 'sorted_nonnegative') {
    $ra = series($bits, $na, 0, 0);
    $rb = series($bits, $nb, int($na/2), 0);
    return ($ra, $rb);
  } else {
    die "unknown profile $profile\n";
  }

  $ra = rotate($ra, $serial) if $profile ne 'late_miss';
  $rb = rotate($rb, $serial+3) if $profile ne 'late_miss';
  return (decorate($ra, $serial % 4),
          decorate($rb, int($serial/3) % 4));
}

sub check_case {
  my ($label, $ra, $rb, $serial) = @_;

  my $seta = $M{toset}->(@$ra);
  my $setb = $M{toset}->(@$rb);
  check_values("$label toset A", $ra, $rb, $seta, $P{toset}->(@$ra));
  check_values("$label toset B", $ra, $rb, $setb, $P{toset}->(@$rb));

  my %result;
  for my $name (qw/setunion setintersect setminus setdelta/) {
    $result{$name} = $M{$name}->($ra, $rb);
    check_values("$label $name", $ra, $rb, $result{$name},
                 $P{$name}->($ra, $rb));
  }

  check_values("$label delta identity", $ra, $rb, $result{setdelta},
               setunion(setminus($ra,$rb), setminus($rb,$ra)));
  check_values("$label union partition", $ra, $rb, $result{setunion},
               setunion($result{setintersect}, $result{setdelta}));

  for my $name (qw/
    set_is_equal set_is_disjoint
    set_is_subset set_is_proper_subset
    set_is_superset set_is_proper_superset
    set_is_proper_intersection
  /) {
    check_scalar("$label $name", $ra, $rb,
                 $M{$name}->($seta,$setb), $P{$name}->($seta,$setb));
  }

  my $subset = $M{set_is_subset}->($seta,$setb);
  my $disjoint = $M{set_is_disjoint}->($seta,$setb);
  for my $impl ([MPU => \%M], [PP => \%P]) {
    my ($iname, $functions) = @$impl;
    check_scalar("$label $iname setcontains ref", $ra, $rb,
                 $functions->{setcontains}->($seta,$setb), $subset);
    check_scalar("$label $iname setcontains list", $ra, $rb,
                 $functions->{setcontains}->($seta,reverse(@$rb)), $subset);
    check_scalar("$label $iname setcontainsany ref", $ra, $rb,
                 $functions->{setcontainsany}->($seta,$setb),
                 $disjoint ? 0 : 1);
    check_scalar("$label $iname setcontainsany list", $ra, $rb,
                 $functions->{setcontainsany}->($seta,reverse(@$rb)),
                 $disjoint ? 0 : 1);
  }

  my %mutator_result = (
    setinsert => $result{setunion},
    setremove => $result{setminus},
    setinvert => $result{setdelta},
  );
  for my $impl ([MPU => \%M], [PP => \%P]) {
    my ($iname, $functions) = @$impl;
    for my $name (qw/setinsert setremove setinvert/) {
      my $a = [@$seta];
      my $oldsize = @$a;
      my $got = $serial % 2
              ? $functions->{$name}->($a,$setb)
              : $functions->{$name}->($a,reverse(@$rb));
      my $wanted = $mutator_result{$name};
      my $wanted_count = $name eq 'setremove'
                       ? $oldsize - @$wanted
                       : @$wanted - $oldsize;
      check_scalar("$label $iname $name count", $ra, $rb,
                   $got, $wanted_count);
      check_values("$label $iname $name result", $ra, $rb, $a, $wanted);
    }
  }

  if (@$seta * @$setb <= 100 && $serial % 3 == 0) {
    my $sum = $M{sumset}->($ra,$rb);
    check_values("$label sumset", $ra, $rb, $sum, $P{sumset}->($ra,$rb));
    my $madd = $M{setbinop}->(sub { addint($a,$b) },$ra,$rb);
    my $padd = $P{setbinop}->(sub { addint($a,$b) },$ra,$rb);
    check_values("$label setbinop add", $ra, $rb, $madd, $padd);
    check_values("$label sumset identity", $ra, $rb, $sum, $madd);
    check_values(
      "$label setbinop subtract", $ra, $rb,
      $M{setbinop}->(sub { subint($a,$b) },$ra,$rb),
      $P{setbinop}->(sub { subint($a,$b) },$ra,$rb)
    );
  }

  if (@$seta <= 24) {
    for my $name (qw/is_sidon_set is_sumfree_set/) {
      check_scalar("$label $name", $ra, $rb,
                   $M{$name}->($seta), $P{$name}->($seta));
    }
  }
}

print "seed $seed\n";
# Exhaust every pair of input sizes through 30.  Cycle the data shape,
# signedness, and scalar representation instead of repeating one random shape.
for my $wi (0 .. $#widths) {
  my $bits = $widths[$wi];
  print "bits $bits:";
  for my $na (0 .. 30) {
    for my $nb (0 .. 30) {
      my $serial = ($wi*31 + $na)*31 + $nb;
      my $mode = $serial % 8;
      my $profile =
          $na == 0 && $nb == 0                 ? 'empty'
        : $mode == 2                           ? 'duplicates'
        : $mode == 3                           ? 'overlap'
        : $mode == 4                           ? 'disjoint'
        : $mode == 5 && $na >= $nb && $nb > 0 ? 'b_subset_a'
        : $mode == 5 && $nb >  $na && $na > 0 ? 'a_subset_b'
        : $mode == 6 && $na == $nb && $na > 0 ? 'late_miss'
        : $mode == 7 && $na == $nb && $na > 0 ? 'equal'
        :                                         'random';
      my ($ra, $rb) = make_case($bits,$profile,$na,$nb,$serial);
      check_case("bits=$bits profile=$profile sizes=$na/$nb",
                 $ra,$rb,$serial);
      $ncases++;
    }
    print '.';
  }
  print "\n";
}

# These sizes straddle implementation changes in XS and PP.  Testing one
# representative case on each side is more useful than sweeping large sets.
my @transitions = (
  [xs_insert_4       => 31,    40,   4, 'disjoint'           ],
  [xs_insert_5       => 31,    40,   5, 'disjoint'           ],
  [xs_first_20       => 31,    20,   6, 'b_subset_a'         ],
  [xs_first_21       => 31,    21,   6, 'b_subset_a'         ],
  [pp_contains_150_2 => 31,   150,   2, 'b_subset_a'         ],
  [pp_contains_151_2 => 31,   151,   2, 'b_subset_a'         ],
  [pp_contains_250_3 => 31,   250,   3, 'b_subset_a'         ],
  [pp_contains_251_3 => 31,   251,   3, 'b_subset_a'         ],
  [pp_remove_34_3    => 31,    34,   3, 'b_subset_a'         ],
  [pp_remove_35_3    => 31,    35,   3, 'b_subset_a'         ],
  [pp_invert_100     => 31,   200, 100, 'b_subset_a'         ],
  [pp_invert_101     => 31,   200, 101, 'b_subset_a'         ],
  [pp_insert_400     => 31,   450, 400, 'disjoint'           ],
  [pp_insert_401     => 31,   450, 401, 'disjoint'           ],
  [xs_remove_499     => 31,   600, 499, 'b_subset_a'         ],
  [xs_remove_500     => 31,   600, 500, 'b_subset_a'         ],
  [xs_fast_merge_31  => 31,    30,  25, 'sorted_nonnegative' ],
  [xs_fast_merge_63  => 63,    30,  25, 'sorted_nonnegative' ],
  [bigint_transition => 65,    30,  25, 'sorted_nonnegative' ],
);

print "transitions:";
for my $i (0 .. $#transitions) {
  my ($name,$bits,$na,$nb,$profile) = @{$transitions[$i]};
  my $serial = 100000 + $i;
  my ($ra,$rb) = make_case($bits,$profile,$na,$nb,$serial);
  check_case("$name bits=$bits sizes=$na/$nb", $ra,$rb,$serial);
  $ncases++;
  print '.';
}
print "\n";
print "PASS: $ncases diverse cases, $nchecks checks\n";
