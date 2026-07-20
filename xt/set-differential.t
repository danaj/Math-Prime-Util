#!/usr/bin/env perl
use strict;
use warnings;

use Config;
use Math::BigInt lib => 'Calc';
use Test::More;
use Math::Prime::Util qw/
  prime_get_config addint subint
  toset setunion setintersect setminus setdelta
  sumset setbinop
  setcontains setcontainsany
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
  sumset setbinop
  setcontains setcontainsany
  setinsert setremove setinvert
  set_is_disjoint set_is_equal
  set_is_subset set_is_proper_subset
  set_is_superset set_is_proper_superset
  set_is_proper_intersection
  is_sidon_set is_sumfree_set
/;

my %mpu = map { $_ => Math::Prime::Util->can($_) } @names;
my %pp  = map { $_ => Math::Prime::Util::PP->can($_) } @names;
my @impl = ([MPU => \%mpu], [PP => \%pp]);

diag('set implementations: ',
     prime_get_config()->{xs} ? 'XS and PP' : 'PP-only build');

my %number;

sub number {
  my ($s) = @_;
  $number{$s} = Math::BigInt->new($s) if !exists $number{$s};
  return $number{$s};
}

sub canon {
  my ($v) = @_;
  my $s = Math::BigInt->new("$v")->bstr;
  number($s);
  return $s;
}

sub sorted_values {
  return [sort { number($a)->bcmp(number($b)) } @_];
}

sub model_set {
  my %seen;
  my @set;
  for my $v (@_) {
    my $s = canon($v);
    push @set, $s if !$seen{$s}++;
  }
  return sorted_values(@set);
}

sub sequence {
  my ($a) = @_;
  return [map { canon($_) } @$a];
}

sub sequence_key {
  my ($a) = @_;
  return scalar(@$a) . ':' .
         join("\x1E", map { length($_) . ':' . $_ } @$a);
}

sub same_sequence {
  my ($a, $b) = @_;
  return sequence_key($a) eq sequence_key($b);
}

sub show_sequence {
  my ($a) = @_;
  return '[' . join(',', @$a) . ']';
}

sub model_op {
  my ($op, $ra, $rb) = @_;
  my $a = model_set(@$ra);
  my $b = model_set(@$rb);
  my (%a, %b, @r);
  @a{@$a} = (1) x @$a;
  @b{@$b} = (1) x @$b;

  if ($op eq 'union') {
    @r = keys %{{%a, %b}};
  } elsif ($op eq 'intersect') {
    @r = grep { $b{$_} } @$a;
  } elsif ($op eq 'minus') {
    @r = grep { !$b{$_} } @$a;
  } elsif ($op eq 'delta') {
    @r = ((grep { !$b{$_} } @$a), (grep { !$a{$_} } @$b));
  } else {
    die "unknown model set operation $op";
  }
  return sorted_values(@r);
}

sub model_add {
  my ($a, $b) = map { canon($_) } @_;
  return number($a)->copy->badd(number($b))->bstr;
}

sub model_sub {
  my ($a, $b) = map { canon($_) } @_;
  return number($a)->copy->bsub(number($b))->bstr;
}

sub model_binop {
  my ($op, $ra, $rb) = @_;
  my @r;
  for my $a (@{model_set(@$ra)}) {
    for my $b (@{model_set(@$rb)}) {
      push @r, $op eq 'add' ? model_add($a,$b) : model_sub($a,$b);
    }
  }
  return model_set(@r);
}

sub invoke {
  my ($fail, $checks, $label, $code) = @_;
  my ($ok, $ret, $error);
  $$checks++;
  {
    local $@;
    $ok = eval {
      $ret = $code->();
      1;
    };
    $error = $@;
  }
  if (!$ok) {
    $error =~ s/\s+\z//;
    push @$fail, "$label died: $error";
    return (0, undef);
  }
  return (1, $ret);
}

sub check_sequence_call {
  my ($fail, $checks, $label, $expected, $code) = @_;
  my ($ok, $got) = invoke($fail, $checks, $label, $code);
  return if !$ok;
  if (ref($got) ne 'ARRAY') {
    push @$fail, "$label returned " . (ref($got) || 'a scalar') .
                 ', expected an array reference';
    return;
  }
  my $cgot = sequence($got);
  push @$fail, "$label got " . show_sequence($cgot) .
               ', expected ' . show_sequence($expected)
    if !same_sequence($cgot, $expected);
}

sub check_scalar_call {
  my ($fail, $checks, $label, $expected, $code) = @_;
  my ($ok, $got) = invoke($fail, $checks, $label, $code);
  return if !$ok;
  my $sgot = defined($got) ? "$got" : '<undef>';
  push @$fail, "$label got $sgot, expected $expected"
    if $sgot ne "$expected";
}

sub check_preserved {
  my ($fail, $checks, $label, $got, $expected) = @_;
  $$checks++;
  my $cgot = sequence($got);
  push @$fail, "$label changed input to " . show_sequence($cgot) .
               ', expected ' . show_sequence($expected)
    if !same_sequence($cgot, $expected);
}

sub finish_checks {
  my ($description, $checks, $fail) = @_;
  ok(!@$fail, "$description ($checks checks)");
  if (@$fail) {
    my $last = $#$fail < 9 ? $#$fail : 9;
    diag($_) for @$fail[0 .. $last];
    diag('additional failures: ' . (@$fail-10)) if @$fail > 10;
  }
}

my $bits = 8 * $Config{uvsize};
my $two = Math::BigInt->new(2);
my $ivmin = $two->copy->bpow($bits-1)->bneg->bstr;
my $ivmax = $two->copy->bpow($bits-1)->bdec->bstr;
my $uvmax = $two->copy->bpow($bits)->bdec->bstr;
my $uvplus1 = Math::BigInt->new($uvmax)->binc->bstr;
my $widepos = $two->copy->bpow(130)->badd(12345)->bstr;
my $wideneg = Math::BigInt->new($widepos)->bneg->bstr;

my @cases = (
  [ [], [], 'empty' ],
  [ [1], [1], 'equal singleton' ],
  [ [-9,-7,-5], [2,4,6], 'disjoint signed' ],
  [ [-5 .. 5], [-3 .. 3], 'dense subset' ],
  [ [0 .. 40], [20 .. 60], 'dense overlap' ],
  [ [map { 101*$_ } 0 .. 20], [map { 103*$_ } 0 .. 18],
    'sparse overlap' ],
  [ [0 .. 199], [0 .. 198, 201], 'late mismatch' ],
  [ [$ivmin,-1,0,1,$ivmax,$uvmax],
    [$ivmin,0,$ivmax,$uvmax], 'native boundaries' ],
  [ [$wideneg,$ivmin,0,$uvplus1,$widepos],
    [$wideneg,-1,$uvplus1], 'wide integers' ],
  [ [Math::BigInt->new($wideneg), Math::BigInt->new(0),
     Math::BigInt->new($widepos)],
    [Math::BigInt->new($widepos), Math::BigInt->new(7)],
    'Math::BigInt objects' ],
);

for my $round (0 .. 47) {
  my $na = ($round * 7 + 3) % 23;
  my $nb = ($round * 11 + 5) % 21;
  my @a = map {
    (($round+3) * ($_+5) * 17 + $_*$_*11) % 161 - 80
  } 0 .. $na-1;
  my @b = map {
    (($round+7) * ($_+2) * 23 + $_*$_*13) % 173 - 86
  } 0 .. $nb-1;
  push @a, $a[0] if @a && $round % 3 == 0;
  push @b, $b[-1] if @b && $round % 4 == 0;

  if ($round % 8 == 0) {
    @b = @a;
  } elsif ($round % 8 == 1) {
    @b = grep { $_ % 2 } @a;
  } elsif ($round % 8 == 2) {
    @a = map { $_-500 } @a;
    @b = map { $_+500 } @b;
  }
  push @cases, [\@a, \@b, "generated $round"];
}

subtest 'pure operations and XS/PP differential' => sub {
  my (@fail, $checks);
  my %opname = (
    setunion     => 'union',
    setintersect => 'intersect',
    setminus     => 'minus',
    setdelta     => 'delta',
  );

  for my $case (@cases) {
    my ($ra, $rb, $name) = @$case;
    my $sa = sequence($ra);
    my $sb = sequence($rb);

    for my $impl (@impl) {
      my ($label, $func) = @$impl;

      check_sequence_call(\@fail, \$checks, "$label toset A: $name",
                          model_set(@$ra),
                          sub { $func->{toset}->(@$ra) });

      for my $fn (qw/setunion setintersect setminus setdelta/) {
        my $a = [@$ra];
        my $b = [@$rb];
        check_sequence_call(
          \@fail, \$checks, "$label $fn: $name",
          model_op($opname{$fn}, $ra, $rb),
          sub { $func->{$fn}->($a,$b) }
        );
        check_preserved(\@fail, \$checks, "$label $fn first input: $name",
                        $a, $sa);
        check_preserved(\@fail, \$checks, "$label $fn second input: $name",
                        $b, $sb);
      }

      if (@$ra * @$rb <= 600) {
        my $add_expected = model_binop('add', $ra, $rb);
        {
          my $ca = [@$ra];
          my $cb = [@$rb];
          check_sequence_call(\@fail, \$checks, "$label sumset: $name",
                              $add_expected,
                              sub { $func->{sumset}->($ca,$cb) });
          check_preserved(\@fail, \$checks,
                          "$label sumset first input: $name", $ca, $sa);
          check_preserved(\@fail, \$checks,
                          "$label sumset second input: $name", $cb, $sb);
        }
        {
          my $ca = [@$ra];
          my $cb = [@$rb];
          check_sequence_call(
            \@fail, \$checks, "$label setbinop add: $name", $add_expected,
            sub {
              $func->{setbinop}->(
                sub { Math::Prime::Util::addint($a,$b) }, $ca, $cb
              );
            }
          );
          check_preserved(\@fail, \$checks,
                          "$label setbinop add first input: $name", $ca, $sa);
          check_preserved(\@fail, \$checks,
                          "$label setbinop add second input: $name", $cb, $sb);
        }
        {
          my $ca = [@$ra];
          my $cb = [@$rb];
          check_sequence_call(
            \@fail, \$checks, "$label setbinop subtract: $name",
            model_binop('sub', $ra, $rb),
            sub {
              $func->{setbinop}->(
                sub { Math::Prime::Util::subint($a,$b) }, $ca, $cb
              );
            }
          );
          check_preserved(
            \@fail, \$checks,
            "$label setbinop subtract first input: $name", $ca, $sa
          );
          check_preserved(
            \@fail, \$checks,
            "$label setbinop subtract second input: $name", $cb, $sb
          );
        }
      }
    }
  }
  finish_checks('pure set operations match independent models',
                $checks, \@fail);
  done_testing();
};

subtest 'same-reference and one-argument forms' => sub {
  my (@fail, $checks);
  my @same_cases = (
    [ [], 'empty' ],
    [ [-8,-3,0,4,11], 'signed native' ],
    [ [$ivmin,0,$ivmax,$uvmax], 'native boundaries' ],
    [ [$wideneg,0,$uvplus1,$widepos], 'wide integers' ],
  );
  my %expected_op = (
    setunion     => 'same',
    setintersect => 'same',
    setminus     => 'empty',
    setdelta     => 'empty',
  );

  for my $case (@same_cases) {
    my ($raw,$name) = @$case;
    my $set = model_set(@$raw);
    for my $impl (@impl) {
      my ($label,$func) = @$impl;

      for my $fn (qw/setunion setintersect setminus setdelta/) {
        my $same = [@$set];
        my $expected = $expected_op{$fn} eq 'same' ? $set : [];
        check_sequence_call(
          \@fail,\$checks,"$label $fn same reference: $name",$expected,
          sub { $func->{$fn}->($same,$same) }
        );
        check_preserved(\@fail,\$checks,
                        "$label $fn same input: $name",$same,$set);
      }

      {
        my $same = [@$set];
        my $expected = model_binop('add',$set,$set);
        check_sequence_call(
          \@fail,\$checks,"$label sumset one argument: $name",$expected,
          sub { $func->{sumset}->($same) }
        );
        check_preserved(\@fail,\$checks,
                        "$label sumset one input: $name",$same,$set);
      }
      {
        my $same = [@$set];
        my $expected = model_binop('add',$set,$set);
        check_sequence_call(
          \@fail,\$checks,"$label setbinop one argument: $name",$expected,
          sub {
            $func->{setbinop}->(
              sub { Math::Prime::Util::addint($a,$b) }, $same
            );
          }
        );
        check_preserved(\@fail,\$checks,
                        "$label setbinop one input: $name",$same,$set);
      }
    }
  }
  finish_checks('same-reference and one-argument forms are correct',
                $checks,\@fail);
  done_testing();
};

subtest 'algebraic identities' => sub {
  my (@fail, $checks);

  for my $case (@cases) {
    my ($ra, $rb, $name) = @$case;
    my ($uok,$union) = invoke(\@fail,\$checks,"union identity: $name",
                              sub { setunion($ra,$rb) });
    my ($iok,$inter) = invoke(\@fail,\$checks,"intersection identity: $name",
                              sub { setintersect($ra,$rb) });
    my ($dok,$delta) = invoke(\@fail,\$checks,"delta identity: $name",
                              sub { setdelta($ra,$rb) });
    my ($aok,$ab) = invoke(\@fail,\$checks,"A-B identity: $name",
                           sub { setminus($ra,$rb) });
    my ($bok,$ba) = invoke(\@fail,\$checks,"B-A identity: $name",
                           sub { setminus($rb,$ra) });
    next if !$uok || !$iok || !$dok || !$aok || !$bok;

    check_sequence_call(
      \@fail, \$checks, "delta is union of directional differences: $name",
      sequence($delta), sub { setunion($ab,$ba) }
    );
    check_sequence_call(
      \@fail, \$checks, "union partitions into intersection and delta: $name",
      sequence($union), sub { setunion($inter,$delta) }
    );
    check_sequence_call(
      \@fail, \$checks, "intersection and delta are disjoint: $name",
      [], sub { setintersect($inter,$delta) }
    );
    check_sequence_call(
      \@fail, \$checks, "union is commutative: $name",
      sequence($union), sub { setunion($rb,$ra) }
    );
    check_sequence_call(
      \@fail, \$checks, "intersection is commutative: $name",
      sequence($inter), sub { setintersect($rb,$ra) }
    );
    check_sequence_call(
      \@fail, \$checks, "delta is commutative: $name",
      sequence($delta), sub { setdelta($rb,$ra) }
    );
  }
  finish_checks('set algebra identities hold', $checks, \@fail);
  done_testing();
};

my @relation_cases = map {
  my ($ra,$rb,$name) = @$_;
  my $a = model_set(@$ra);
  my $b = model_set(@$rb);
  $a = [reverse @$a] if @$a > 1;
  $b = [@$b[1 .. $#$b], $b->[0]] if @$b > 1;
  [$a,$b,$name];
} @cases;

push @relation_cases,
  [ [0 .. 1000], [0,499,997,1000], 'large sparse subset' ],
  [ [0 .. 1000], [0,499,997,1001], 'large late subset miss' ],
  [ [map { 3*$_ } 0 .. 1000], [map { 3*$_ } 900 .. 1000],
    'large tail subset' ],
  [ [map { 3*$_ } 0 .. 1000], [map { 3*$_+1 } 900 .. 1000],
    'large disjoint tail' ];

sub relation_model {
  my ($ra, $rb) = @_;
  my $a = model_set(@$ra);
  my $b = model_set(@$rb);
  my (%a, %b);
  @a{@$a} = (1) x @$a;
  @b{@$b} = (1) x @$b;
  my $inter = scalar grep { $b{$_} } @$a;
  my $b_in_a = !grep { !$a{$_} } @$b;
  my $a_in_b = !grep { !$b{$_} } @$a;
  return {
    set_is_equal               => @$a == @$b && $inter == @$a,
    set_is_disjoint            => $inter == 0,
    set_is_subset              => $b_in_a,
    set_is_proper_subset       => $b_in_a && @$a > @$b,
    set_is_superset            => $a_in_b,
    set_is_proper_superset     => $a_in_b && @$b > @$a,
    set_is_proper_intersection => $inter > 0 &&
                                  $inter < @$a && $inter < @$b,
    setcontains                => $b_in_a,
    setcontainsany             => $inter > 0,
  };
}

subtest 'relations and containment' => sub {
  my (@fail, $checks);
  my @relations = qw/
    set_is_equal set_is_disjoint
    set_is_subset set_is_proper_subset
    set_is_superset set_is_proper_superset
    set_is_proper_intersection
  /;

  for my $case (@relation_cases) {
    my ($rawa,$rawb,$name) = @$case;
    my $model = relation_model($rawa,$rawb);
    my $seta = model_set(@$rawa);
    my $setb = model_set(@$rawb);

    for my $impl (@impl) {
      my ($label,$func) = @$impl;
      for my $fn (@relations) {
        my $a = [@$rawa];
        my $b = [@$rawb];
        check_scalar_call(\@fail,\$checks,"$label $fn: $name",
                          $model->{$fn} ? 1 : 0,
                          sub { $func->{$fn}->($a,$b) });
      }

      check_scalar_call(\@fail,\$checks,"$label setcontains ref: $name",
                        $model->{setcontains} ? 1 : 0,
                        sub { $func->{setcontains}->($seta,$setb) });
      check_scalar_call(\@fail,\$checks,"$label setcontains list: $name",
                        $model->{setcontains} ? 1 : 0,
                        sub {
                          $func->{setcontains}->(
                            $seta, reverse(@$setb), @$setb ? $setb->[0] : ()
                          );
                        });
      check_scalar_call(\@fail,\$checks,"$label setcontainsany ref: $name",
                        $model->{setcontainsany} ? 1 : 0,
                        sub { $func->{setcontainsany}->($seta,$setb) });
      check_scalar_call(\@fail,\$checks,"$label setcontainsany list: $name",
                        $model->{setcontainsany} ? 1 : 0,
                        sub {
                          $func->{setcontainsany}->(
                            $seta, reverse(@$setb), @$setb ? $setb->[0] : ()
                          );
                        });
    }
  }

  for my $impl (@impl) {
    my ($label,$func) = @$impl;
    check_scalar_call(\@fail,\$checks,"$label setcontains empty list",
                      1, sub { $func->{setcontains}->([1,2,3]) });
    check_scalar_call(\@fail,\$checks,"$label setcontainsany empty list",
                      0, sub { $func->{setcontainsany}->([1,2,3]) });
  }

  finish_checks('relations and containment match independent models',
                $checks, \@fail);
  done_testing();
};

my @mutation_cases = (
  [ [], [], 'empty' ],
  [ [1,3,5], [0,2,4,6], 'interleaved' ],
  [ [-20 .. 20], [-30,-20,-10,0,10,20,30], 'signed overlap' ],
  [ [$ivmin,0,$ivmax,$uvmax], [$ivmin,1,$uvmax], 'native boundaries' ],
  [ [$wideneg,0,$widepos], [$wideneg,7,$uvplus1], 'wide integers' ],
  [ [map { 2*$_ } 0 .. 699], [map { 3*$_ } 0 .. 549],
    'large merge paths' ],
);

subtest 'mutators agree with pure operations' => sub {
  my (@fail, $checks);

  for my $case (@mutation_cases) {
    my ($rawa,$rawb,$name) = @$case;
    my $a0 = model_set(@$rawa);
    my $b0 = model_set(@$rawb);
    my @listb = (reverse(@$b0), @$b0 ? ($b0->[0],$b0->[-1]) : ());
    my %expected = (
      setinsert => model_op('union', $a0, $b0),
      setremove => model_op('minus', $a0, $b0),
      setinvert => model_op('delta', $a0, $b0),
    );

    for my $impl (@impl) {
      my ($label,$func) = @$impl;
      for my $fn (qw/setinsert setremove setinvert/) {
        for my $form (qw/reference list/) {
          my $a = [@$a0];
          my $b = [@$b0];
          my $oldsize = @$a;
          my ($ok,$ret) = invoke(
            \@fail,\$checks,"$label $fn $form: $name",
            sub {
              return $form eq 'reference'
                   ? $func->{$fn}->($a,$b)
                   : $func->{$fn}->($a,@listb);
            }
          );
          next if !$ok;

          my $exp = $expected{$fn};
          my $expected_ret = $fn eq 'setremove'
                           ? $oldsize - @$exp
                           : @$exp - $oldsize;
          push @fail, "$label $fn $form return: $name got $ret, " .
                      "expected $expected_ret"
            if "$ret" ne "$expected_ret";
          my $ca = sequence($a);
          push @fail, "$label $fn $form result: $name got " .
                      show_sequence($ca) . ', expected ' .
                      show_sequence($exp)
            if !same_sequence($ca,$exp);
          check_preserved(\@fail,\$checks,
                          "$label $fn second input: $name", $b, $b0)
            if $form eq 'reference';
        }
      }
    }
  }

  for my $impl (@impl) {
    my ($label,$func) = @$impl;
    my $base = model_set(-5,-1,0,4,9);

    my $insert = [@$base];
    check_scalar_call(\@fail,\$checks,"$label same-ref setinsert",0,
                      sub { $func->{setinsert}->($insert,$insert) });
    check_preserved(\@fail,\$checks,"$label same-ref setinsert result",
                    $insert,$base);

    my $remove = [@$base];
    check_scalar_call(\@fail,\$checks,"$label same-ref setremove",
                      scalar(@$base),
                      sub { $func->{setremove}->($remove,$remove) });
    check_preserved(\@fail,\$checks,"$label same-ref setremove result",
                    $remove,[]);

    my $invert = [@$base];
    check_scalar_call(\@fail,\$checks,"$label same-ref setinvert",
                      -scalar(@$base),
                      sub { $func->{setinvert}->($invert,$invert) });
    check_preserved(\@fail,\$checks,"$label same-ref setinvert result",
                    $invert,[]);
  }

  finish_checks('mutators match union, difference, and symmetric difference',
                $checks, \@fail);
  done_testing();
};

sub model_sidon {
  my ($ra) = @_;
  my $set = model_set(@$ra);
  return 0 if @$set && number($set->[0])->is_neg;
  my %sum;
  for my $i (0 .. $#$set) {
    for my $j ($i .. $#$set) {
      my $v = model_add($set->[$i],$set->[$j]);
      return 0 if $sum{$v}++;
    }
  }
  return 1;
}

sub model_sumfree {
  my ($ra) = @_;
  my $set = model_set(@$ra);
  my %in;
  @in{@$set} = (1) x @$set;
  for my $i (0 .. $#$set) {
    for my $j ($i .. $#$set) {
      return 0 if $in{model_add($set->[$i],$set->[$j])};
    }
  }
  return 1;
}

my @property_cases = (
  [], [0], [-1], [1], [1,1,3], [0,1,3], [0,1,4,6],
  [2,7,12], [-5,8,-13,5],
  [$uvmax], [$uvmax,$uvplus1], [$widepos],
);
for my $round (0 .. 47) {
  my $len = $round % 8;
  my @v = map { (($round+5)*($_+3)*7 + $_*$_) % 31 - 10 }
          0 .. $len-1;
  push @v, $v[0] if @v && $round % 5 == 0;
  push @property_cases, \@v;
}

subtest 'Sidon and sum-free predicates' => sub {
  my (@fail, $checks);
  for my $i (0 .. $#property_cases) {
    my $set = $property_cases[$i];
    my $sidon = model_sidon($set);
    my $sumfree = model_sumfree($set);
    for my $impl (@impl) {
      my ($label,$func) = @$impl;
      check_scalar_call(\@fail,\$checks,"$label is_sidon_set case $i",
                        $sidon, sub { $func->{is_sidon_set}->($set) });
      check_scalar_call(\@fail,\$checks,"$label is_sumfree_set case $i",
                        $sumfree, sub { $func->{is_sumfree_set}->($set) });
    }
  }
  finish_checks('Sidon and sum-free predicates match independent models',
                $checks, \@fail);
  done_testing();
};

subtest 'setbinop callback behavior' => sub {
  for my $impl (@impl) {
    my ($label,$func) = @$impl;

    my (@a,@b);
    @a = (1,2);
    @b = (10,20);
    my $mutated = $func->{setbinop}->(
      sub {
        my $v = $a+$b;
        $a = 0;
        $b = 0;
        @a = ();
        @b = ();
        return $v;
      },
      \@a,\@b
    );
    is_deeply(sequence($mutated), [qw/11 12 21 22/],
              "$label snapshots arrays and resets block values");

    @a = (1,2);
    @b = (10,20);
    my $ok;
    {
      local $@;
      $ok = eval {
        $func->{setbinop}->(sub { die "set callback failure\n" },
                            \@a,\@b);
        1;
      };
      ok(!$ok && $@ =~ /set callback failure/,
         "$label propagates callback exceptions");
    }
    is_deeply(\@a,[1,2],"$label exception preserves first input");
    is_deeply(\@b,[10,20],"$label exception preserves second input");
    is_deeply(sequence($func->{setbinop}->(sub { $a+$b },\@a,\@b)),
              [qw/11 12 21 22/],
              "$label remains usable after callback exception");
  }
  done_testing();
};

{
  package MPUSetTieArray;
  require Tie::Array;
  use vars qw(@ISA);
  @ISA = qw(Tie::StdArray);
}

sub tied_array {
  my @a;
  tie @a, 'MPUSetTieArray';
  @a = @_;
  return \@a;
}

subtest 'sparse and magical arrays' => sub {
  for my $impl (@impl) {
    my ($label,$func) = @$impl;

    my $ta = tied_array(3,1,3);
    my $tb = tied_array(2,3);
    is_deeply(sequence($func->{setunion}->($ta,$tb)), [qw/1 2 3/],
              "$label setunion accepts tied arrays");

    $ta = tied_array(3,1,2);
    $tb = tied_array(2,3,1);
    is($func->{set_is_equal}->($ta,$tb),1,
       "$label relation accepts tied arrays");

    $ta = tied_array(1,3,5);
    $tb = tied_array(1,5);
    is($func->{setcontains}->($ta,$tb),1,
       "$label containment accepts tied arrays in set form");

    $ta = tied_array(1,3,5);
    $tb = tied_array(2,4);
    is($func->{setinsert}->($ta,$tb),2,
       "$label mutator returns count for tied array");
    is_deeply(sequence($ta),[qw/1 2 3 4 5/],
              "$label mutates tied array correctly");
    is_deeply(sequence($tb),[qw/2 4/],
              "$label preserves tied second input");

    my (@sparse1,@sparse2);
    $sparse1[0] = 1;
    $sparse1[2] = 3;
    $sparse2[0] = 1;
    $sparse2[2] = 3;
    my @bad = (
      [setunion => sub { $func->{setunion}->(\@sparse1,[2]) }],
      [sumset => sub { $func->{sumset}->(\@sparse1,[2]) }],
      [setbinop => sub {
         $func->{setbinop}->(sub { $a+$b },\@sparse1,[2])
       }],
      [set_is_equal => sub {
         $func->{set_is_equal}->(\@sparse1,\@sparse2)
       }],
      [is_sidon_set => sub { $func->{is_sidon_set}->(\@sparse1) }],
      [is_sumfree_set => sub { $func->{is_sumfree_set}->(\@sparse1) }],
    );
    for my $test (@bad) {
      my ($name,$code) = @$test;
      my $ok;
      {
        local $SIG{__WARN__} = sub {};
        $ok = eval { $code->(); 1 };
      }
      ok(!$ok,"$label $name rejects sparse arrays");
    }
  }
  done_testing();
};

done_testing();
