#!/usr/bin/env perl
use strict;
use warnings;

use B ();
use Math::Prime::Util qw/:all/;

sub classify_call {
  my ($op, $name, $seen) = @_;
  return 'OTHER' unless defined($op) && $$op;
  return 'OTHER' if $seen->{$$op}++;

  my $opname = eval { $op->name } // '';
  return 'OPCODE' if ref($op) eq 'B::OP' && $opname eq $name;
  return 'SUB'    if $opname eq 'entersub';

  for my $childmeth (qw/first other sibling/) {
    my $child = eval { $op->$childmeth() };
    my $kind = classify_call($child, $name, $seen);
    return $kind if $kind ne 'OTHER';
  }
  return 'OTHER';
}

sub compiled_kind {
  my ($name, $src) = @_;
  my $cv = eval $src;
  die "eval failed for $name source '$src': $@" if $@ || ref($cv) ne 'CODE';
  my $root = B::svref_2object($cv)->ROOT;
  return classify_call($root, $name, {});
}

my @cases = (
  [ 'addint',        q{ my $y = 56; sub { addint($y, $_) } },      'OPCODE', '$_ ok' ],
  [ 'subint',        q{ my $y = 56; sub { subint($y, $_) } },      'OPCODE', '$_ ok' ],
  [ 'mulint',        q{ my $y = 56; sub { mulint($y, $_) } },      'OPCODE', '$_ ok' ],
  [ 'divint',        q{ my $y = 56; sub { divint($y, $_) } },      'OPCODE', '$_ ok' ],
  [ 'modint',        q{ my $y = 56; sub { modint($y, $_) } },      'OPCODE', '$_ ok' ],
  [ 'cdivint',       q{ my $y = 56; sub { cdivint($y, $_) } },     'OPCODE', '$_ ok' ],
  [ 'powint',        q{ my $y = 56; sub { powint($y, $_) } },      'OPCODE', '$_ ok' ],
  [ 'is_odd',        q{sub { is_odd($_) for 1..10 }},              'OPCODE', '$_ ok' ],
  [ 'is_even',       q{sub { is_even($_) for 1..10 }},             'OPCODE', '$_ ok' ],
  [ 'signint',       q{sub { signint(56) }},                       'OPCODE', 'simple only' ],
  [ 'signint',       q{ my $x = 56; sub { signint($x) } },         'SUB',    'simple only; lexical' ],
  [ 'signint',       q{my $x=0; sub{$x+=signint($_) for -5..5}},   'SUB',    '$_ ok' ],
  [ 'cmpint',        q{sub { cmpint(56, 0) }},                     'OPCODE', 'simple only' ],
  [ 'cmpint',        q{ my ($a,$b) = (56,0); sub { cmpint($a,$b) } }, 'SUB', 'simple only; lexical' ],

  [ 'is_power',      q{sub { is_power(56) }},                      'SUB',    'candidate' ],
  [ 'is_prime',      q{sub { is_prime(56) }},                      'SUB',    'candidate' ],
  [ 'is_square',     q{sub { is_square(56) }},                     'SUB',    'candidate' ],
  [ 'is_square_free',q{sub { is_square_free(56) }},                'SUB',    'candidate' ],
  [ 'is_congruent',  q{sub { is_congruent(10,1,3) }},              'SUB',    'candidate' ],
  [ 'kronecker',     q{sub { kronecker(5,7) }},                    'SUB',    'candidate' ],
  [ 'urandomm',      q{sub { urandomm(10) }},                      'SUB',    'candidate' ],
  [ 'urandomr',      q{sub { urandomr(1,10) }},                    'SUB',    'candidate' ],
  [ 'addmod',        q{sub { addmod(5,7,11) }},                    'SUB',    'candidate' ],
  [ 'mulmod',        q{sub { mulmod(5,7,11) }},                    'SUB',    'candidate' ],
  [ 'setcontains',   q{sub { setcontains([1,2,3],2) }},            'SUB',    'candidate' ],
  [ 'setcontainsany',q{sub { setcontainsany([1,2,3],2) }},         'SUB',    'candidate' ],
  [ 'absint',        q{sub { absint(56) }},                        'OPCODE', 'candidate' ],
  [ 'absint',        q{ my $x = 56; sub { absint($x) } },          'SUB',    'candidate' ],
  [ 'negint',        q{sub { negint(56) }},                        'OPCODE', 'candidate' ],
  [ 'negint',        q{ my $x = 56; sub { negint($x) } },          'SUB',    'candidate' ],
  [ 'sqrtint',       q{sub { sqrtint(56) }},                       'OPCODE', 'candidate' ],
  [ 'sqrtint',       q{ my $x = 56; sub { sqrtint($x) } },         'SUB',    'candidate' ],
);

printf "%-14s %-7s %-7s %s\n", "Function", "Got", "Want", "Notes";
printf "%-14s %-7s %-7s %s\n", "-" x 14, "-" x 7, "-" x 7, "-" x 24;

my $lname = 'xxxxxx';
for my $case (@cases) {
  my ($name, $src, $want, $note) = @$case;
  my $got = compiled_kind($name, $src);
  ($lname,$name) = ($name, $name eq $lname ? '' : $name);
  printf "%-14s %-7s %-7s %s\n", $name, $got, $want, $note;
}

