#!/usr/bin/perl
use warnings;
use strict;
use lib 'lib';

# Copied and adapted from Set::Tiny 0.06.
# Non-representative benchmark of different Set:: modules

use Benchmark qw( cmpthese );

use Set::Tiny;
#use Set::Scalar;
#use Set::Object;
use Math::Prime::Util qw/toset setinsert setcontains setremove setinvert
                         setintersect setunion setminus setdelta vecequal/;

#my @a = 1 .. 100;
#my @b = 51 .. 150;
my @a = 1 .. 1000;
my @b = 501 .. 2500;

my $s_t1 = Set::Tiny->new(@a);
my $s_t2 = Set::Tiny->new(@b);
my $s_m1 = [toset(\@a)];
my $s_m2 = [toset(\@b)];
#my $s_s1 = Set::Scalar->new(@a);
#my $s_s2 = Set::Scalar->new(@b);
#my $s_o1 = Set::Object->new(@a);
#my $s_o2 = Set::Object->new(@b);

my %tests = (
    A_new => {
        t => sub { Set::Tiny->new(@a) },
        m => sub { [toset(\@a)] },
        #s => sub { Set::Scalar->new(@a) },
        #o => sub { Set::Object->new(@a) },
    },

    A_clone => {
        t => sub { $s_t1->clone },
        m => sub { [@$s_m1] },
        #s => sub { $s_s1->clone },
        #o => sub { },
    },
    A_insert => {
        t => sub { Set::Tiny->new->insert(@a) },
        m => sub { setinsert([],\@a); },
        #s => sub { Set::Scalar->new->insert(@a) },
        #o => sub { Set::Object->new->insert(@a) },
    },
    A_delete => {
        t => sub { Set::Tiny->new(@a)->delete(@b) },
        m => sub { setremove([toset(\@a)],\@b); },
        #s => sub { Set::Scalar->new(@a)->delete(@b) },
        #o => sub { Set::Object->new(@a)->delete(@b) },
    },
    A_invert => {
        t => sub { Set::Tiny->new(@a)->invert(@b) },
        m => sub { setinvert([toset(\@a)],\@b); },
        #s => sub { Set::Scalar->new(@a)->invert(@b) },
        #o => sub { Set::Object->new(@a)->invert(@b) },
    },
    C_is_equal => {
        t => sub { $s_t1->is_equal($s_t2) },
        #m => sub { vecequal($s_m1,$s_m2) },    # Probably much faster
        m => sub { Math::Prime::Util::set_is_equal($s_m1,$s_m2) },
        #s => sub { $s_s1->is_equal($s_s2) },
        #o => sub { $s_o1->equal($s_o2) },
    },

#  Set::Tiny  $s->is_subset($t)     is $s a subset of $t?
#  MPU:       set_is_subset($s,$t)  is $t a subset of $s?
    C_is_subset => {
        t => sub { $s_t1->is_subset($s_t2) },
        #m => sub { setcontains($s_m1, $s_m2) },    # Probably faster
        m => sub { Math::Prime::Util::set_is_subset($s_m2,$s_m1); },
        #s => sub { $s_s1->is_subset($s_s2) },
        #o => sub { $s_o1->subset($s_o2) },
    },
    C_is_proper_subset => {
        t => sub { $s_t1->is_proper_subset($s_t2) },
        m => sub { Math::Prime::Util::set_is_proper_subset($s_m2,$s_m1); },
        #s => sub { $s_s1->is_proper_subset($s_s2) },
        #o => sub { $s_o1->proper_subset($s_o2) },
    },
    C_is_superset => {
        t => sub { $s_t1->is_superset($s_t2) },
        m => sub { Math::Prime::Util::set_is_superset($s_m2,$s_m1); },
        #s => sub { $s_s1->is_superset($s_s2) },
        #o => sub { $s_o1->superset($s_o2) },
    },
    C_is_proper_superset => {
        t => sub { $s_t1->is_proper_superset($s_t2) },
        m => sub { Math::Prime::Util::set_is_proper_superset($s_m2,$s_m1); },
        #s => sub { $s_s1->is_proper_superset($s_s2) },
        #o => sub { $s_o1->proper_superset($s_o2) },
    },
    C_is_disjoint => {
        t => sub { $s_t1->is_disjoint($s_t2) },
        m => sub { Math::Prime::Util::set_is_disjoint($s_m1,$s_m2); },
        #s => sub { $s_s1->is_disjoint($s_s2) },
        #o => sub { $s_o1->is_disjoint($s_o2) },
    },

    # The $set->contains(@elements) methods are not identical:
    # MPU, Set::Tiny, Set::Object return true if $set contains *all* elements.
    # Set::Scalar returns true if $set contains *any* elements.

    B_contains => {
        t => sub { $s_t1->contains(@b) },
        m => sub { setcontains($s_m1,\@b) },
        #s => sub { $s_s1->contains(@b) },
        #o => sub { $s_o1->contains(@b) },
    },
#  Set::Tiny  $s->difference($t)     $s minus $t
#  MPU:       set_is_subset($s,$t)   $s minus $t
    B_difference => {
        t => sub { $s_t1->difference($s_t2) },
        m => sub { setminus($s_m1,$s_m2) },
        #s => sub { $s_s1->difference($s_s2) },
        #o => sub { $s_o1->difference($s_o2) },
    },
    B_union => {
        t => sub { $s_t1->union($s_t2) },
        m => sub { setunion($s_m1,$s_m2) },
        #s => sub { $s_s1->union($s_s2) },
        #o => sub { $s_o1->union($s_o2) },
    },
    B_intersection => {
        t => sub { $s_t1->intersection($s_t2) },
        m => sub { setintersect($s_m1,$s_m2) },
        #s => sub { $s_s1->intersection($s_s2) },
        #o => sub { $s_o1->intersection($s_o2) },
    },
    B_symmetric_difference => {
        t => sub { $s_t1->symmetric_difference($s_t2) },
        m => sub { setdelta($s_m1,$s_m2) },
        #s => sub { $s_s1->symmetric_difference($s_s2) },
        #o => sub { $s_o1->symmetric_difference($s_o2) },
    },
);

print "running benchmarks with sets of size ",
    scalar @a, " and ", scalar @b, "\n";
for my $test ( sort keys %tests ) {
    print "\n$test:\n";
    cmpthese(
        -1,
        {
            'Set::Tiny'   => $tests{$test}{t},
            #'Set::Scalar' => $tests{$test}{s},
            #'Set::Object' => $tests{$test}{o},
            'MPU' => $tests{$test}{m},
        }
    );
}

