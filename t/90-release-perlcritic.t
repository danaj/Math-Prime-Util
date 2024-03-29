#!/usr/bin/env perl
use strict;
use warnings;
use Test::More;

BEGIN {
  unless ($ENV{RELEASE_TESTING}) {
    plan( skip_all => 'these tests are for release candidate testing' );
  }
}

#---------------------------------------------------------------------


eval { require Test::Perl::Critic; };
plan skip_all => "Test::Perl::Critic required for testing PBP compliance" if $@;

Test::Perl::Critic->import(
        -verbose => 10,
        -severity => 'gentle',   # default
        -force => 0,             # default (allow ## no critic)
        # We probably shouldn't do this, but this is tiresome
        -exclude => [qw/ProhibitExplicitReturnUndef/],
       );

all_critic_ok();
