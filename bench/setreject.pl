#!perl
use strict;
use warnings;
use feature 'say';

# See https://perlmonks.org/?node_id=11165809

#no warnings "experimental::keyword_any";
#use experimental 'keyword_any';

use ntheory qw/shuffle vecany setcontains setcontainsany set_is_subset/;
use Set::Tiny;
use Benchmark 'cmpthese';

my $data = '';
for my $r ( shuffle 0 .. 31 ) {
    for my $c ( shuffle 0 .. 31 ) {
        $data .= "$c $r whatever\n"
    }
}
my @skip = ( 0, 15, 16, 31 );

# If we do this, we have many values to check.
# The vecany method (or List::Util any) slows down drastically.
# regex and index slow down too, though not as much.
# setissubset and setcontainsany slow down a very small amount
#my @skip = map { 2*$_ } 0..200;

my %skiphash;
@skiphash{@skip} = undef;

my $skipstr = join " ",@skip;
my $skipidxstr = " $skipstr ";

my $skipset = Set::Tiny->new(@skip);

cmpthese -2, {
    # Perl 5.42, performance is approximately equal to vecany
    #any => sub {
    #    while ( $data =~ /^(\d+) (\d+)/mg ) {
    #        next if any { $1 == $_ || $2 == $_ } @skip;
    #    }
    #    return 1
    #},
    # This is identical to List::Util any
    vecany => sub {
        while ( $data =~ /^(\d+) (\d+)/mg ) {
            next if vecany { $1 == $_ || $2 == $_ } @skip;
        }
        return 1
    },
    # This version is handing in magic variables, so parsing is slow.
    setcontains => sub {
        while ( $data =~ /^(\d+) (\d+)/mg ) {
            next if setcontains(\@skip,$1) || setcontains(\@skip,$2);
        }
        return 1
    },
    Lsetcontains => sub {
        my @matches = ($data =~ /^(\d+) (\d+)/mg);
        while (my ($c,$r) = splice @matches,0,2) {
            next if setcontains(\@skip,$c) || setcontains(\@skip,$r);
        }
        return 1
    },
    # Here we force the input into a numerical value so parsing is very fast.
    setcontains0 => sub {
        while ( $data =~ /^(\d+) (\d+)/mg ) {
            next if setcontains(\@skip,0+$1) || setcontains(\@skip,0+$2);
        }
        return 1
    },
    # Putting $1 etc into an array ref de-magics it (see with Devel::Peek)
    # It might be very slightly faster with 0+
    setcontainsany => sub {
        while ( $data =~ /^(\d+) (\d+)/mg ) {
            next if setcontainsany(\@skip,[$1,$2]);
        }
        return 1
    },
    setissubset => sub {
        while ( $data =~ /^(\d+) (\d+)/mg ) {
            next if set_is_subset(\@skip,[$1]) || set_is_subset(\@skip,[$2]);
        }
        return 1
    },
    # Can use regex to search a string list
    regex => sub {
        while ( $data =~ /^(\d+) (\d+)/mg ) {
            my($s,$t)=($1,$2);
            next if $skipstr =~ /\b$s\b/ || $skipstr =~ /\b$t\b/;
        }
        return 1
    },
    # We're doing things the Perl 4 way now.  :)
    index => sub {
        while ( $data =~ /^(\d+) (\d+)/mg ) {
            my($s,$t)=(" $1 "," $2 ");
            next if index($skipidxstr,$s) >= 0 || index($skipidxstr,$t) >= 0;
        }
        return 1
    },
    # Rather obvious hash solution
    hash => sub {
        while ( $data =~ /^(\d+) (\d+)/mg ) {
            next if exists $skiphash{$1} or exists $skiphash{$2};
        }
        return 1;
    },
    # Set::Tiny uses hashes underneath.  This looks pretty.
    settiny => sub {
        while ( $data =~ /^(\d+) (\d+)/mg ) {
            next if $skipset->contains($1) || $skipset->contains($2);
        }
        return 1;
    },
};

__END__

                 Rate regex  any vecany setissubset settiny setcontainsany setcontains index setcontains0 hash
regex           721/s    -- -55%   -58%        -58%    -61%           -68%        -68%  -70%         -72% -83%
any            1594/s  121%   --    -7%         -8%    -15%           -29%        -29%  -34%         -37% -63%
vecany         1715/s  138%   8%     --         -1%     -8%           -23%        -23%  -29%         -32% -60%
setissubset    1730/s  140%   9%     1%          --     -7%           -22%        -23%  -28%         -32% -59%
settiny        1866/s  159%  17%     9%          8%      --           -16%        -16%  -23%         -26% -56%
setcontainsany 2231/s  209%  40%    30%         29%     20%             --         -0%   -8%         -12% -48%
setcontains    2235/s  210%  40%    30%         29%     20%             0%          --   -8%         -12% -48%
index          2418/s  235%  52%    41%         40%     30%             8%          8%    --          -5% -43%
setcontains0   2536/s  252%  59%    48%         47%     36%            14%         13%    5%           -- -40%
hash           4257/s  490% 167%   148%        146%    128%            91%         90%   76%          68%   --
