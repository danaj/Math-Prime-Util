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
# say $data; die;
my @skip = ( 0, 15, 16, 31 );

my %skiphash;
@skiphash{@skip} = undef;

my $skipstr = join " ",@skip;
my $skipidxstr = " $skipstr ";

my $skipset = Set::Tiny->new(@skip);

cmpthese -2, {
    # Perl 5.42
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

                 Rate regex vecany setissubset settiny setcontainsany index setcontains setcontains0 hash
regex           680/s    --   -58%        -60%    -61%           -65%  -70%        -70%         -73% -82%
vecany         1618/s  138%     --         -4%     -7%           -17%  -28%        -29%         -37% -58%
setissubset    1680/s  147%     4%          --     -3%           -13%  -25%        -27%         -34% -56%
settiny        1732/s  155%     7%          3%      --           -11%  -23%        -24%         -32% -55%
setcontainsany 1939/s  185%    20%         15%     12%             --  -13%        -15%         -24% -49%
index          2239/s  229%    38%         33%     29%            15%    --         -2%         -12% -42%
setcontains    2291/s  237%    42%         36%     32%            18%    2%          --         -10% -40%
setcontains0   2550/s  275%    58%         52%     47%            32%   14%         11%           -- -33%
hash           3829/s  463%   137%        128%    121%            97%   71%         67%          50%   --
