use ExtUtils::MakeMaker;

my $broken64 = (18446744073709550592 == ~0);
if ($broken64) {
  warn <<EOW;

Your Perl has a broken 64-bit implementation.
Arithmetic operations on numbers larger than 2^53 are wrong.

Much core functionality will work, but some functions, especially
random primes, will not work for large inputs.

I highly recommend upgrading to a newer version of Perl.

EOW
}

WriteMakefile1(
    NAME         => 'Math::Prime::Util',
    ABSTRACT     => 'Utilities related to prime numbers, including fast sieves and factoring',
    VERSION_FROM => 'lib/Math/Prime/Util.pm',
    LICENSE      => 'perl',
    AUTHOR       => 'Dana A Jacobsen <dana@acm.org>',

    OBJECT       => 'cache.o '    .
                    'factor.o '   .
                    'primality.o '.
                    'aks.o '      .
                    'lehmer.o '   .
                    'lmo.o '      .
                    'sieve.o '    .
                    'util.o '     .
                    'XS.o',
    LIBS         => ['-lm'],

    EXE_FILES    => ['bin/primes.pl', 'bin/factor.pl'],

    TEST_REQUIRES=> {
                      'Test::More'       => '0.45',
                      'bignum'           => '0.22',   # 'use bigint' in tests
                    },
    PREREQ_PM    => {
                      'Exporter'         => '5.562',
                      'XSLoader'         => '0.01',
                      'Carp'             => ($] < 5.008) ? '1.17' : 0,
                      'Tie::Array'       => 0,
                      'base'             => 0,
                      'constant'         => 0,
                      'Config'           => 0,
                      # 1.99 fixes the FastCalc SvUV bug, we work around it.
                      'Math::BigInt'     => '1.88',
                      'Math::BigFloat'   => '1.59',
                      'Bytes::Random::Secure' => '0.23',
                    },
    META_MERGE   => {
               'meta-spec'  => {
                  version => '2',
                  url     => 'http://search.cpan.org/perldoc?CPAN::Meta::Spec',
               },
               resources    => {
                  license    => [ 'http://dev.perl.org/licenses/' ],
                  homepage   => 'https://github.com/danaj/Math-Prime-Util',
                  repository => {
                    url  => 'https://github.com/danaj/Math-Prime-Util',
                  },
               },
               prereqs => {
                 runtime => {
                   recommends => {
                     'Math::Prime::Util::GMP' => 0.20,
                     'Math::BigInt::GMP'      => 0,
                     'Math::MPFR'             => 2.03,
                     'Digest::SHA'            => 5.87,
                   },
                 },
                 test  => {
                   suggests => {
                     'Test::Warn' => 0,
                   },
                 },
               },
             },

    MIN_PERL_VERSION => 5.006002,
);

sub WriteMakefile1 {   # Cribbed from eumm-upgrade by Alexandr Ciornii
  my %params = @_;
  my $eumm_version = $ExtUtils::MakeMaker::VERSION;
  $eumm_version = eval $eumm_version;

  if ($params{TEST_REQUIRES} and $eumm_version < 6.6303) {
      $params{BUILD_REQUIRES}={ %{$params{BUILD_REQUIRES} || {}} , %{$params{TEST_REQUIRES}} };
      delete $params{TEST_REQUIRES};
  }
  if ($params{BUILD_REQUIRES} and $eumm_version < 6.5503) {
      #EUMM 6.5502 has problems with BUILD_REQUIRES
      $params{PREREQ_PM}={ %{$params{PREREQ_PM} || {}} , %{$params{BUILD_REQUIRES}} };
      delete $params{BUILD_REQUIRES};
  }
  delete $params{CONFIGURE_REQUIRES} if $eumm_version < 6.52;
  delete $params{MIN_PERL_VERSION} if $eumm_version < 6.48;
  delete $params{META_MERGE} if $eumm_version < 6.46;
  delete $params{META_ADD} if $eumm_version < 6.46;
  delete $params{LICENSE} if $eumm_version < 6.31;
  delete $params{AUTHOR} if $] < 5.005;
  delete $params{ABSTRACT_FROM} if $] < 5.005;
  delete $params{BINARY_LOCATION} if $] < 5.005;

  WriteMakefile(%params);
}
