package ntheory;
use strict;
use warnings;

BEGIN {
  $ntheory::AUTHORITY = 'cpan:DANAJ';
  $ntheory::VERSION = '0.43';
}

BEGIN {
  require Math::Prime::Util; 
  *ntheory:: = *Math::Prime::Util::;
  $INC{"ntheory"}++;
}

1;

__END__


# ABSTRACT: Number theory utilities

=pod

=head1 NAME

ntheory - Number theory utilities

=head1 VERSION

Version 0.43

=head1 SEE

See L<Math::Prime::Util> for complete documentation.

=head1 COPYRIGHT

Copyright 2011-2014 by Dana Jacobsen E<lt>dana@acm.orgE<gt>

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut
