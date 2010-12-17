package AlignDB::Window;
# ABSTRACT: Split integer spans into a series of windows

use Moose;
use Carp;

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::IntSpan;

=attr sw_size

sliding windows' size, default is 100

=cut

has 'sw_size' => ( is => 'rw', isa => 'Int', default => 100, );

=attr min_interval

mininal indel interval length, default is 11

=cut

has 'min_interval' => ( is => 'rw', isa => 'Int', default => 11, );

=attr max_out_distance

maximal outside distance, default is 10

=cut

has 'max_out_distance' => ( is => 'rw', isa => 'Int', default => 10, );

=attr max_in_distance

maximal inside distance, default is 5

=cut

has 'max_in_distance' => ( is => 'rw', isa => 'Int', default => 5, );

=method interval_window

      Usage : $self->interval_window(
            :     $comparable_set, $interval_start, $interval_end,
            :     $window_size, $minimal_interval,
            : );
    Purpose : split an interval to windows
            : length of windows are variable, but all positions of
            :   the interval are counted
    Returns : @interval_windows
            : each member is a hash_ref
 Parameters : $comparable_set   : AlignDB::IntSpan object
            : $interval_start   : start position of the interval
            : $interval_end     : end position of the interval
            : $window_size      : size of windows
            : $minimal_interval : minimal size of intervals
     Throws : no exceptions
   Comments : none 
   See Also : interval_window_2

=cut

sub interval_window {
    my ( $self, $comparable_set, $interval_start, $interval_end, $window_size,
        $minimal_interval )
        = @_;

    # if undefined, set to default values
    $window_size      ||= $self->sw_size;
    $minimal_interval ||= $self->min_interval;

    my $comparable_number = $comparable_set->cardinality;

    my $interval_set = AlignDB::IntSpan->new("$interval_start-$interval_end");
    $interval_set = $interval_set->intersect($comparable_set);
    my $interval_length = $interval_set->cardinality;

    my $density = int( $interval_length / $window_size ) - 1;

    # used for mark windows
    my ( $tmp_start, $tmp_end ) = ( 1, $interval_length );

    my @interval_windows;

    # More windows will be submitted in the following for-loop
    #   and if section
    if ( $interval_length < $minimal_interval ) {
        return @interval_windows;
    }
    elsif ( $interval_length < $window_size ) {

        # $interval_length between $minimal_interval and 99 bp
        my %window_info;
        $window_info{set}      = $interval_set;
        $window_info{distance} = -1;
        $window_info{density}  = $density;
        $window_info{type}     = 'S';
        push @interval_windows, \%window_info;

        return @interval_windows;
    }
    elsif ( $interval_length < $window_size * 2 ) {

        # $interval_length between 100 and 199 bp
        # LO between 50 and 100 bp
        # R0 between 50 and 99 bp
        # left windows
        my %l_window_info;
        my $l_start = $tmp_start;
        my $l_end   = $l_start + int( $interval_length / 2 ) - 1;
        if ( $interval_length % 2 ) {
            $l_end++;
        }
        $l_window_info{set}      = $interval_set->slice( $l_start, $l_end );
        $l_window_info{distance} = 0;
        $l_window_info{density}  = $density;
        $l_window_info{type}     = 'L';
        push @interval_windows, \%l_window_info;

        # right windows
        my %r_window_info;
        my $r_end   = $tmp_end;
        my $r_start = $r_end - int( $interval_length / 2 ) + 1;
        $r_window_info{set}      = $interval_set->slice( $r_start, $r_end );
        $r_window_info{distance} = 0;
        $r_window_info{density}  = $density;
        $r_window_info{type}     = 'R';
        push @interval_windows, \%r_window_info;

        return @interval_windows;
    }
    else {

        # LO => 50 bp, R0 => 50 bp
        # left windows
        my %l_window_info;
        my $l_start = $tmp_start;
        my $l_end   = $l_start + $window_size / 2 - 1;
        $l_window_info{set}      = $interval_set->slice( $l_start, $l_end );
        $l_window_info{distance} = 0;
        $l_window_info{density}  = $density;
        $l_window_info{type}     = 'L';
        push @interval_windows, \%l_window_info;
        $tmp_start = $l_end + 1;

        # right windows
        my %r_window_info;
        my $r_end   = $tmp_end;
        my $r_start = $r_end - $window_size / 2 + 1;
        $r_window_info{set}      = $interval_set->slice( $r_start, $r_end );
        $r_window_info{distance} = 0;
        $r_window_info{density}  = $density;
        $r_window_info{type}     = 'R';
        push @interval_windows, \%r_window_info;
        $tmp_end = $r_start - 1;
    }

    # regular 100 bp windows
    # When density IN (1, 2), this loop will be bypassed.
    my $half_windows_number = int( ( $density - 1 ) / 2 );
    for ( my $i = 1; $i <= $half_windows_number; $i++ ) {

        # left windows
        my %l_window_info;
        my $l_start = $tmp_start;
        my $l_end   = $l_start + $window_size - 1;
        $l_window_info{set}      = $interval_set->slice( $l_start, $l_end );
        $l_window_info{distance} = $i;
        $l_window_info{density}  = $density;
        $l_window_info{type}     = 'L';
        push @interval_windows, \%l_window_info;
        $tmp_start = $l_end + 1;

        # right windows
        my %r_window_info;
        my $r_end   = $tmp_end;
        my $r_start = $r_end - $window_size + 1;
        $r_window_info{set}      = $interval_set->slice( $r_start, $r_end );
        $r_window_info{distance} = $i;
        $r_window_info{density}  = $density;
        $r_window_info{type}     = 'R';
        push @interval_windows, \%r_window_info;
        $tmp_end = $r_start - 1;
    }

    # Three special windwos:
    #   1. Ln, range from 100 to 150 bp when even density number
    #   2. Rn, range from 100 to 149 bp when even density number
    #   3. Ln or Rn, range from 100 to 199 bp when odd density number
    if ( $density % 2 ) {

        # middle windows
        my %window_info;
        $window_info{set}      = $interval_set->slice( $tmp_start, $tmp_end );
        $window_info{distance} = $half_windows_number + 1;
        $window_info{density}  = $density;
        $window_info{type}     = rand > 0.5 ? "L" : "R";
        push @interval_windows, \%window_info;
    }
    else {
        my $remain_length = $tmp_end - $tmp_start + 1;

        # left windows
        my %l_window_info;
        my $l_start = $tmp_start;
        my $l_end   = $l_start + int( $remain_length / 2 ) - 1;
        if ( $remain_length % 2 ) {
            $l_end++;
        }
        $l_window_info{set}      = $interval_set->slice( $l_start, $l_end );
        $l_window_info{distance} = $half_windows_number + 1;
        $l_window_info{density}  = $density;
        $l_window_info{type}     = 'L';
        push @interval_windows, \%l_window_info;

        # right windows
        my %r_window_info;
        my $r_end   = $tmp_end;
        my $r_start = $r_end - int( $remain_length / 2 ) + 1;
        $r_window_info{set}      = $interval_set->slice( $r_start, $r_end );
        $r_window_info{distance} = $half_windows_number + 1;
        $r_window_info{density}  = $density;
        $r_window_info{type}     = 'R';
        push @interval_windows, \%r_window_info;
    }

    return @interval_windows;
}

=method interval_window_2

      Usage : $self->interval_window_2(
            :     $comparable_set, $interval_start, $interval_end,
            :     $window_size, $minimal_interval,
            : );
    Purpose : split an interval to windows
            : all windows are 100 bp length
    Returns : @interval_windows
            : each member is a hash_ref
 Parameters : $comparable_set   : AlignDB::IntSpan object
            : $interval_start   : start position of the interval
            : $interval_end     : end position of the interval
            : $window_size      : size of windows
            : $minimal_interval : minimal size of intervals
     Throws : no exceptions
   Comments : none 
   See Also : interval_window

=cut

sub interval_window_2 {
    my ( $self, $comparable_set, $interval_start, $interval_end, $window_size,
        $minimal_interval )
        = @_;

    # if undefined, set to default values
    $window_size      ||= $self->sw_size;
    $minimal_interval ||= $self->min_interval;

    my $comparable_number = $comparable_set->cardinality;

    ## Don't need these
    #my $start = &ceil_nearest($comparable_set,$interval_start);
    #my $end = &floor_nearest($comparable_set,$interval_end);

    my $interval_set = AlignDB::IntSpan->new("$interval_start-$interval_end");
    $interval_set = $interval_set->intersect($comparable_set);
    my $interval_length = $interval_set->cardinality;

    my $density = int( $interval_length / $window_size );

    # used for mark windows
    my ( $tmp_start, $tmp_end ) = ( 1, $interval_length );

    my @interval_windows;

    # More windows will be submitted in the following for-loop
    #   and if section
    if ( $interval_length < $minimal_interval ) {
        return @interval_windows;
    }
    elsif ( $interval_length < 2 * $window_size ) {

        # $interval_length between $minimal_interval and 199 bp
        my %window_info;
        $window_info{set}      = $interval_set;
        $window_info{distance} = 0;
        $window_info{density}  = $density;
        $window_info{type}     = 'S';
        push @interval_windows, \%window_info;

        return @interval_windows;
    }

    # regular 100 bp windows
    my $half_windows_number = int( ( $density - 1 ) / 2 );
    for my $i ( 1 .. $half_windows_number ) {

        # left windows
        my %l_window_info;
        my $l_start = $tmp_start;
        my $l_end   = $l_start + $window_size - 1;
        $l_window_info{set}      = $interval_set->slice( $l_start, $l_end );
        $l_window_info{distance} = $i;
        $l_window_info{density}  = $density;
        $l_window_info{type}     = 'L';
        push @interval_windows, \%l_window_info;
        $tmp_start = $l_end + 1;

        # right windows
        my %r_window_info;
        my $r_end   = $tmp_end;
        my $r_start = $r_end - $window_size + 1;
        $r_window_info{set}      = $interval_set->slice( $r_start, $r_end );
        $r_window_info{distance} = $i;
        $r_window_info{density}  = $density;
        $r_window_info{type}     = 'R';
        push @interval_windows, \%r_window_info;
        $tmp_end = $r_start - 1;
    }

    # Special middle windwo:
    #   Ln or Rn when odd density number
    if ( $density % 2 ) {
        my %window_info;
        $window_info{type} = rand > 0.5 ? "L" : "R";
        my ( $window_start, $window_end );
        if ( $window_info{type} eq "L" ) {
            $window_start = $tmp_start;
            $window_end   = $window_start + $window_size - 1;
        }
        else {
            $window_end   = $tmp_end;
            $window_start = $window_end - $window_size + 1;
        }

        $window_info{set}
            = $interval_set->slice( $window_start, $window_end );
        $window_info{distance} = $half_windows_number + 1;
        $window_info{density}  = $density;

        push @interval_windows, \%window_info;
    }

    return @interval_windows;
}

=method outside_window

      Usage : $self->outside_window(
            :     $comparable_set, $interval_start, $interval_end,
            :     $window_size, $maximal_distance,
            : );
    Purpose : draw outside windows from a internal region
            : all windows are 100 bp length
            : start from 1 and end to $maximal_distance
    Returns : @outside_windows
            : each member is a hash_ref
 Parameters : $comparable_set   : AlignDB::IntSpan object
            : $interval_start   : start position of the interval
            : $interval_end     : end position of the interval
            : $window_size      : size of windows
            : $maximal_distance : maximal distance
     Throws : no exceptions
   Comments : none 
   See Also : interval_window, interval_window_2, inside_window

=cut

sub outside_window {
    my ( $self, $comparable_set, $internal_start, $internal_end, $window_size,
        $maximal_distance )
        = @_;

    # if undefined, set to default values
    $window_size      ||= $self->sw_size;
    $maximal_distance ||= $self->max_out_distance;

    my $comparable_number = $comparable_set->cardinality;

    my @outside_windows;

    foreach my $sw_type (qw/L R/) {

        # $sw_start and $sw_end are both index of $target_set
        my ( $sw_start, $sw_end );

        # If gene_set containing promoter or terminator sequences,
        #   use transcript start or end as gene start or end
        if ( $sw_type eq 'R' ) {
            $sw_start = $comparable_set->index($internal_end) + 1;
            $sw_end   = $sw_start + $window_size - 1;
        }
        elsif ( $sw_type eq 'L' ) {
            $sw_end   = $comparable_set->index($internal_start) - 1;
            $sw_start = $sw_end - $window_size + 1;
        }

        # $gsw_distance is from 1 to $sw_max_distance
    OUTSIDESW: foreach my $sw_distance ( 1 .. $maximal_distance ) {
            last if $sw_start < 1;
            last if $sw_end > $comparable_number;
            my $sw_set = $comparable_set->slice( $sw_start, $sw_end );
            my $sw_set_member_number = $sw_set->cardinality;
            if ( $sw_set_member_number < $window_size ) {
                last OUTSIDESW;
            }

            my %window_info;
            $window_info{type}     = $sw_type;
            $window_info{set}      = $sw_set;
            $window_info{distance} = $sw_distance;

            push @outside_windows, \%window_info;

            if ( $sw_type eq 'R' ) {
                $sw_start = $sw_end + 1;
                $sw_end   = $sw_start + $window_size - 1;
            }
            elsif ( $sw_type eq 'L' ) {
                $sw_end   = $sw_start - 1;
                $sw_start = $sw_end - $window_size + 1;
            }
        }
    }

    return @outside_windows;
}

=method outside_window_2

      Usage : $self->outside_window_2(
            :     $comparable_set, $interval_start, $interval_end,
            :     $window_size, $maximal_distance,
            : );
    Purpose : draw outside windows from a internal region
            : the first window is 50 bp and all others are 100 bp length
            : start from 0 and end to $maximal_distance
    Returns : @outside_windows
            : each member is a hash_ref
 Parameters : $comparable_set   : AlignDB::IntSpan object
            : $interval_start   : start position of the interval
            : $interval_end     : end position of the interval
            : $window_size      : size of windows
            : $maximal_distance : maximal distance
     Throws : no exceptions
   Comments : none 
   See Also : interval_window, interval_window_2, inside_window

=cut

sub outside_window_2 {
    my ( $self, $comparable_set, $internal_start, $internal_end, $window_size,
        $maximal_distance )
        = @_;

    # if undefined, set to default values
    $window_size      ||= $self->sw_size;
    $maximal_distance ||= $self->max_out_distance;

    my $window0_size = int( $window_size / 2 );

    my $comparable_number = $comparable_set->cardinality;

    my @outside_windows;

    foreach my $sw_type (qw/L R/) {

        # $sw_start and $sw_end are both index of $target_set
        my ( $sw_start, $sw_end );

        # If gene_set containing promoter or terminator sequences,
        #   use transcript start or end as gene start or end
        if ( $sw_type eq 'R' ) {
            $sw_start = $comparable_set->index($internal_end) + 1;
            $sw_end   = $sw_start + $window0_size - 1;
        }
        elsif ( $sw_type eq 'L' ) {
            $sw_end   = $comparable_set->index($internal_start) - 1;
            $sw_start = $sw_end - $window0_size + 1;
        }

        # distance is from 0 to $maximal_distance
    OUTSIDESW2: foreach my $sw_distance ( 0 .. $maximal_distance ) {
            last if $sw_start < 1;
            last if $sw_end > $comparable_number;
            my $sw_set = $comparable_set->slice( $sw_start, $sw_end );
            my $sw_set_member_number = $sw_set->cardinality;
            if ( $sw_set_member_number < $window_size ) {
                last OUTSIDESW2;
            }

            my %window_info;
            $window_info{type}     = $sw_type;
            $window_info{set}      = $sw_set;
            $window_info{distance} = $sw_distance;

            push @outside_windows, \%window_info;

            if ( $sw_type eq 'R' ) {
                $sw_start = $sw_end + 1;
                $sw_end   = $sw_start + $window_size - 1;
            }
            elsif ( $sw_type eq 'L' ) {
                $sw_end   = $sw_start - 1;
                $sw_start = $sw_end - $window_size + 1;
            }
        }
    }

    return @outside_windows;
}

=method inside_window

      Usage : $self->inside_window(
            :     $comparable_set, $interval_start, $interval_end,
            :     $window_size, $maximal_distance,
            : );
    Purpose : draw inside windows from a internal region
            : all windows are 100 bp length
            : start counting from the edges
    Returns : @inside_windows
            : each member is a hash_ref
 Parameters : $comparable_set   : AlignDB::IntSpan object
            : $interval_start   : start position of the interval
            : $interval_end     : end position of the interval
            : $window_size      : size of windows
            : $maximal_distance : maximal distance
     Throws : no exceptions
   Comments : none 
   See Also : interval_window, interval_window_2, outside_window

=cut

sub inside_window {
    my ( $self, $comparable_set, $internal_start, $internal_end, $window_size,
        $maximal_distance )
        = @_;

    # if undefined, set to default values
    $window_size      ||= $self->sw_size;
    $maximal_distance ||= $self->max_in_distance;

    my $comparable_number = $comparable_set->cardinality;

    my @inside_windows;

    foreach my $sw_type (qw/l r/) {

        # $sw_start and $sw_end are both index of $target_set
        my ( $sw_start, $sw_end );

        my $working_set
            = $comparable_set->intersect("$internal_start-$internal_end");
        my $working_length = $working_set->cardinality;
        last if $working_length < $window_size;

        # the 'r' and 'l' windows may overlap each other
        if ( $sw_type eq 'r' ) {
            $sw_end   = $working_set->cardinality;
            $sw_start = $sw_end - $window_size + 1;
        }
        elsif ( $sw_type eq 'l' ) {
            $sw_start = 1;
            $sw_end   = $sw_start + $window_size - 1;
        }

        my $available_distance
            = int( $working_length / ( $window_size * 2 ) );
        my $max_distance = min( $available_distance, $maximal_distance );

        # sw_distance is from -1 to -max_distance
    INSIDESW: foreach my $i ( 1 .. $max_distance ) {
            my $sw_set = $working_set->slice( $sw_start, $sw_end );
            my $sw_set_member_number = $sw_set->cardinality;
            if ( $sw_set_member_number < $window_size ) {
                last INSIDESW;
            }
            my $sw_distance = -$i;

            my %window_info;
            $window_info{type}     = $sw_type;
            $window_info{set}      = $sw_set;
            $window_info{distance} = $sw_distance;

            push @inside_windows, \%window_info;

            if ( $sw_type eq 'r' ) {
                $sw_end   = $sw_start - 1;
                $sw_start = $sw_end - $window_size + 1;
            }
            elsif ( $sw_type eq 'l' ) {
                $sw_start = $sw_end + 1;
                $sw_end   = $sw_start + $window_size - 1;
            }
        }
    }

    return @inside_windows;
}

=method inside_window2

      Usage : $self->inside_window2(
            :     $comparable_set, $interval_start, $interval_end,
            :     $window_size, $maximal_distance,
            : );
    Purpose : draw inside windows from a internal region
            : all windows are 100 bp length
            : start counting from the center
    Returns : @inside_windows
            : each member is a hash_ref
 Parameters : $comparable_set   : AlignDB::IntSpan object
            : $interval_start   : start position of the interval
            : $interval_end     : end position of the interval
            : $window_size      : size of windows
            : $maximal_distance : maximal distance
     Throws : no exceptions
   Comments : none 
   See Also : interval_window, interval_window_2, outside_window

=cut

sub inside_window2 {
    my ( $self, $comparable_set, $internal_start, $internal_end, $window_size,
        $maximal_distance )
        = @_;

    # if undefined, set to default values
    $window_size      ||= $self->sw_size;
    $maximal_distance ||= $self->max_in_distance;

    my $comparable_number = $comparable_set->cardinality;

    my @inside_windows;

    foreach my $sw_type (qw/l r/) {

        # $sw_start and $sw_end are both index of $target_set
        my ( $sw_start, $sw_end );

        my $working_set
            = $comparable_set->intersect("$internal_start-$internal_end");
        my $working_length = $working_set->cardinality;
        last if $working_length < $window_size;

        # the windows start from the center
        if ( $sw_type eq 'r' ) {
            $sw_start = int( $working_length / 2 ) + 1;
            $sw_end   = $sw_start + $window_size - 1;
        }
        elsif ( $sw_type eq 'l' ) {
            $sw_end   = int( $working_length / 2 );
            $sw_start = $sw_end - $window_size + 1;
        }

        my $available_distance
            = int( $working_length / ( $window_size * 2 ) );
        my $max_distance = min( $available_distance, $maximal_distance );

        # sw_distance is from -90 to -90 + max_distance - 1
    INSIDESW2: foreach my $i ( -90 .. ( -90 + $max_distance - 1 ) ) {
            my $sw_set = $working_set->slice( $sw_start, $sw_end );
            my $sw_set_member_number = $sw_set->cardinality;
            if ( $sw_set_member_number < $window_size ) {
                last INSIDESW2;
            }
            my $sw_distance = $i;

            my %window_info;
            $window_info{type}     = $sw_type;
            $window_info{set}      = $sw_set;
            $window_info{distance} = $sw_distance;

            push @inside_windows, \%window_info;

            if ( $sw_type eq 'r' ) {
                $sw_start = $sw_end + 1;
                $sw_end   = $sw_start + $window_size - 1;
            }
            elsif ( $sw_type eq 'l' ) {
                $sw_end   = $sw_start - 1;
                $sw_start = $sw_end - $window_size + 1;
            }
        }
    }

    return @inside_windows;
}

1;    # Magic true value required at end of module

__END__

=head1 DESCRIPTION

AlignDB::Window is a part of AlignDB package. It provides methods to split
integer spans, including interval, inside area and outside area, into a series
of windows.

For more information, see L<AlignDB>.
