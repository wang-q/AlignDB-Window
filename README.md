# NAME

AlignDB::Window - Split integer spans into a series of windows

# DESCRIPTION

AlignDB::Window provides methods to split integer spans, including interval, inside area and outside area, into a series of windows.

# ATTRIBUTES

## sw\_size

sliding windows' size, default is 100

## min\_interval

mininal indel interval length, default is 11

## max\_out\_distance

maximal outside distance, default is 10

## max\_in\_distance

maximal inside distance, default is 5

# METHODS

## interval\_window

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

## interval\_window\_2

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

## outside\_window

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

## outside\_window\_2

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

## inside\_window

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

## inside\_window2

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

## center\_window

         Usage : $self->outside_window_2(
               :     $comparable_set, $interval_start, $interval_end,
               :     $window_size, $maximal_distance,
               : );
       Purpose : draw windows for a certain region, center is 0, and
               : the first window is 50 bp and all others are 100 bp length
               : start from 0 and end to $maximal_distance
       Returns : @outside_windows
               : each member is a hash_ref
    Parameters : $comparable_set   : AlignDB::IntSpan object
               : $interval_start   : start position of the interval
               : $interval_end     : end position of the interval
               : $window_size      : size of windows
               : $maximal_distance : maximal distance

## center\_intact\_window

         Usage : $self->outside_window_2(
               :     $comparable_set, $interval_start, $interval_end,
               :     $window_size, $maximal_distance,
               : );
       Purpose : draw windows for a certain region, center is 0, and
               : the first window is 50 bp and all others are 100 bp length
               : start from 0 and end to $maximal_distance
       Returns : @outside_windows
               : each member is a hash_ref
    Parameters : $comparable_set   : AlignDB::IntSpan object
               : $interval_start   : start position of the interval
               : $interval_end     : end position of the interval
               : $window_size      : size of windows
               : $maximal_distance : maximal distance

## strand\_window

         Usage : $self->strand_window(
               :     $comparable_set, $interval_start, $interval_end,
               :     $window_size, $strand,
               : );
       Purpose : draw windows for a certain region
       Returns : @windows
               : each member is a hash_ref
    Parameters : $comparable_set   : AlignDB::IntSpan object
               : $interval_start   : start position of the interval
               : $interval_end     : end position of the interval
               : $window_size      : size of windows
               : $strand           : '+' or '-'

# AUTHOR

Qiang Wang &lt;wang-q@outlook.com>

# COPYRIGHT AND LICENSE

This software is copyright (c) 2008 by Qiang Wang.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.
