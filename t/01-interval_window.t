use strict;
use warnings;

use Test::More;
use AlignDB::IntSpan;

use AlignDB::Window;

{
    print "#interval_window\n";

    my $maker = AlignDB::Window->new;

    my @data = (
        [ [ AlignDB::IntSpan->new->add_pair( 1, 10 ), 1, 99 ], [], ],
        [   [ AlignDB::IntSpan->new->add_pair( 1, 99 ), 1, 99 ],
            [   {   density  => -1,
                    distance => -1,
                    set      => "1-99",
                    type     => "S",
                },
            ],
        ],
        [   [ AlignDB::IntSpan->new->add_pair( 1, 99 ), 1, 9999 ],
            [   {   density  => -1,
                    distance => -1,
                    set      => "1-99",
                    type     => "S",
                },
            ],
        ],
        [   [ AlignDB::IntSpan->new->add_pair( 1, 100 ), 1, 9999 ],
            [   {   density  => 0,
                    distance => 0,
                    set      => "1-50",
                    type     => "L",
                },
                {   density  => 0,
                    distance => 0,
                    set      => "51-100",
                    type     => "R",
                },
            ],
        ],
        [   [ AlignDB::IntSpan->new->add_pair( 1, 111 ), 1, 9999 ],
            [   {   density  => 0,
                    distance => 0,
                    set      => "1-56",
                    type     => "L",
                },
                {   density  => 0,
                    distance => 0,
                    set      => "57-111",
                    type     => "R",
                },
            ],
        ],
        [   [ AlignDB::IntSpan->new->add_pair( 1, 200 ), 1, 9999 ],
            [   {   density  => 1,
                    distance => 0,
                    set      => "1-50",
                    type     => "L",
                },
                {   density  => 1,
                    distance => 0,
                    set      => "151-200",
                    type     => "R",
                },
                {   density  => 1,
                    distance => 1,
                    set      => "51-150",
                    type     => "L",
                },
            ],
        ],
        [   [ AlignDB::IntSpan->new->add_pair( 1, 605 ), 1, 9999 ],
            [   {   density  => 5,
                    distance => 0,
                    set      => "1-50",
                    type     => "L",
                },
                {   density  => 5,
                    distance => 0,
                    set      => "556-605",
                    type     => "R",
                },
                {   density  => 5,
                    distance => 1,
                    set      => "51-150",
                    type     => "L",
                },
                {   density  => 5,
                    distance => 1,
                    set      => "456-555",
                    type     => "R",
                },
                {   density  => 5,
                    distance => 2,
                    set      => "151-250",
                    type     => "L",
                },
                {   density  => 5,
                    distance => 2,
                    set      => "356-455",
                    type     => "R",
                },
                {   density  => 5,
                    distance => 3,
                    set      => "251-355",
                    type     => "L",
                },
            ],
        ],
    );

    for my $i ( 0 .. $#data ) {
        my ( $input_ref, $except_ref ) = @{ $data[$i] };

        my @results = $maker->interval_window( @{$input_ref} );
        $_->{set} = $_->{set}->runlist for @results;
        is_deeply( \@results, $except_ref, "interval window $i" );
    }
}

done_testing();
