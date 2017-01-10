use strict;
use warnings;

use Test::More;
use AlignDB::IntSpan;

use AlignDB::Window;

{
    print "#interval_window\n";

    my $maker = AlignDB::Window->new( max_out_distance => 2, );

    my @data = (
        [   [ AlignDB::IntSpan->new->add_pair( 1, 9999 ), 500, 500 ],
            [   {   distance => 1,
                    set      => "400-499",
                    type     => "L",
                },
                {   distance => 2,
                    set      => "300-399",
                    type     => "L",
                },
                {   distance => 1,
                    set      => "501-600",
                    type     => "R",
                },
                {   distance => 2,
                    set      => "601-700",
                    type     => "R",
                },
            ],
        ],
        [   [ AlignDB::IntSpan->new->add_pair( 1, 9999 ), 500, 600 ],
            [   {   distance => 1,
                    set      => "400-499",
                    type     => "L",
                },
                {   distance => 2,
                    set      => "300-399",
                    type     => "L",
                },
                {   distance => 1,
                    set      => "601-700",
                    type     => "R",
                },
                {   distance => 2,
                    set      => "701-800",
                    type     => "R",
                },
            ],
        ],
        [   [ AlignDB::IntSpan->new->add_pair( 1, 9999 ), 101, 101 ],
            [   {   distance => 1,
                    set      => "1-100",
                    type     => "L",
                },
                {   distance => 1,
                    set      => "102-201",
                    type     => "R",
                },
                {   distance => 2,
                    set      => "202-301",
                    type     => "R",
                },
            ],
        ],
    );

    for my $i ( 0 .. $#data ) {
        my ( $input_ref, $except_ref ) = @{ $data[$i] };

        my @results = $maker->outside_window( @{$input_ref} );
        $_->{set} = $_->{set}->runlist for @results;
        is_deeply( \@results, $except_ref, "outside window $i" );
    }
}

done_testing();
