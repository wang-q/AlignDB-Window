use strict;
use warnings;

use Test::More;
use AlignDB::IntSpan;

use AlignDB::Window;

{
    print "#interval_window\n";

    my $maker = AlignDB::Window->new;

    my @data = (
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
    );

    for my $i ( 0 .. $#data ) {
        my ( $input_ref, $except_ref ) = @{ $data[$i] };

        my @results = $maker->interval_window( @{$input_ref} );
        $_->{set} = $_->{set}->runlist for @results;
#        print YAML::Syck::Dump \@results;
        is_deeply( \@results, $except_ref, "interval window $i" );
    }
}

done_testing();
