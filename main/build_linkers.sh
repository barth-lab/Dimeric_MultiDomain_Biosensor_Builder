    $R/source/bin/loopmodel.linuxgccdebug \
        -database $R/database/  \
        -in:file:s $I/c.0.0.pdb \
        -loop_file $I/loopfile \
        -loops:frag_sizes 9 3 \
        -loops:frag_files input_scaffold${r}/frags.200.9mers input_scaffold${r}/frags.200.3mers \
        -loops:remodel quick_ccd \
        -loops:refine refine_ccd \
        -loops:relax no \
        -constant_seed \
        -nstruct 1 \
        -out:suffix _loop \
        -out:path:all $O \
        -out:pdb \
        -ignore_zero_occupancy false \
        > $O/log \
        2> $O/err

