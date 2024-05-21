This is a small test programme for the 1d Boyd filter.
This programme is intended to filtered a 1d function elementwisely. The field on collocation points are transfered into polynomial space and truncation (filtering) is performed on higher polynomial orders. (see test_filter%trnsfr(test_filter%nx-0) = 0.0 etc.)

Boyd filter is used in Dynamic Smagorinsky model (or any other dynamic models). This programme shows the effects of such filter.
Other elementwise filter based on filtering in polynomial space could be also test based on this programme for a basic 1d case.