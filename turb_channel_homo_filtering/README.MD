This is a user .f90 file for filtering in homogeneous directions (i.e., x- and z-direction) in the example of turb_channel.
The intention of this filtering is to smear out the unphysical oscillations close to the edges of elements and potentially contribute to solve the issue of uneven distribution of wavenumbers in SEM method due to the finer resolution close to the edges of elements.
Application simulation: Wall-resolved LES
For wall-modelled LES, one could refer to the Master thesis of Lorenzo Luca Donati about the atmospheric boundary layer simulation, and it turns out the field in such simulations is rather smooth, so it might not need such filtering technique.
Beside, it's worthy discussion about the sequence of 1. time stepping, 2. filtering and 3. field output
