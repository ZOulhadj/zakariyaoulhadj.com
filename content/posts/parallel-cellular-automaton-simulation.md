+++
title = "Parallel Cellular Automaton Simulation"
author = ["Zakariya Oulhadj"]
publishDate = 2024-11-27T00:00:00+00:00
tags = ["school", "c", "mpi"]
draft = false
+++

As part of the Message-Passing Programming ([EPCC11002](http://www.drps.ed.ac.uk/24-25/dpt/cxepcc11002.htm)) course, I developed a
parallel implementation of a 2D-decomposed cellular automaton with periodic
boundary conditions on the \\(i^{th} \\) dimension. The boundary conditions for \\(
\frac{2}{3} \\)​ of the \\( j^{th} \\) dimension are set to alive cells. A
termination condition is imposed for the simulation in which the program
terminates if the number of living cells is below 3/4 or greater than 4/3 of the
initial living cells. The implementation uses MPI and a cartesian virtual
topology to decompose the grid into two dimensions where each process receives a
subsection of the grid. Communication between processes is performed using
halo-swapping via non-blocking point-to-point communication (MPI_Isend and
MPI_Irecv).

{{< figure src="/img/posts/mpp-cellular-automaton.png" >}}
