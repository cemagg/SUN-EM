## Change Log
All notable changes to this project will be documented in this file. 

Note the following list of abbreviations:
	WIP : Work in Progress

===================================================================================================

## 1.27 - 2018-08-08
- Added now also PEC fine mesh example in Example-10

## 1.26 - 2018-08-01
- Added now a suggestion for calling the sequential and parallel C++ binaries

## 1.25 - 2018-07-31
- Added now a new example (PEC plate, example 10) for a skripsie student to start working on C++
  Z matrix filling. Also added some place-holders for where this Z-matrix C++ filling should be called.

## 1.24 - 2018-06-20
- Working on Adaptive MBFs -- vivaldi array (example-4) not working.

## 1.23 - 2018-06-20
- Removed now also windowing method from adaptive MBF approach (convergence achieved now for 6x1 bow-tie array
  see example-6).

## 1.22 - 2018-06-20
- Removed now windowing method from iterative Jacobi approach (convergence achieved now for 6x1 bow-tie array
  see example-6).
- TO-DO: Expand this now for the CBFM-enhanced Jacobi (adaptive MBFs, i.e. IFBalg 14)

## 1.21 - 2018-06-14
- Added linear 10x1 bow-tie array (results used for ICEAA'18)

## 1.20 - 2018-06-14
- Added now DGFM (+ iDGFM). Generated sufficient results for linear array (example-8).

## 1.19 - 2018-06-14
- Added now Equivalent Dipole Method (EDM) for accelerated Z matrix filling. Still a bit of tweaking to
  do, but at least the basic feature is there (threshold distance quite a parameter). Need to add now 
  DGFM (+ i-DGFM support).

.0.,
## 1.18 - 2018-06-13
- Removed some debug info and also added a relativeErrorNorm for a matrix (the latter is showing a fairly
  large error for the PEC plate example ~ 42%, even thought the curernt solution is ~ 10.6%)

## 1.17 - 2018-06-13
- I think I fixed the issue now wrt the Z calculation (edge fill) using [DBD2011]. The 180 deg. phase
  shift issue is a result of the Amn+ and Amn- (and Phi_mn+ and Phi_mn-) sign allocation in the routine
  FillZMatrixByEdge. The reason that I think the signs should be reversed is that the sign convention
  stemming from [RWG82, Eq. 24] is not interpreted correctly in the former reference. Results are 
  promising - for example-7/ the rectangular PEC plate mesh results in a rel. error norm % of about
  10.64% (with SING=true) and about 12% with SING=false (i.e. without singularity treatment).

## 1.16 - 2018-06-12
- Some more progress wrt the SUN-EM MoM solver. Very close now. A -1 factor still out. The implementation
  in [DBD2011] is not correct. This -1 sign issue has to do with the magnetic vector potential calculation
  most likely.

## 1.15 - 2018-06-11
- Added some improvements for our internal MoM solver (adapted [DBD2011] implementation to work with 
the Solver_setup data structures.) Not tested yet.

## 1.14 - 2018-06-09
- Added now the dBMoMinterface, i.e. essentially the 3D MoM implementation by David B. Davidson,
  detailed in Chapter 6 of his book, Computational Electromagnetics for RF and Microwave Engineering,
  2nd Edition.

## 1.13 - 2018-06-07
- Promoted now another example for the iterative Jacobi method (example-6, a larger bow-tie array).

## 1.12 - 2018-06-07
- WIP: Primary MBFs + windowing is now working correctly. For a 5 x 1 bow-tie array (example-3), 
  I get now a relative error of 0.82%.
- Continue now first with the Jacobi iterations.

## 1.11 - 2018-06-07
- WIP: Generated now the primary MBFs on the generating sub-arrays (type 1,2 and 3). We need to
  map them now to the correct array elements. Then we can test at least the primary MBFs for the 
  CBFM. For the Jacobi, we actually do not need the secondary effects using the individual primaries, 
  as we only work with a single primary. We can perhaps do a CBFM solution for the 0th current
  distribution (i.e. to get the starting vector a bit more promising).

## 1.10 - 2018-06-07
- WIP: Interconnected domains - sub-arrays for generating better primary MBFs (plotting for
  debug purposes now sorted out). Starting work in runMBFgenerator for the sub-arrays.

## 1.9 - 2018-06-07
- WIP: Interconnected domains - sub-arrays for generating better primary MBFs

## 1.8 - 2018-06-07
- Backing up before adding the primary generating sub-arrays for the CBFM.

## 1.7 - 2018-06-06
- IFBMoMSolver now working with the connected domains (Jacobi - Alg. 7).

## 1.6 - 2018-06-06
- IFBMoMSolver now working with the Adaptive MBF solver again (Alg. 14)
  for disconnected array elements.

## 1.5 - 2018-06-05
- Work in progress: Got now the IFBMoMsolver + JackitSolver going for
  disconnected domains at least.

## 1.4 - 2018-06-05
- Work in progress: Added now example-3 (radiating bow-tie array example).

## 1.3 - 2018-06-05
- Something not quite right still ...

## 1.2 - 2018-06-04
- Working on interconnected CBFM solver. Promote before biggish refactoring.

## 1.1 - 2018-06-01
- Did quite a lot of work on interconnected domains (parsing *.out file for geometry
  and basis function setup).
- TO-DO: Still working on treatment of interconnected structures in CBFM solver.

## 1.0 - 2018-05-31
- First promote - took some of the existing work done for /feko-connect and /harp
  (excluding solver details from the latter)
