## Change Log
All notable changes to this project will be documented in this file. 

Note the following list of abbreviations:
	WIP : Work in Progress

===================================================================================================

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
