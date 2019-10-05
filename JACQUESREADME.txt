- Run any of the new examples uploaded (the ones not numbered, but named)
- Look at the "orthCMAMBFs" function in the "CMAMBFinterface" folder, unsure of what to use for orthonormalization. The 'U' and the "orth" function gives different amount of MBFs
- Look at "runCMA_MBFgenerator" function in the "CMAMBFinterface" folder, I cannot identify any particular problems.
- I believe the "runCMACBFM" and "runCMADGFMsolver" functions are correct, but it wouldn't hurt to give them a scan.

At the moment in the examples, only CMA is activated

You'll see the results looking more or less like this:

Using 18 MBFs there is an error norm of 1.261240 percent compared to FEKO sol.
Using 19 MBFs there is an error norm of 1.260086 percent compared to FEKO sol.
Using 20 MBFs there is an error norm of 1.230537 percent compared to FEKO sol.
Using 21 MBFs there is an error norm of 1.135204 percent compared to FEKO sol.
Using 22 MBFs there is an error norm of 1.143585 percent compared to FEKO sol.
Using 23 MBFs there is an error norm of 1.125797 percent compared to FEKO sol.
Using 24 MBFs there is an error norm of 1.123772 percent compared to FEKO sol.
Using 25 MBFs there is an error norm of 0.542017 percent compared to FEKO sol.
Using 26 MBFs there is an error norm of 0.534458 percent compared to FEKO sol.
Using 27 MBFs there is an error norm of 0.000000 percent compared to FEKO sol.


