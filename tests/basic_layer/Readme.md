# basic_layer tests

These tests address the basic functionality of MASKE, with only particle deletion and nucleation. The examples all assume nearly-fixed solution with nearly-constant concentrations of ions, thus constant saturation index $\beta$

They all refer to a thick layer of calcium carbonate, dissolving or growing depending on $\beta$. 

## Layout

| File or folder | Description                |
|--------|----------------------------|
| benchmarking.xlsx    | Support MS Excel spreadsheet to easily plot results               |
| prepare_inputs.xlsx    | MS Excel spreadsheet to help prepare inputs if needed       |
| Test_details.pdf  | A pdf with explanation of inputs    |
| sim_net    | Simulations with net rates        |
| sim_ki0.5 | Simulations with straight rates and energy penalty scaled by  *$\ki=0.5$*     |


## Running the tests

The first step is to obtain simulation results, running the tests in the *sim_net* and *sim_ki0.5* folders. Then you can compare your results with ours, either by looking directly at the *thermo* files, or using the *benchmarking.xlsx* spreadsheet provided. The other files, *prepare_inputs.xlsx* and *Test_details.pdf* are to support you in case you want to build furhter on the example; their content and usage is not explained here.

### sim_net

Here you can run two simulations for benchmarking: one with $\beta = 0.1$, causing dissolution, the other one with $\beta = 1.1$, causing precipitation.

The folder contains the following files:
* *chemDB.dat* : MASKE's chemsitry database
* *input.dat* : MASKE's input file
* *thermo_b0.1_20210205.txt* : thermo output for $\beta = 0.1$
* *thermo_b1.1_20210205.txt* : thermo output for $\beta = 1.1$

To set $\beta=0.1$, open ``input.dat`` , find the command ``sol_start`` and make sure that only the  ``sol_start`` under the title ``# beta = 0.1`` is not commented, viz. remove the # symbol from before ``sol_start``. Do not remove the # from ``# beta = 0.1``; that on is just a title. Do not leave other ``sol_start`` command lines uncommented, or that might cause the simulation to run with a different $\beta$.

Run the simulation from the command line. The default option in ``input.dat`` is to use 2 processors, thus use:

```
mpiexec -np 2 ./maske input.dat
```

If you want to run on a single processor, you need to change the line in ``input.dat`` reading as follows:

```
subcomm 1 Antonello 2 lmp_yes 123
```

replacing the 2 with a 1 after *Antonello*, then running mpiexec with ``-np 1`` .

You can then compare the results in the *thermo.txt* produced by the simulation with the results in our thermo file for benchmarking. You can also copy-paste the content of your *thermo.txt* into the provided *benchmarking.xlsx* to have the results plotted against the benchmark.

The procedure to run the simulation with $\beta = 1.1$ is the same, except that in ``input.dat`` you must uncomment the ``sol_start`` command corresponding to the case ``# beta = 1.1`` and make sure every other ``sol_start`` is commented out. 

These simulations should be quite fast (< 1 minute).

If you want to see a video of the output, the trajectory is saves as LAMMPS dump file in ``dump.all.Antonello``.

The other output files are not explained here.


### sim_ki0.5

Follow the same steps as for the *sim_net* examples. 

In this case, the simulations can be lengthy (a few minutes). I suggest you stop them after a few hundreds of steps, by pressing   ``ctrl c`` in the terminal. You can see how many steps have been completed from the first column in  *thermo.txt*, while it is being written during the simulation.
