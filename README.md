# MSEsim

MSEsim is an IDL code designed to simulate the preformance of an MSE diagnostic by modelling the emission from a tokamak including the effects of a finite viewing volume and finite collection optics f-number. It was written by Maarten DeBock with support from Euratom, CCFE and TuE. It is released under the MIT Open Source licence.

## Installation & Setup

1. Obtain the files from this repository with the Git command:

  `git clone https://github.com/euratom-software/msesim.git`

1. Add the code diretories to your IDL path before attempting to run any IDL
programs. eg 

  `export IDL_PATH=$IDL_PATH+":+/<your_directory>/msesim"`

1. Change to the `runs` directory and create a subdirectory to work in. You
should give this the name of your tokamak. There is a `test` directory
which you can use as a template. It should contain an output, and a settings folder.

1. Edit the file `batch.xml` to describe what you want the code to do (add
filters, show spectrum filtered or unfiltered, merge the spectrum etc) there's an
example of this in the runs/test folder.

  **_I'm not at all sure that the comments below are accurate_**

  The program needs various files to describe the plasma, the optics and the beam
  geometries. There are some tools to help prepare these and to write them in the
  correct format.

  * Within `msesim/equi` there is a file `create_equig.pro` which when run
  
  `create_equig, <name_of_geqdsk_file>`
  
  will read the given EFIT geqdsk file and write an IDL save file (xdr/sav format) for input to MSEsim

  * The program `coll_optics\create_coll_param.pro` is a useful template for
    creating the optics definition `.xml` file
	
  Beam geometries, beamlet layouts, energy fractions, current waveforms are read
  from an `xml` file in the ```beam``` directory. See ```beam/SouthNBI_2018.xml```
  as an example.

1. Run the code with:

   `IDL`  
   `> msesim, 'name of directory within <runs>'`
   
## Issues
