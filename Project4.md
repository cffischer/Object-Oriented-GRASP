## Introduction
Project 4 is developing a version of the GRASP2018 `rangular` program that introduces several very important data structures.
* orbitals (ORB)
* configuration state functions (CSF)
* Linked  Lists of Angular data

The `Report-WaveFunctionExpansion` in the present repository comments on various aspects of expansions in the application of GRASP to heavy elements and indicates the cases that should be benchmarked in this work.  Similarly, the `Report-DataStructures` reviews various objects in atomic structure calculations and suggests some user defined data types but they have not been validated in practice and hence are only a starting point in this work.

The objective here is not to change the way angular data is computed but to make the program object oriented.  The first attempt is to make minimal changes and to totally change the way information is stored. In particular, GRASP2018 memory management system for ALLOC, DALLOC, and RALLOC will not be used.    

Project 3 introduced data structures for ORBs  and CSFs but these do not play an essential role in the algorithms generating expansions.  Though a consistent definition is needed for these two projects,  the program RCSFGENERATE_shorten is already available for producing input data for the current project. Thus this project can proceed without the completion of Project 3.

## The Main program
The main program of the current RANGULAR version must remove all the USE statements that refer to variables in COMMON or variables that are replaced in new modules associated with data types or other new modules.
### Module Const

This could be the Zconst module of DBSR_HF  with addition parameters associated with our GRASP.  These include  the following:
  * Kind- constants ( we should define both the GRASP kinds as well as the DBSR_HF kinds but with the intention of ultimately using only the latter)
  * Constants  (used for readability, such as `half` rather than `0.5_dp`)
  * Physical constants (conversion factors for units)
  * Parameters for our codes  (like those in our parameter_def  file and others)
  
Let's start with the zconst version and add to this module as we see fit.  It is my recollection that parameters like factorials are needed in atomic structure. Maybe they should go into this module. It would be used by every program in the GRASP collection.
* Module Atoms_const
This could be our version of the DBSR_HF  module atoms_par but modified to our view of things.  In particular, our new rnucleus module could be included here   (CONTAINS). 

### Library  Angular

Unfortunately the  GRASP Lib92 and RANGULAR are not clearly documented.  

These new modules would be USED In the main program.  Some CALLS can be removed since, at this stage, we do not want to think about debugging and user defaults.  I think we would also move in the DBSR_HF directions and have `atoms.inp` file that has defaults as well as arguments on the command line.  So that means we need to think about a set of basic keywords that are recognized on the command line, but maybe this is not something we should worry about initially. 

### The algorithm

The current version generates all the angular data and then sorts the data into some order.  Sorting can be avoided which is desirable when lists get very large.  A  faster procedure is to generate the list of integrals from the orbital set as in `genintrk.f90` of RCI, except that the 'Value' array is not computed by RANGULAR.  There could by a Logical array to keep track of which integrals actually appear in the final list. This list of Integrals and whether they are used could be stored and passed on the RMCDHF so that the latter knows which integrals need to be evaluated.  So the basic algorithm is:

  * Read an atom.c file that contains the orbitals and their expansion
  * Generate the integral list
  * For each block
      For each column  j
         For each row  i <= j 
            Compute contributions to the data list in the form of [ coefficient, index to integral, I ] if the matrix element is not zero.
            
Also needed is an array indicating where the data for each row ends and another that indicates where the data for each column ends.  Note that the matrix is sparse so there may not be entries for all (i,j) but there always data for i=j.  The data is collected in a `chunk`, size =1000 (a parameter) and the chunks linked together.  Actually, maybe a linked list is not needed if information about size is stored in a separate file ( atom.h header file). This would allow RMCDHF to read the header, allocate the needed memory and read (or leave on disk if not available). In parallel version, the columns are each distributed cyclically to a processor.   So maybe with the new algorithm we do not need to introduce a linked list only for smaller, in memory cases.  The linked list would be needed if RANGULAR and RMCDHF were merged with no writing to disk.  For parallel execution this would mean we have no restart capabilities. Possibly, initially, we could omit the linked list.

