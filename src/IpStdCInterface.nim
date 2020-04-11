## ************************************************************************
##    Copyright (C) 2004, 2010 International Business Machines and others.
##    All Rights Reserved.
##    This code is published under the Eclipse Public License.
##
##    $Id: IpStdCInterface.h 2082 2012-02-16 03:00:34Z andreasw $
##
##    Authors:  Carl Laird, Andreas Waechter     IBM    2004-09-02
## ***********************************************************************

{.deadCodeElim: on.}
when defined(windows):
  const
    #libName* = "lib/libipopt.lib"
    libName* = "lib/win64_dll/libipopt-1.dll"  # Downloaded from: https://github.com/JuliaOpt/IpoptBuilder/releases
elif defined(macosx):
  const
    libName* = "libsoloud.dylib"
else:
  const
    libName* = "libipopt.so"
## * Type for all number.  We need to make sure that this is
##       identical with what is defined in Common/IpTypes.hpp
#{.link: libName .}



## * enum to indicate the mode in which the algorithm is

type
  AlgorithmMode* = enum
    RegularMode = 0, RestorationPhaseMode = 1

type
  Number* = cdouble

## * Type for all incides.  We need to make sure that this is
##       identical with what is defined in Common/IpTypes.hpp

type
  Index* = cint

## * Type for all integers.  We need to make sure that this is
##       identical with what is defined in Common/IpTypes.hpp

type
  Int* = cint

##  This includes the SolverReturn enum type

type
  ApplicationReturnStatus* = enum
    Internal_Error = -199, Insufficient_Memory = -102,
    NonIpopt_Exception_Thrown = -101, Unrecoverable_Exception = -100,
    Invalid_Number_Detected = -13, Invalid_Option = -12,
    Invalid_Problem_Definition = -11, Not_Enough_Degrees_Of_Freedom = -10,
    Maximum_CpuTime_Exceeded = -4, Error_In_Step_Computation = -3,
    Restoration_Failed = -2, Maximum_Iterations_Exceeded = -1, Solve_Succeeded = 0,
    Solved_To_Acceptable_Level = 1, Infeasible_Problem_Detected = 2,
    Search_Direction_Becomes_Too_Small = 3, Diverging_Iterates = 4,
    User_Requested_Stop = 5, Feasible_Point_Found = 6

## * Structure collecting all information about the problem
##   definition and solve statistics etc.  This is defined in the
##   source file.

type
  IpoptProblemInfo* {.bycopy.} = object


## * Pointer to a Ipopt Problem.

type
  IpoptProblem* = ptr IpoptProblemInfo

## * define a boolean type for C

type
  Bool* = cint

const
  TRUE* :Bool = 1
  FALSE*:Bool = 0

## * A pointer for anything that is to be passed between the called
##   and individual callback function

type
  UserDataPtr* = pointer

## * Type defining the callback function for evaluating the value of
##   the objective function.  Return value should be set to false if
##   there was a problem doing the evaluation.

type
  Eval_F_CB* = proc (n: Index; x: ptr Number; new_x: Bool; obj_value: ptr Number;
                  user_data: UserDataPtr): Bool {.cdecl.}

## * Type defining the callback function for evaluating the gradient of
##   the objective function.  Return value should be set to false if
##   there was a problem doing the evaluation.

type
  Eval_Grad_F_CB* = proc (n: Index; x: ptr Number; new_x: Bool; grad_f: ptr Number;
                       user_data: UserDataPtr): Bool {.cdecl.}

## * Type defining the callback function for evaluating the value of
##   the constraint functions.  Return value should be set to false if
##   there was a problem doing the evaluation.

type
  Eval_G_CB* = proc (n: Index; x: ptr Number; new_x: Bool; m: Index; g: ptr Number;
                  user_data: UserDataPtr): Bool {.cdecl.}

## * Type defining the callback function for evaluating the Jacobian of
##   the constrant functions.  Return value should be set to false if
##   there was a problem doing the evaluation.

type
  Eval_Jac_G_CB* = proc (n: Index; x: ptr Number; new_x: Bool; m: Index; nele_jac: Index;
                      iRow: ptr Index; jCol: ptr Index; values: ptr Number;
                      user_data: UserDataPtr): Bool {.cdecl.}

## * Type defining the callback function for evaluating the Hessian of
##   the Lagrangian function.  Return value should be set to false if
##   there was a problem doing the evaluation.

type
  Eval_H_CB* = proc (n: Index; x: ptr Number; new_x: Bool; obj_factor: Number; m: Index;
                  lambda: ptr Number; new_lambda: Bool; nele_hess: Index;
                  iRow: ptr Index; jCol: ptr Index; values: ptr Number;
                  user_data: UserDataPtr): Bool {.cdecl.}

## * Type defining the callback function for giving intermediate
##   execution control to the user.  If set, it is called once per
##   iteration, providing the user with some information on the state
##   of the optimization.  This can be used to print some
##   user-defined output.  It also gives the user a way to terminate
##   the optimization prematurely.  If this method returns false,
##   Ipopt will terminate the optimization.

type
  Intermediate_CB* = proc (alg_mod: Index; iter_count: Index; obj_value: Number;
                        inf_pr: Number; inf_du: Number; mu: Number; d_norm: Number;
                        regularization_size: Number; alpha_du: Number;
                        alpha_pr: Number; ls_trials: Index; user_data: UserDataPtr): Bool {.
      cdecl.}                 ##  0 is regular, 1 is resto

## * Function for creating a new Ipopt Problem object.  This function
##   returns an object that can be passed to the IpoptSolve call.  It
##   contains the basic definition of the optimization problem, such
##   as number of variables and constraints, bounds on variables and
##   constraints, information about the derivatives, and the callback
##   function for the computation of the optimization problem
##   functions and derivatives.  During this call, the options file
##   PARAMS.DAT is read as well.
##
##   If NULL is returned, there was a problem with one of the inputs
##   or reading the options file.

proc CreateIpoptProblem*(n: Index; ## * Number of optimization variables
                        x_L: ptr Number; ## * Lower bounds on variables. This array of
                                      ##                               size n is copied internally, so that the
                                      ##                               caller can change the incoming data after
                                      ##                               return without that IpoptProblem is
                                      ##                               modified.  Any value less or equal than
                                      ##                               the number specified by option
                                      ##                               'nlp_lower_bound_inf' is interpreted to
                                      ##                               be minus infinity.
                        x_U: ptr Number; ## * Upper bounds on variables. This array of
                                      ##                               size n is copied internally, so that the
                                      ##                               caller can change the incoming data after
                                      ##                               return without that IpoptProblem is
                                      ##                               modified.  Any value greater or equal
                                      ##                               than the number specified by option
                                      ##                               'nlp_upper_bound_inf' is interpreted to
                                      ##                               be plus infinity.
                        m: Index; ## * Number of constraints.
                        g_L: ptr Number; ## * Lower bounds on constraints. This array of
                                      ##                               size m is copied internally, so that the
                                      ##                               caller can change the incoming data after
                                      ##                               return without that IpoptProblem is
                                      ##                               modified.  Any value less or equal than
                                      ##                               the number specified by option
                                      ##                               'nlp_lower_bound_inf' is interpreted to
                                      ##                               be minus infinity.
                        g_U: ptr Number; ## * Upper bounds on constraints. This array of
                                      ##                               size m is copied internally, so that the
                                      ##                               caller can change the incoming data after
                                      ##                               return without that IpoptProblem is
                                      ##                               modified.  Any value greater or equal
                                      ##                               than the number specified by option
                                      ##                               'nlp_upper_bound_inf' is interpreted to
                                      ##                               be plus infinity.
                        nele_jac: Index; ## * Number of non-zero elements in constraint
                                       ##                               Jacobian.
                        nele_hess: Index; ## * Number of non-zero elements in Hessian of
                                        ##                               Lagrangian.
                        index_style: Index; ## * indexing style for iRow & jCol,
                                          ## 				 0 for C style, 1 for Fortran style
                        eval_f: Eval_F_CB; ## * Callback function for evaluating
                                         ##                               objective function
                        eval_g: Eval_G_CB; ## * Callback function for evaluating
                                         ##                               constraint functions
                        eval_grad_f: Eval_Grad_F_CB; ## * Callback function for evaluating gradient
                                                   ##                               of objective function
                        eval_jac_g: Eval_Jac_G_CB; ## * Callback function for evaluating Jacobian
                                                 ##                               of constraint functions
                        eval_h: Eval_H_CB): IpoptProblem {.cdecl,
    importc: "CreateIpoptProblem", dynlib: libName.}
  ## * Callback function for evaluating Hessian
  ##                               of Lagrangian function
## * Method for freeing a previously created IpoptProblem.  After
##       freeing an IpoptProblem, it cannot be used anymore.

proc FreeIpoptProblem*(ipopt_problem: IpoptProblem) {.cdecl,
    importc: "FreeIpoptProblem", dynlib: libName.}
## * Function for adding a string option.  Returns FALSE the option
##   could not be set (e.g., if keyword is unknown)

proc AddIpoptStrOption*(ipopt_problem: IpoptProblem; keyword: cstring; val: cstring): Bool {.
    cdecl, importc: "AddIpoptStrOption", dynlib: libName.}
## * Function for adding a Number option.  Returns FALSE the option
##   could not be set (e.g., if keyword is unknown)

proc AddIpoptNumOption*(ipopt_problem: IpoptProblem; keyword: cstring; val: Number): Bool {.
    cdecl, importc: "AddIpoptNumOption", dynlib: libName.}
## * Function for adding an Int option.  Returns FALSE the option
##   could not be set (e.g., if keyword is unknown)

proc AddIpoptIntOption*(ipopt_problem: IpoptProblem; keyword: cstring; val: Int): Bool {.
    cdecl, importc: "AddIpoptIntOption", dynlib: libName.}
## * Function for opening an output file for a given name with given
##   printlevel.  Returns false, if there was a problem opening the
##   file.

proc OpenIpoptOutputFile*(ipopt_problem: IpoptProblem; file_name: cstring;
                         print_level: Int): Bool {.cdecl,
    importc: "OpenIpoptOutputFile", dynlib: libName.}
## * Optional function for setting scaling parameter for the NLP.
##   This corresponds to the get_scaling_parameters method in TNLP.
##   If the pointers x_scaling or g_scaling are NULL, then no scaling
##   for x resp. g is done.

proc SetIpoptProblemScaling*(ipopt_problem: IpoptProblem; obj_scaling: Number;
                            x_scaling: ptr Number; g_scaling: ptr Number): Bool {.
    cdecl, importc: "SetIpoptProblemScaling", dynlib: libName.}
## * Setting a callback function for the "intermediate callback"
##   method in the TNLP.  This gives control back to the user once
##   per iteration.  If set, it provides the user with some
##   information on the state of the optimization.  This can be used
##   to print some user-defined output.  It also gives the user a way
##   to terminate the optimization prematurely.  If the callback
##   method returns false, Ipopt will terminate the optimization.
##   Calling this set method to set the CB pointer to NULL disables
##   the intermediate callback functionality.

proc SetIntermediateCallback*(ipopt_problem: IpoptProblem;
                             intermediate_cb: Intermediate_CB): Bool {.cdecl,
    importc: "SetIntermediateCallback", dynlib: libName.}
## * Function calling the Ipopt optimization algorithm for a problem
##       previously defined with CreateIpoptProblem.  The return
##       specified outcome of the optimization procedure (e.g., success,
##       failure etc).
##

proc IpoptSolve*(ipopt_problem: IpoptProblem; 
                x: ptr Number;
                g: ptr Number;
                obj_val: ptr Number; 
                mult_g: ptr Number; 
                mult_x_L: ptr Number;
                mult_x_U: ptr Number;
                user_data: UserDataPtr): ApplicationReturnStatus {.cdecl,
    importc: "IpoptSolve", dynlib: libName.}
  ##[
  - `ipopt_problem`: Problem that is to be optimized.  Ipopt will use the options previously specified with AddIpoptOption (etc) for this problem.
  
  - `x`: 

    - Input:  Starting point
    - Output: Optimal solution

  - `g`: Values of constraint at final point (output only - ignored if set to NULL)
  - `obj_val`: Final value of objective function (output only - ignored if set to NULL)

  - `mult_g`:

    - Input: Initial values for the constraint multipliers (only if warm start option is chosen)
    - Output: Final multipliers for constraints (ignored if set to NULL)

  - `mult_x_L`:
 
    - Input: Initial values for the multipliers for lower variable bounds (only if warm start option is chosen)
    - Output: Final multipliers for lower variable bounds (ignored if set to NULL)

  - `mult_x_U`: 

    - Input: Initial values for the multipliers for upper variable bounds (only if warm start opon is chosen)  
    - Output: Final multipliers for upper variable  bounds (ignored if set to NULL)

  - `user_data`: Pointer to user data.  This will be passed unmodified to the callback functions.
  ]##