#{.passL:"lib/libipopt.lib".}
import IpStdCInterface
import strformat

# https://github.com/JuliaOpt/Ipopt.jl/blob/master/src/Ipopt.jl
#[
                        eval_grad_f: Eval_Grad_F_CB; ## * 
                        eval_jac_g: Eval_Jac_G_CB; ## * 
                        eval_h: Eval_H_CB)
]#


# This is an example how user_data can be used. 
type 
   MyUserData {.bycopy.} = object
     g_offset:array[2,Number]    # This is an offset for the constraints.
   #UserDataPtr = ptr MyUserData
#[
Bool intermediate_cb(
   Index       alg_mod,
   Index       iter_count,
   Number      obj_value,
   Number      inf_pr,
   Number      inf_du,
   Number      mu,
   Number      d_norm,
   Number      regularization_size,
   Number      alpha_du,
   Number      alpha_pr,
   Index       ls_trials,
   UserDataPtr user_data
);


]#
proc eval_f( n:Index,
             x:ptr Number, 
             new_x:Bool,
             obj_value:ptr Number, 
             user_data:UserDataPtr):Bool  {.cdecl,exportc.} =
  ## Callback function for evaluating objective function
  assert(n == 4)
  let buf = cast[ptr UncheckedArray[Number]](x) 
  
  let tmp = buf[0]*buf[3]*(buf[0]+buf[1]+buf[2]) + buf[2]
  obj_value[] = tmp
  return TRUE


proc eval_grad_f(n: Index; 
                 xx: ptr Number; # Input
                 new_x: Bool; 
                 grad_ff: ptr Number; # Output
                 user_data: UserDataPtr): Bool {.cdecl, exportc.} =
  ## Callback function for evaluating gradient of objective function
  assert(n == 4)
  let x     = cast[ptr UncheckedArray[Number]](xx)
  let grad_f = cast[ptr UncheckedArray[Number]](grad_ff)
  grad_f[0] = x[0] * x[3] + x[3] * (x[0] + x[1] + x[2])
  grad_f[1] = x[0] * x[3]
  grad_f[2] = x[0] * x[3] + 1.0
  grad_f[3] = x[0] * (x[0] + x[1] + x[2])
  return TRUE


proc eval_g( n: Index; 
             xx: ptr Number; 
             new_x: Bool; 
             m: Index; 
             gg: ptr Number; 
             user_data: UserDataPtr): Bool  {.cdecl, exportc.} =
  ##Callback function for evaluating constraint functions
  assert(n == 4)
  assert(m == 2)
  
  let my_data = cast[ptr MyUserData](user_data)
  let x       = cast[ptr UncheckedArray[Number]](xx)
  let g       = cast[ptr UncheckedArray[Number]](gg)

  g[0] = x[0] * x[1] * x[2] * x[3] + my_data.g_offset[0]
  g[1] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3] + my_data.g_offset[1]
  return TRUE


proc eval_jac_g( n: Index; 
                 xx: ptr Number; 
                 new_x: Bool; 
                 m: Index; 
                 nele_jac: Index;
                 iiRow: ptr Index; 
                 jjCol: ptr Index; 
                 vvalues: ptr Number;
                 user_data: UserDataPtr): Bool {.cdecl, exportc.} =
  ## Callback function for evaluating Jacobian of constraint functions
  if vvalues == nil:
    # return the structure of the jacobian 
    # this particular jacobian is dense 
    let iRow = cast[ptr UncheckedArray[Number]](iiRow)
    let jCol = cast[ptr UncheckedArray[Number]](jjCol)
    iRow[0] = 0.0
    jCol[0] = 0.0
    iRow[1] = 0.0
    jCol[1] = 1.0
    iRow[2] = 0.0
    jCol[2] = 2.0
    iRow[3] = 0.0
    jCol[3] = 3.0
    iRow[4] = 1.0
    jCol[4] = 0.0
    iRow[5] = 1.0
    jCol[5] = 1.0
    iRow[6] = 1.0
    jCol[6] = 2.0
    iRow[7] = 1.0
    jCol[7] = 3.0
  else:
    # return the values of the jacobian of the constraints 
    let x      = cast[ptr UncheckedArray[Number]](xx)
    let values = cast[ptr UncheckedArray[Number]](vvalues)
    values[0] = x[1] * x[2] * x[3]   # 0,0
    values[1] = x[0] * x[2] * x[3]   # 0,1
    values[2] = x[0] * x[1] * x[3]   # 0,2 
    values[3] = x[0] * x[1] * x[2]   # 0,3

    values[4] = 2.0 * x[0]             # 1,0 
    values[5] = 2.0 * x[1]             # 1,1 
    values[6] = 2.0 * x[2]             # 1,2 
    values[7] = 2.0 * x[3]             # 1,3 
  return TRUE


proc eval_h( n: Index; 
             xx: ptr Number; 
             new_x: Bool; 
             obj_factor: Number; 
             m: Index;
             llambda: ptr Number; 
             new_lambda: Bool; 
             nele_hess: Index;
             iiRow: ptr Index; 
             jjCol: ptr Index; 
             vvalues: ptr Number;
             user_data: UserDataPtr): Bool {.cdecl, exportc.} =
  
  var idx:Index = 0 # nonzero element counter
  var row:Index = 0 # row counter for loop 
  var col:Index = 0 # col counter for loop

  if vvalues == nil:
    # return the structure. This is a symmetric matrix, fill the lower left triangle only. 
    # the hessian for this problem is actually dense
    let iRow = cast[ptr UncheckedArray[Number]](iiRow)
    let jCol = cast[ptr UncheckedArray[Number]](jjCol)
    for row in 0..<4:
      for col in 0..row:
        iRow[idx] = row.float
        jCol[idx] = col.float
        idx += 1
    assert(idx == nele_hess) 

  else:
    # return the values. This is a symmetric matrix, fill the lower left triangle only
    
    # fill the objective portion
    let x      = cast[ptr UncheckedArray[Number]](xx)
    let values = cast[ptr UncheckedArray[Number]](vvalues)

    values[0] = obj_factor * (2 * x[3])  # 0,0 

    values[1] = obj_factor * (x[3])      # 1,0 
    values[2] = 0                        # 1,1 

    values[3] = obj_factor * (x[3])      # 2,0 
    values[4] = 0                        # 2,1 
    values[5] = 0                        # 2,2 

    values[6] = obj_factor * (2 * x[0] + x[1] + x[2])   # 3,0 
    values[7] = obj_factor * (x[0])                     # 3,1 
    values[8] = obj_factor * (x[0])                     # 3,2 
    values[9] = 0                                       # 3,3 

    # add the portion for the first constraint
    let lambda = cast[ptr UncheckedArray[Number]](llambda)
    values[1] += lambda[0] * (x[2] * x[3])   # 1,0 

    values[3] += lambda[0] * (x[1] * x[3])   # 2,0 
    values[4] += lambda[0] * (x[0] * x[3])   # 2,1 

    values[6] += lambda[0] * (x[1] * x[2])   # 3,0 
    values[7] += lambda[0] * (x[0] * x[2])   # 3,1 
    values[8] += lambda[0] * (x[0] * x[1])   # 3,2 

    # add the portion for the second constraint
    values[0] += lambda[1] * 2.0     # 0,0 

    values[2] += lambda[1] * 2.0     # 1,1 

    values[5] += lambda[1] * 2.0     # 2,2 

    values[9] += lambda[1] * 2.0     # 3,3 

  return TRUE


#------------------------------------------------------------------
when isMainModule:
    #var nlp:IpoptProblem                 # IpoptProblem 
    #var status:ApplicationReturnStatus   # Solve return code 

    #
    # Define the problem
    #

    #var i:Index                          # generic counter 

    # Variable limits (bounds)
    let x_L = [1.0,1.0,1.0,1.0]   # lower bounds on x 
    let x_U = [5.0,5.0,5.0,5.0]   # upper bounds on x 
    assert(x_L.len == 4)
    assert(x_U.len == 4)


    # Constrains bounds
    # set the number of constraints and allocate space for the bounds 
    let g_L = [25.0,40.0]         # lower bounds on g 
    let g_U = [2e19,40.0]         # upper bounds on g 

    assert(g_L.len == 2)
    assert(g_U.len == 2)

    # set the number of nonzeros in the Jacobian and Hessian 
    let nele_jac = 8              # number of nonzeros in the Jacobian of the constraints 
    let nele_hess = 10            # number of nonzeros in the Hessian of the Lagrangian (lower or upper triangular part only) 

    # set the indexing style to C-style (start counting of rows and column indices at 0)
    let index_style = 0           # indexing style for matrices 

    # create the IpoptProblem
    let nlp:IpoptProblem = CreateIpoptProblem( 
          x_L.len.Index, 
          cast[ptr Number](unsafeAddr(x_L)), 
          cast[ptr Number](unsafeAddr(x_U)), 
          g_L.len.Index, 
          cast[ptr Number](unsafeAddr(g_L)), 
          cast[ptr Number](unsafeAddr(g_U)),
          nele_jac.Index, 
          nele_hess.Index, 
          index_style.Index,
          eval_f, 
          eval_g, 
          eval_grad_f,
          eval_jac_g, 
          eval_h )
    #nlp = ipopt.problem(
    #            n=len(x0),
    #            m=len(cl),
    #            problem_obj=hs071(),
    #            lb=lb,
    #            ub=ub,
    #            cl=cl,
    #            cu=cu
    #            )                         

    # NOTE: the values are copied internally so we could free them

    if nlp == nil:
      raise newException(ValueError, "CreateIpoptProblem failed - there was a problem with one of the inputs")

    # Set some options.  Note the following ones are only examples,
    # they might not be suitable for your problem.

    var ret = AddIpoptNumOption(nlp, "tol".cstring, 1e-7.Number)
    assert( ret == TRUE, "NumOption" )
    ret = AddIpoptStrOption(nlp, "mu_strategy", "adaptive")
    assert( ret == TRUE, "StrOption" )
    ret = AddIpoptStrOption(nlp, "output_file", "ipopt.out")
    assert( ret == TRUE, "setting option `output_file`failed" )

    # SOLVE THE PROBLEM
    # -----------------

    # starting point and solution vector (allocate space for the initial point and set the values)
    var x = [1.0, 5.0, 5.0, 1.0]

    # allocate space to store the bound multipliers at the solution 
    var
      mult_g:array[2, Number]         # constraint multipliers at the solution 
      mult_x_L:array[4, Number]       # lower bound multipliers at the solution
      mult_x_U:array[4, Number]       # upper bound multipliers at the solution 
      obj:Number                      # objective value

    # Initialize the user data 
    var user_data = MyUserData(g_offset:[0.0,0.0])             # our user data for the function evaluations 

    # Set the callback method for intermediate user-control.
    # This is not required, just gives you some intermediate control in
    # case you need it.
    # /* SetIntermediateCallback(nlp, intermediate_cb); */

    # solve the problem   
    #[
    let status = IpoptSolve( nlp, # Problem that is to be optimized
                    cast[ptr Number](unsafeAddr(x[0])),  # Input:  Starting point
                    nil, # Values of constraint at final point
                    cast[ptr Number](unsafeAddr(obj)), # Final value of objective function
                    cast[ptr Number](unsafeAddr(mult_g[0])),   # Input: Initial values for the constraint
                    cast[ptr Number](unsafeAddr(mult_x_L[0])), # Input: Initial values for the multipliers
                    cast[ptr Number](unsafeAddr(mult_x_U[0])), # Input: Initial values for the multipliers
                    cast[ptr MyUserData](unsafeAddr(user_data)) )
    ]#
    #let px = cast[ptr Number](x.addr)
    let status = IpoptSolve( nlp, # Problem that is to be optimized
                    cast[ptr Number](x.addr),  # Input:  Starting point
                    cast[ptr Number](nil), # Values of constraint at final point
                    cast[ptr Number](obj.addr), # Final value of objective function
                    cast[ptr Number](mult_g.addr),   # Input: Initial values for the constraint
                    cast[ptr Number](mult_x_L.addr), # Input: Initial values for the multipliers
                    cast[ptr Number](addr mult_x_U), # Input: Initial values for the multipliers
                    cast[pointer](addr user_data) )

    if status == Solve_Succeeded:
      echo "\n\nSolution of the primal variables, x\n"
      for i in 0..<x_L.len:
        echo &"x[{i}] = {x[i]}\n"
      
      echo "\n\nSolution of the constraint multipliers, lambda\n"
      for i in 0..g_L.len:
        echo &"lambda[{i}] = {mult_g[i]}\n"
      
      echo "\n\nSolution of the bound multipliers, z_L and z_U\n"
      for i in 0..<x_L.len:
        echo &"z_L[{i}] = {mult_x_L[i]}\n"
      
      for i in 0..<x_L.len:
        echo &"z_U[{i}] = {mult_x_U[i]}\n"
    
      echo "\n\nObjective value\nf(x*) = {obj}\n"
   
    else:  
      echo "\n\nERROR OCCURRED DURING IPOPT OPTIMIZATION.\n"


#[
    
    def objective(self, x):
        #
        # The callback for calculating the objective
        #
        return x[0] * x[3] * np.sum(x[0:3]) + x[2]
        
    def gradient(self, x):
        #
        # The callback for calculating the gradient
        #
        return np.array([
                    x[0] * x[3] + x[3] * np.sum(x[0:3]), 
                    x[0] * x[3],
                    x[0] * x[3] + 1.0,
                    x[0] * np.sum(x[0:3])
                    ])
    
    def constraints(self, x):
        #
        # The callback for calculating the constraints
        #
        return np.array((np.prod(x), np.dot(x, x)))
    
    def jacobian(self, x):
        #
        # The callback for calculating the Jacobian
        #
        return np.concatenate((np.prod(x) / x, 2*x))
    
    def hessianstructure(self):
        #
        # The structure of the Hessian
        # Note:
        # The default hessian structure is of a lower triangular matrix. Therefore
        # this function is redundant. I include it as an example for structure
        # callback.
        # 
        global hs
        
        hs = sps.coo_matrix(np.tril(np.ones((4, 4))))
        return (hs.col, hs.row)
    
    def hessian(self, x, lagrange, obj_factor):
        #
        # The callback for calculating the Hessian
        #
        H = obj_factor*np.array((
                (2*x[3], 0, 0, 0),
                (x[3],   0, 0, 0),
                (x[3],   0, 0, 0),
                (2*x[0]+x[1]+x[2], x[0], x[0], 0)))
                
        H += lagrange[0]*np.array((
                (0, 0, 0, 0),
                (x[2]*x[3], 0, 0, 0),
                (x[1]*x[3], x[0]*x[3], 0, 0),
                (x[1]*x[2], x[0]*x[2], x[0]*x[1], 0)))
                
        H += lagrange[1]*2*np.eye(4)
    
        #
        # Note:
        # 
        #
        return H[hs.row, hs.col]

    def intermediate(
            self, 
            alg_mod,
            iter_count,
            obj_value,
            inf_pr,
            inf_du,
            mu,
            d_norm,
            regularization_size,
            alpha_du,
            alpha_pr,
            ls_trials
            ):

        #
        # Example for the use of the intermediate callback.
        #
        print("Objective value at iteration #%d is - %g" % (iter_count, obj_value))
    
    
def main():


    #
    # Set solver options
    #
    #nlp.addOption('derivative_test', 'second-order')
    nlp.addOption(b'mu_strategy', b'adaptive')
    nlp.addOption(b'tol', 1e-7)

    #
    # Scale the problem (Just for demonstration purposes)
    #
    nlp.setProblemScaling(
        obj_scaling=2,
        x_scaling=[1, 1, 1, 1]
        )
    nlp.addOption(b'nlp_scaling_method', b'user-scaling')
    
    #
    # Solve the problem
    #
    x, info = nlp.solve(x0)
    
    print("Solution of the primal variables: x=%s\n" % repr(x))
    
    print("Solution of the dual variables: lambda=%s\n" % repr(info['mult_g']))
    
    print("Objective=%s\n" % repr(info['obj_val']))


if __name__ == '__main__':
    main()    
]#