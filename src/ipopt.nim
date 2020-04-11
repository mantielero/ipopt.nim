include "IpStdCInterface"
import strformat

proc createProblem*( x_L, x_U, g_L, g_U:  openArray[Number],
                     nele_jac, nele_hess: int,
                     eval_f:      Eval_F_CB,
                     eval_g:      Eval_G_CB,
                     eval_grad_f: Eval_Grad_F_CB,
                     eval_jac_g:  Eval_Jac_G_CB, 
                     eval_h:      Eval_H_CB ):IpoptProblem =
  ##[
  - `eval_f`: Callback function for evaluating objective function
  - `eval_g`: Callback function for evaluating constraint functions
  - `eval_grad_f`: Callback function for evaluating gradient of objective function
  - `eval_jac_g`: Callback function for evaluating Jacobian of constraint functions
  - `eval_h`: 
  ]##
  doAssert( x_L.len == x_U.len, "Both `x_L` and `x_U` shall have the same length" )
  doAssert( g_L.len == g_U.len, "Both `g_L` and `g_U` shall have the same length" )

  # set the indexing style to C-style (start counting of rows and column indices at 0)
  let index_style = 0.Index  # indexing style for matrices (starts in 0)
  
  result = CreateIpoptProblem( 
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
  doAssert(result != nil, "CreateIpoptProblem failed - there was a problem with one of the inputs")

proc `[]=`*(ipopt_problem: IpoptProblem, key:string, val:string) =
  let ret = AddIpoptStrOption(ipopt_problem, key.cstring, val.cstring)
  doAssert( ret == TRUE, &"failed setting \"{key}\" with value: \"{val}\"" )

proc `[]=`*(ipopt_problem: IpoptProblem, key:string, val:Number) =
  let ret = AddIpoptNumOption(ipopt_problem, key.cstring, val )
  doAssert( ret == TRUE, &"failed setting \"{key}\" with value: {val}" )  



#[


proc solve(problem:IpoptProblem, init:openArray[Number]) =
    var
      #x = [1.0, 5.0, 5.0, 1.0]        # starting point and solution vector (allocate space for the initial point and set the values)
      mult_g:array[2, Number]         # constraint multipliers at the solution 
      mult_x_L:array[4, Number]       # lower bound multipliers at the solution
      mult_x_U:array[4, Number]       # upper bound multipliers at the solution 
      obj:Number                      # objective value

  let status = IpoptSolve( nlp, # Problem that is to be optimized
                  cast[ptr Number](x.addr),  # Input:  Starting point
                  cast[ptr Number](nil), # Values of constraint at final point
                  cast[ptr Number](obj.addr), # Final value of objective function
                  cast[ptr Number](mult_g.addr),   # Input: Initial values for the constraint
                  cast[ptr Number](mult_x_L.addr), # Input: Initial values for the multipliers
                  cast[ptr Number](addr mult_x_U), # Input: Initial values for the multipliers
                  cast[ptr MyUserData](addr user_data) )  

]#

#template objective_function()