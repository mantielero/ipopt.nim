import ipopt
import macros

#dumpAstGen:

expandMacros:
  problem:
    data:
      g_offset: [0.0,0.0]
    lower_limits: [1.0,1.0,1.0,1.0]
    upper_limits: [5.0,5.0,5.0,5.0]
    objective: x[0]*x[3]*(x[0]+x[1]+x[2]) + x[2]
    grad: 
      [ x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]),
        x[0] * x[3],
        x[0] * x[3] + 1.0,
        x[0] * (x[0] + x[1] + x[2])  
      ]
    hess:
      @[ @[ 2.0 * x[3]],
         @[ x[3],                     0.0 ],
         @[ x[3],                     0.0,  0.0 ],
         @[ 2.0 * x[0] + x[1] + x[2], x[0], x[0], 0.0]
      ]
    constrain:
      lower_limit: 25.0
      upper_limit: high(Number)
      function: x[0] * x[1] * x[2] * x[3] #+ data.g_offset[0]
      grad:
        [ x[1] * x[2] * x[3],
          x[0] * x[2] * x[3],
          x[0] * x[1] * x[3],
          x[0] * x[1] * x[2]
        ]
      hess:
        @[ @[0.0],
           @[x[2] * x[3], 0.0],
           @[x[1] * x[3], x[0] * x[3], 0.0 ],
           @[x[1] * x[2], x[0] * x[2], x[0] * x[1], 0.0 ]           
        ]

    constrain:
      lower_limit: 40.0
      upper_limit: 40.0
      function: x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3] #+ data.g_offset[1]
      grad:
        [ 2.0 * x[0], 
          2.0 * x[1], 
          2.0 * x[2],
          2.0 * x[3]
        ]
      hess:
        @[ @[2.0],
           @[0.0, 2.0],
           @[0.0, 0.0, 2.0],
           @[0.0, 0.0, 0.0, 2.0]
        ]
    options:
      tol: 1e-7
      mu_strategy: "adaptive"
      output_file: "ipopt.out"


#nlp["tol"] = 1e-7
#nlp["mu_strategy"] = "adaptive"
#nlp["output_file"] = "ipopt.out"

var
    x = [1.0, 5.0, 5.0, 1.0]        # starting point and solution vector (allocate space for the initial point and set the values)
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
                nil)
                #cast[ptr MyUserData](addr user_data) )

echo "Objective: ", obj