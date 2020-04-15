import ipopt
import macros
import strformat
#dumpAstGen:

proc prettyPrint(obj:Number, sol, mult_g, mult_x_L, mult_x_U:seq[Number]) =
  echo "\n\nSolution of the primal variables, x\n"
  for i, xi in sol:
    echo &"x[{i}] = {xi}"
  
  echo "\n\nSolution of the constraint multipliers, lambda\n"
  for i,xi in mult_g:
    echo &"lambda[{i}] = {xi}"
  
  echo "\n\nSolution of the bound multipliers, z_L and z_U\n"
  for i, xi in mult_x_L:
    echo &"z_L[{i}] = {xi}"

  for i, xi in mult_x_U:
    echo &"z_U[{i}] = {xi}"

  echo &"\n\nObjective value\nf(x*) = {obj}\n"

var g_offset = @[0.0,0.0]

#expandMacros:
problem("nlp"):
    lower_limits: @[1.0,1.0,1.0,1.0]
    upper_limits: @[5.0,5.0,5.0,5.0]
    objective: x[0]*x[3]*(x[0]+x[1]+x[2]) + x[2]
    grad: 
      @[ x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]),
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
      function: x[0] * x[1] * x[2] * x[3] + g_offset[0]
      grad:
        @[ x[1] * x[2] * x[3],
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
      function: x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3] + g_offset[1]
      grad:
        @[ 2.0 * x[0], 
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


var ini:seq[Number] = @[1.0, 5.0, 5.0, 1.0]
var (obj, sol, mult_g, mult_x_L, mult_x_U) = nlp.solve(ini)
#echo repr obj
pretty_print(obj, sol, mult_g, mult_x_L, mult_x_U)

# Now a warm calculation
g_offset = @[0.2,0.0]
nlp["bound_push"] = 1e-5
nlp["bound_frac"] = 1e-5
(obj, sol, mult_g, mult_x_L, mult_x_U) = nlp.solveWarm( sol, mult_g, mult_x_L, mult_x_U)
pretty_print(obj, sol, mult_g, mult_x_L, mult_x_U)
