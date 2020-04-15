include "IpStdCInterface"
import strformat
import macros
import options
#TODO: en solve que reciba como parámetro el número de constrains

#import sugar
#[
type
  Problem = object
    x_L, x_U, g_L, g_U:  openArray[Number]
    nele_jac, nele_hess: int
    eval_f:      Eval_F_CB
    eval_g:      Eval_G_CB
    eval_grad_f: Eval_Grad_F_CB
    eval_jac_g:  Eval_Jac_G_CB 
    eval_h:      Eval_H_CB
    nlp: IpoptProblem
]#

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


proc solve*(problem: IpoptProblem, ini:seq[Number]):tuple[obj:Number, sol:seq[Number], mult_g:seq[Number], mult_x_L:seq[Number],  mult_x_U:seq[Number]] =
  ##[
  - `problem`: constains the problem to be solved
  - `ini`: starting point and solution vector (allocate space for the initial point and set the values)
  ]##
  var
    x = ini
    mult_g = newSeq[Number](2)#array[2, Number]         # constraint multipliers at the solution 
    mult_x_L = newSeq[Number](ini.len) #:array[4, Number]       # lower bound multipliers at the solution
    mult_x_U = newSeq[Number](ini.len)       # upper bound multipliers at the solution 
    obj:Number   

  let status = IpoptSolve( problem, # Problem that is to be optimized
                 cast[ptr Number](x[0].addr),  # Input:  Starting point
                 cast[ptr Number](nil), # Values of constraint at final point
                 cast[ptr Number](obj.unsafeAddr), # Final value of objective function
                 cast[ptr Number](mult_g[0].addr),   # Input: Initial values for the constraint
                 cast[ptr Number](mult_x_L[0].addr), # Input: Initial values for the multipliers
                 cast[ptr Number](mult_x_U[0].addr), # Input: Initial values for the multipliers
                 nil)
                #cast[ptr MyUserData](addr user_data) )
  doAssert( status == Solve_Succeeded, "ERROR OCCURRED DURING IPOPT OPTIMIZATION")
  return (obj, x, mult_g, mult_x_L, mult_x_U)

proc solveWarm*(problem: IpoptProblem, ini, mult_gg, mult_x_LL, mult_x_UU:seq[Number]):tuple[obj:Number, sol:seq[Number], mult_g:seq[Number], mult_x_L:seq[Number],  mult_x_U:seq[Number]] =
  var
    x = ini
    mult_g = mult_gg     #array[2, Number]         # constraint multipliers at the solution 
    mult_x_L = mult_x_LL #:array[4, Number]       # lower bound multipliers at the solution
    mult_x_U = mult_x_UU # upper bound multipliers at the solution 
    obj:Number 
  problem["warm_start_init_point"] = "yes"
  let status = IpoptSolve( problem, # Problem that is to be optimized
                 cast[ptr Number](x[0].addr),  # Input:  Starting point
                 cast[ptr Number](nil), # Values of constraint at final point
                 cast[ptr Number](obj.unsafeAddr), # Final value of objective function
                 cast[ptr Number](mult_g[0].addr),   # Input: Initial values for the constraint
                 cast[ptr Number](mult_x_L[0].addr), # Input: Initial values for the multipliers
                 cast[ptr Number](mult_x_U[0].addr), # Input: Initial values for the multipliers
                 nil) 

  doAssert( status == Solve_Succeeded, "ERROR OCCURRED DURING IPOPT OPTIMIZATION")
  return (obj, x, mult_g, mult_x_L, mult_x_U)



#------------------------
# Macro
#------------------------
#echo item.astGenRepr

proc parseLimits(body:NimNode):(NimNode,int) =
  # Parses the upper/lower limits for the variables. Returns the number of elements.
  var counter:int = 0
  assert( body[1][0][0].strVal == "@", "a sequence shall be provided" )
  let tmp = body[1][0]
  return (tmp, tmp[1].len)



#[
proc procToLambda(body:NimNode):NimNode = 
  var data = newTree(nnkLambda)
  for i in body:
    data.add i
  data[0] = newEmptyNode() #.astGenRepr
  #var data2 = newTree(nnkLambda)
  return data

proc parseObjective(body:NimNode, fname:NimNode):NimNode =
  let tmp = body[1][0]
  let x = newIdentNode("x")
  var data = quote do:  ##`fname`
      proc `fname`( n:Index, xx:ptr Number, new_x:Bool, obj_value:ptr Number, user_data:UserDataPtr ):Bool  {.cdecl,exportc.} =
        ## Callback function for evaluating objective function
        #assert(n == 4)
        let `x` = cast[ptr UncheckedArray[Number]](xx) 

        obj_value[] = `tmp`
        return TRUE
  return procToLambda(data)
  #echo result.astGenRepr
]#
proc parseObjective(body:NimNode, fname:NimNode):NimNode =
  let tmp = body[1][0]
  let x = newIdentNode("x")
  return quote do:  ##`fname`
      proc `fname`( n:Index, xx:ptr Number, new_x:Bool, obj_value:ptr Number, user_data:UserDataPtr ):Bool  {.cdecl,exportc.} =
        ## Callback function for evaluating objective function
        #assert(n == 4)
        let `x` = cast[ptr UncheckedArray[Number]](xx) 

        obj_value[] = `tmp`
        return TRUE  

proc parseGrad(body:NimNode, fname:NimNode):NimNode =
  let tmp = body[1][0]
  let x = newIdentNode("x")
  #let n = x_L.len
  return quote do:
    proc `fname`(n: Index; 
              xx: ptr Number; # Input
              new_x: Bool; 
              grad_ff: ptr Number; # Output
              user_data: UserDataPtr): Bool {.cdecl, exportc.} =
        ## Callback function for evaluating objective function
        #assert(n == 4)
        let `x` = cast[ptr UncheckedArray[Number]](xx) 
        let grad_f = cast[ptr UncheckedArray[Number]](grad_ff)
        for i in 0..<`tmp`.len:
          grad_f[i] = `tmp`[i]
        #let tmp = `tmp` #buf[0]*buf[3]*(buf[0]+buf[1]+buf[2]) + buf[2]
        #obj_value[] = `tmp`
        return TRUE  

proc parseObjHess(body:NimNode ):NimNode =
  assert( body[1][0][0].strVal == "@", "a sequence shall be provided" )
  return body[1][0]


proc parseConstrain(body:NimNode ):tuple[lower, upper, function, grad, hess:NimNode] =
  let x = newIdentNode("x")
  var 
    lower, upper, function:NimNode
    grad:NimNode
    hess:NimNode

  for item in body[1]:
    #echo repr item[0]
    if eqIdent(item[0], "lower_limit"):
      lower = item[1]
      
    if eqIdent(item[0], "upper_limit"):
      upper = item[1]
    if eqIdent(item[0], "function"):
      function = item[1]    
    if eqIdent(item[0], "grad"):
      assert( item[1][0][0].strVal == "@", "a sequence shall be provided" )
      grad = item[1][0]
      #for i in item[1][0]:
      #  grad &= i

    if eqIdent(item[0], "hess"):
      assert( item[1][0][0].strVal == "@", "a sequence shall be provided" )
      hess = item[1][0]       

  return (lower, upper, function, grad, hess)        

proc genJacG(fname:NimNode, g_grad:seq[NimNode] ):tuple[data:NimNode, n_constrains:seq[int]] =
  # Jacobian
  let x = newIdentNode("x")
  var n_constrains:seq[int]
  for i in 0..<g_grad.len:
    n_constrains &= g_grad[i][1].len
  echo repr n_constrains
  let data = quote do:
    proc `fname`( n: Index; 
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
        let iRow = cast[ptr UncheckedArray[Index]](iiRow)
        let jCol = cast[ptr UncheckedArray[Index]](jjCol)
        var idx = 0
        for row in 0..<`n_constrains`.len:
          for col in 0..<`n_constrains`[row]:
            iRow[idx] = row.Index
            jCol[idx] = col.Index
            idx += 1
      else:
        # return the values of the jacobian of the constraints 
        let `x`    = cast[ptr UncheckedArray[Number]](xx)
        let values = cast[ptr UncheckedArray[Number]](vvalues)
        var idx = 0
        for row in 0..<`n_constrains`.len:   # For each constrain
          for col in 0..<`n_constrains`[row]:
            #let tmp = `g_grad`[row][col]
            values[idx] = `g_grad`[row][col]
            idx += 1

      return TRUE
  return (data, n_constrains)

proc genHessian(fname:NimNode, tmpObjHess:NimNode, tmpConsHess:seq[NimNode], nn:int, n_constrains:seq[int] ):NimNode =
  let x = newIdentNode("x")
  return quote do:
    proc `fname`( n: Index; 
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

      if vvalues == nil:
        # return the structure. This is a symmetric matrix, fill the lower left triangle only. 
        # the hessian for this problem is actually dense
        let iRow = cast[ptr UncheckedArray[Index]](iiRow)
        let jCol = cast[ptr UncheckedArray[Index]](jjCol)
        var idx = 0
        for row in 0..<`nn`:
          for col in 0..row:
            iRow[idx] = row.Index
            jCol[idx] = col.Index
            #assert(idx < nele_hess)
            idx += 1
            #echo idx
        #assert(idx == nele_hess) 

      else:
        # return the values. This is a symmetric matrix, fill the lower left triangle only
        
        # fill the objective portion
        let `x`    = cast[ptr UncheckedArray[Number]](xx)
        let values = cast[ptr UncheckedArray[Number]](vvalues)
        var idx = 0

        for row in 0..<`nn`:
          for col in 0..row:
            values[idx] = obj_factor * `tmpObjHess`[row][col]
            idx += 1
    
        # add the portion for the constraints
        var nConst = 0
        let lambda = cast[ptr UncheckedArray[Number]](llambda)
        
        for constrain in 0..<`n_constrains`.len:
          idx = 0
          for row in 0..<`nn`:
            for col in 0..row:
              #let tmp = `tmpObjHess`[constrain][row]
              values[idx] += lambda[constrain] * `tmpConsHess`[constrain][row][col]
              idx += 1

      return TRUE

proc genConstrains(fname:NimNode, g:seq[NimNode]):NimNode =
  let nG = len(g)
  let x = newIdentNode("x")
  return quote do:
    proc `fname`( n: Index; 
             xx: ptr Number; 
             new_x: Bool; 
             m: Index; 
             gg: ptr Number; 
             user_data: UserDataPtr): Bool  {.cdecl, exportc.} =
      ##Callback function for evaluating constraint functions
      #assert(n == 4)
      #assert(m == 2)
      
      #let my_data = cast[ptr MyUserData](user_data)
      let `x`       = cast[ptr UncheckedArray[Number]](xx)
      let gArray  = cast[ptr UncheckedArray[Number]](gg)

      for i in 0..<`nG`: 
        gArray[i] = `g`[i]   #x[0] * x[1] * x[2] * x[3] + my_data.g_offset[0]
      #g[1] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3] + my_data.g_offset[1]
      return TRUE


proc parseOptions(body, nlp:NimNode ):NimNode =
  result = nnkStmtList.newTree()
  #let nlp = newIdentNode("nlp")
  for item in body[1]:
    echo repr item[0], "====", repr item[1][0].astGenRepr 
    let key = item[0].strVal
    let value = item[1][0]
    result.add quote do:
      `nlp`[`key`] = `value`


proc genSequence(data:seq[NimNode]):NimNode =
  var body:NimNode = newTree(nnkPrefix, newIdentNode("@"))
  body.add newTree(nnkBracket) 

  for i in data:
    body[1].add i[0]
  return body

macro problem*(name:string, mbody:untyped):untyped =
  result = nnkStmtList.newTree()
  var xLower,xUpper:NimNode
  #let x_L = newIdentNode("x_L")
  #let x_U = newIdentNode("x_U")
  let eval_f = genSym(nskProc, ident="eval_ff")   #newIdentNode()  #nskFunc  # genSym(nskLabel)#
  var objective:NimNode
  let eval_g = genSym(nskProc, ident="eval_gg")  
  let eval_grad_f = genSym(nskProc, ident="eval_grad_ff")   #newIdentNode()  #nskFunc  # genSym(nskLabel)#
  let eval_jac_g = genSym(nskProc, ident="eval_jac_gg")
  let eval_h = genSym(nskProc, ident="eval_hh")  
  let g_L = newIdentNode("g_L")  
  let g_U = newIdentNode("g_U")
  let nlp = newIdentNode(name.strVal)


  var x_Lvec:seq[NimNode]
  var x_Uvec:seq[NimNode]
  var g_lower:seq[NimNode]
  var g_upper:seq[NimNode]
  var g:seq[NimNode]
  var g_grad:seq[NimNode]
  var g_hess:seq[ seq[seq[NimNode]] ]
  #var g_constrains:
  var nLower = 0
  var nUpper = 0
  #var m = 0
  var tmpObjHess:NimNode
  var tmpConsHess:seq[NimNode]
  var options:NimNode
  var customData:NimNode

  for body in mbody:
    body.expectKind nnkCall
    var tmp:NimNode
    if eqIdent(body[0], "upper_limits"):
      (xUpper, nUpper) = parseLimits(body)

    if eqIdent(body[0], "lower_limits"):
      (xLower, nLower) = parseLimits(body)          

    if eqIdent(body[0], "objective"):  
      result.add parseObjective(body, eval_f)
      #let obj2 = procToLambda(objective)
      #echo obj2.astGenRepr

    if eqIdent(body[0], "grad"): 
      result.add parseGrad(body, eval_grad_f)

    if eqIdent(body[0], "hess"):      
      #g_hess &= parseHess(body)
      tmpObjhess = parseObjHess(body)

    if eqIdent(body[0], "constrain"):        
      #let tmp = body[1][0]
      let (lower, upper, function, grad, hess) = parseConstrain(body)
      g_lower &= lower
      g_upper &= upper
      g &= function
      g_grad &= grad
      tmpConsHess &= hess

    if eqIdent(body[0], "options"):      
      options = parseOptions(body, nlp)


  let (data, n_constrains) = genJacG(eval_jac_g, g_grad)
  result.add data

  # Hessian
  let hessian = genHessian(eval_h, tmpObjHess, tmpConsHess, nLower, n_constrains)
  result.add hessian

  # ----- eval_g -------------------------
  result.add genConstrains(eval_g, g)

  # set the number of nonzeros in the Jacobian and Hessian 
  var nele_jac = 0       # number of nonzeros in the Jacobian of the constraints 
  for i in n_constrains:
    nele_jac += i


  let nele_hess = (nLower * (nLower + 1)/2).int     # number of nonzeros in the Hessian of the Lagrangian (lower or upper triangular part only) 

  var ggLower = genSequence( g_lower )
  var ggUpper = genSequence( g_upper )

  

  #result.add quote do:
  #  var `eval_f` = `objective`
    
  result.add quote do:
    let `nlp`:IpoptProblem = createProblem( 
              `xLower`, `xUpper`, `ggLower`, `ggUpper`,
              `nele_jac`, `nele_hess`,
              `eval_f`, `eval_g`, `eval_grad_f`, `eval_jac_g`, `eval_h` )

  result.add options

  #echo repr g_lower.astGenRepr
  #result.add quote do:
  #  nlp["tol"] = 1e-7

