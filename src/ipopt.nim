include "IpStdCInterface"
import strformat
import macros
import options
import sugar

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

#------------------------
# Macro
#------------------------
#echo item.astGenRepr

proc parseLimits(body:NimNode, id:NimNode):(NimNode,int) =
  # Parses the upper/lower limits for the variables. Returns the number of elements.
  var counter:int = 0
  let tmp = body[1][0]
  var data = quote do:
    let `id` = `tmp`
  
  for i in tmp:
    counter += 1
  return (data, counter)

proc parseObjective(body:NimNode, fname:NimNode):NimNode =
  let tmp = body[1][0]
  let x = newIdentNode("x")
  #let n = x_L.len
  return quote do:
      proc `fname`( n:Index,
          x:ptr Number, 
          new_x:Bool,
          obj_value:ptr Number, 
          user_data:UserDataPtr ):Bool  {.cdecl,exportc.} =
        ## Callback function for evaluating objective function
        #assert(n == 4)
        let `x` = cast[ptr UncheckedArray[Number]](x) 

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


proc parseConstrain(body:NimNode ):tuple[lower:NimNode, upper:NimNode, function:NimNode, grad:seq[NimNode], hess:NimNode] =
  let x = newIdentNode("x")
  var 
    lower, upper, function:NimNode
    grad:seq[NimNode]
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
      
      for i in item[1][0]:
        grad &= i

    if eqIdent(item[0], "hess"):
      hess = item[1][0]       

  return (lower, upper, function, grad, hess)        

proc genJacG(fname:NimNode, g_grad:seq[seq[NimNode]] ):tuple[data:NimNode, n_constrains:seq[int]] =
  # Jacobian
  let x = newIdentNode("x")
  var n_constrains:seq[int]
  #let n = nLower
  for i in 0..<g_grad.len:
    n_constrains &= g_grad[i].len
  #echo repr n_constrains
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


proc parseOptions(body:NimNode ):NimNode =
  #let x = newIdentNode("x")
  #var 
  #  lower, upper, function:NimNode
  #  grad:seq[NimNode]
  #  hess:NimNode
  result = nnkStmtList.newTree()
  let nlp = newIdentNode("nlp")
  for item in body[1]:
    echo repr item[0], "====", repr item[1][0].astGenRepr 
    let key = item[0].strVal
    let value = item[1][0]
    #echo key
    result.add quote do:
      `nlp`[`key`] = `value`
  echo repr result
    #   `[]=`( `nlp`, `key`, `item`[1] )


    # []=(nlp, "tol", 1e-07)
    #echo repr item[0]
    #if eqIdent(item[0], "tol"):
    #  lower = item[1]
      
    #if eqIdent(item[0], "mu_strategy"):
    #  upper = item[1]
    #if eqIdent(item[0], "output_file"):
    #  function = item[1] 


macro problem*(mbody:untyped):untyped =
  result = nnkStmtList.newTree()
  let x_L = newIdentNode("x_L")
  let x_U = newIdentNode("x_U")
  let eval_f = genSym(nskProc, ident="eval_ff")   #newIdentNode()  #nskFunc  # genSym(nskLabel)#
  let eval_g = genSym(nskProc, ident="eval_gg")  
  let eval_grad_f = genSym(nskProc, ident="eval_grad_ff")   #newIdentNode()  #nskFunc  # genSym(nskLabel)#
  let eval_jac_g = genSym(nskProc, ident="eval_jac_gg")
  let eval_h = genSym(nskProc, ident="eval_hh")  
  let g_L = newIdentNode("g_L")  
  let g_U = newIdentNode("g_U")
  let nlp = newIdentNode("nlp")


  var x_Lvec:seq[NimNode]
  var x_Uvec:seq[NimNode]
  var g_lower:seq[NimNode]
  var g_upper:seq[NimNode]
  var g:seq[NimNode]
  var g_grad:seq[seq[NimNode]]
  var g_hess:seq[ seq[seq[NimNode]] ]
  #var g_constrains:
  var nLower = 0
  var nUpper = 0
  #var m = 0
  var tmpObjHess:NimNode
  var tmpConsHess:seq[NimNode]
  var options:NimNode

  for body in mbody:
    body.expectKind nnkCall
    var tmp:NimNode
    if eqIdent(body[0], "upper_limits"):
      (tmp, nUpper) = parseLimits(body, x_U)
      result.add tmp

    if eqIdent(body[0], "lower_limits"):
      (tmp, nLower) = parseLimits(body, x_L)
      result.add tmp      

    if eqIdent(body[0], "objective"):  
      result.add parseObjective(body, eval_f)


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
      options = parseOptions(body)

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


  result.add quote do:
    let `g_L` = `g_lower` 
    let `g_U` = `g_upper`     
  
  result.add quote do:
    let `nlp`:IpoptProblem = createProblem( 
              `x_L`, `x_U`, `g_L`, `g_U`,  #`g_L`, `g_U`, 
              `nele_jac`, `nele_hess`,
              `eval_f`, `eval_g`, `eval_grad_f`, `eval_jac_g`, `eval_h` )

  result.add options
  #result.add quote do:
  #  nlp["tol"] = 1e-7

