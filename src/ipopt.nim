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

  for body in mbody:
    body.expectKind nnkCall
    if eqIdent(body[0], "upper_limits"):
      let tmp = body[1][0]
      #echo repr "-----", tmp
      result.add quote do:
        let `x_U` = `tmp`
      
      for i in tmp:
        nUpper += 1

    if eqIdent(body[0], "lower_limits"):
      let tmp = body[1][0]
      result.add quote do:
        let `x_L` = `tmp` 

      for i in tmp:
        nLower += 1

    if eqIdent(body[0], "objective"):        
      let tmp = body[1][0]
      let x = newIdentNode("x")
      #let n = x_L.len
      result.add quote do:
          proc `eval_f`( n:Index,
              x:ptr Number, 
              new_x:Bool,
              obj_value:ptr Number, 
              user_data:UserDataPtr ):Bool  {.cdecl,exportc.} =
            ## Callback function for evaluating objective function
            #assert(n == 4)
            let `x` = cast[ptr UncheckedArray[Number]](x) 
  
            #let tmp = `tmp` #buf[0]*buf[3]*(buf[0]+buf[1]+buf[2]) + buf[2]
            obj_value[] = `tmp`
            return TRUE

    if eqIdent(body[0], "grad"):        
      let tmp = body[1][0]
      let x = newIdentNode("x")
      #let n = x_L.len
      result.add quote do:
        proc `eval_grad_f`(n: Index; 
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

    if eqIdent(body[0], "hess"):          

      for i in body[1]:   # Each part of the Hessian  (nnkStmtList)
        #if i[0].strVal != "@": #Ident (@)       
        #  raise newException()
        assert( i[0].strVal == "@", "a sequence shall be provided" )
        #echo i.astGenRepr  # lispSt 
        var tmp1:seq[seq[NimNode]] 
        #echo "==> ", repr i.type
        tmpObjHess = i
        for r in i[1]:
          var tmp2:seq[NimNode]
          for c in r[1]:
            #echo "---> ", repr c, " ", type(c)
            #echo c.astGenRepr
            tmp2 &= c
          
          tmp1 &= tmp2
        g_hess &= tmp1
      #echo ">> ", repr g_hess

    if eqIdent(body[0], "constrain"):        
      #let tmp = body[1][0]
      let x = newIdentNode("x")
      for item in body[1]:
        #echo repr item[0]
        if eqIdent(item[0], "lower_limit"):
          g_lower &= item[1]
         
        if eqIdent(item[0], "upper_limit"):
          g_upper &= item[1]
        if eqIdent(item[0], "function"):
          g &= item[1]    
        if eqIdent(item[0], "grad"):
          var tmp:seq[NimNode]
          #echo item.astGenRepr
          for i in item[1][0]:
            tmp &= i
          
          g_grad &= tmp
        if eqIdent(item[0], "hess"):          
          for i in item[1]:   # Each part of the Hessian
            var tmp1:seq[seq[NimNode]]
            #echo ">>>> ", repr i 
            tmpConsHess &= i
            for r in i[1]:
              var tmp2:seq[NimNode]
              for c in r[1]:  
                #echo repr c 
                tmp2 &= c
              
              tmp1 &= tmp2
            g_hess &= tmp1
          
        #echo "==>", repr g_hess
          #echo repr g_hess
      #let n = x_L.len
      #[
      result.add quote do:
        proc `eval_grad_f`(n: Index; 
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
      ]#

  # Jacobian
  let x = newIdentNode("x")
  var n_constrains:seq[int]
  let n = nLower
  for i in 0..<g_grad.len:
    n_constrains &= g_grad[i].len
  #echo repr n_constrains
  result.add quote do:
    proc `eval_jac_g`( n: Index; 
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

  # Hessian
  #echo "----------"
  #for r in 0..<g_hess[0].len:
  #  for c in 0..<g_hess[0][r].len:
  #    echo ">> ", repr g_hess[0][r][c]

  #echo g_grad.astGenRepr
  let obj_hess = g_hess[0]
  for i in obj_hess:
    echo repr i," ", type(i)
  #echo type(obj_hess)
  echo "============================"
  #echo repr tmpConsHess
  #echo obj_hess.astGenRepr
  #let tmp = obj_hess[0][0]
  #let borrame = quote do:
  #   let a = @[@[1.0],@[1.0,2.0]]
  #echo borrame.astGenRepr
  #let n = x_U.len
  #echo ">> ", n
  for constrain in 0..<n_constrains.len:
    echo repr tmpConsHess[constrain]


  result.add quote do:
    proc `eval_h`( n: Index; 
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
      #var idx:Index = 0 # nonzero element counter
      #var row:Index = 0 # row counter for loop 
      #var col:Index = 0 # col counter for loop

      if vvalues == nil:
        # return the structure. This is a symmetric matrix, fill the lower left triangle only. 
        # the hessian for this problem is actually dense
        let iRow = cast[ptr UncheckedArray[Index]](iiRow)
        let jCol = cast[ptr UncheckedArray[Index]](jjCol)
        var idx = 0
        for row in 0..<`n`:
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

        for row in 0..<`n`:
          for col in 0..row:
            values[idx] = obj_factor * `tmpObjHess`[row][col]
            idx += 1
    
        # add the portion for the constraints
        var nConst = 0
        let lambda = cast[ptr UncheckedArray[Number]](llambda)
        
        for constrain in 0..<`n_constrains`.len:
          idx = 0
          for row in 0..<`n`:
            for col in 0..row:
              #let tmp = `tmpObjHess`[constrain][row]
              values[idx] += lambda[constrain] * `tmpConsHess`[constrain][row][col]
              idx += 1

      return TRUE

  # ----- eval_g -------------------------
  let nG = len(g)
  
  result.add quote do:
    proc `eval_g`( n: Index; 
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

  # set the number of nonzeros in the Jacobian and Hessian 
  let nele_jac = 8     # number of nonzeros in the Jacobian of the constraints 
  let nele_hess = 10     # number of nonzeros in the Hessian of the Lagrangian (lower or upper triangular part only) 

  result.add quote do:
    let `g_L` = `g_lower` 
    let `g_U` = `g_upper`     
  
  result.add quote do:
    let `nlp`:IpoptProblem = createProblem( 
              `x_L`, `x_U`, `g_L`, `g_U`,  #`g_L`, `g_U`, 
              `nele_jac`, `nele_hess`,
              `eval_f`, `eval_g`, `eval_grad_f`, `eval_jac_g`, `eval_h` )


