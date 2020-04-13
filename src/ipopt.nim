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



#[
template objective(x:untyped, body:untyped):untyped =   #openArray[Number]
  proc eval_ff( n:Index,
                x:ptr Number, 
                new_x:Bool,
                obj_value:ptr Number, 
                user_data:UserDataPtr ):Bool  {.cdecl,exportc.} =
            ## Callback function for evaluating objective function
            assert(n == 4)
            let `x` = cast[ptr UncheckedArray[Number]](x) 
  
            #let tmp = `tmp` #buf[0]*buf[3]*(buf[0]+buf[1]+buf[2]) + buf[2]
            #obj_value[] = `body`
            return TRUE
expandMacros:
  objective(x):
    x[0]*x[3]*(x[0]+x[1]+x[2]) + x[2]
]#

#[

]#


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

#def eval_f(x, user_data = None):
#    assert len(x) == 4
#    return x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2]

#[
proc f(x:openArray[Number]):Number = x[0]*x[3]*(x[0]+x[1]+x[2]) + x[2]

#
proc gen_eval_f( vec:openArray[Number],
                 f1: proc(vec:openArray[Number]):Number ):Eval_F_CB = 
  return proc eval_f( n:Index,
               x:ptr Number, 
               new_x:Bool,
               obj_value:ptr Number, 
               user_data:UserDataPtr):Bool  {.cdecl,exportc.} =
    assert n == vec.len
    obj_value[] = f1( cast[ptr UncheckedArray[Number]](x) )
    return TRUE
]#

#eval_f = gen_eval_f()

#[
template objectiveFunction(body:untyped):untyped =
  proc eval_f( n:Index,
               x:ptr Number, 
               new_x:Bool,
               obj_value:ptr Number, 
               user_data:UserDataPtr):Bool  {.cdecl,exportc.} =
    obj_value[] = f( cast[ptr UncheckedArray[Number]](x) )
    return TRUE  


]#
#[
macro objective*( buf:untyped, body:untyped):untyped =
  #let vec = newIdentNode(x & "Vec") 

  result = nnkStmtList.newTree()
  result.add quote do:
    proc eval_f( n:Index,
                 x:ptr Number, 
                 new_x:Bool,
                 obj_value:ptr Number, 
                 user_data:UserDataPtr):Bool  {.cdecl,exportc.} =
      let buf = cast[ptr UncheckedArray[Number]](x)
      `body`      
      obj_value[] = f( cast[ptr UncheckedArray[Number]](x) )
      return TRUE       

objective(x):
  #assert( x.len == n)
  x[0]*x[3]*(x[0]+x[1]+x[2]) + x[2]  
]#
#[
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
]#
#proc gen_objective()

#[
macro hola(x:untyped, body:untyped):untyped =
  
  result = nnkStmtList.newTree()
  result.add quote do:
    proc eval_ff( n:Index,
                  x:ptr Number, 
                  new_x:Bool,
                  obj_value:ptr Number, 
                  user_data:UserDataPtr):Bool  {.cdecl,exportc.} =    
      `body`
      echo `x`


hola(x):
  let prueba = 3
]#
#echo 



#macro problem(body:untyped):untyped =
#  result = nnkStmtList.newTree()

#[
type
  UserData = object
  ObjectiveFunction = proc(x:openArray[Number];user_data:UserData):Number
  Problem = object
    upper_limits:seq[Number]
    lower_limits:seq[Number]
    objective:ObjectiveFunction

var p:Problem
proc setUpperLimits(p:var Problem, x:openArray[Number]) = 
  var s:seq[Number]
  for i in 0..<x.len:
    s &= x[i]
  p.upper_limits = s

proc setLowerLimits(p:var Problem, x:openArray[Number]) = 
  var s:seq[Number]
  for i in 0..<x.len:
    s &= x[i]  
  p.lower_limits = s

p.setUpperLimits( [1.0,1.0,1.0,1.0] )
p.setLowerLimits( [5.0,5.0,5.0,5.0] ) 
p.objective = proc (x:openArray[Number]) = x[0]*x[3]*(x[0]+x[1]+x[2]) + x[2]
]#
#[
macro upper_limits(x:untyped):untyped = #:untyped = 
  let x_U = newIdentNode("x_U")
  result = quote do:
    let `x_U` = `x`

macro lower_limits(x:untyped):untyped = #:untyped = 
  let x_L = newIdentNode("x_L")
  result = quote do:
    let `x_L` = `x`



echo ">>>>>> ", x_L
echo ">>>>>> ", x_U
]#

macro problem(mbody:untyped):untyped =
  result = nnkStmtList.newTree()
  let x_L = newIdentNode("x_L")
  let x_U = newIdentNode("x_U")
  let eval_f = genSym(nskProc, ident="eval_ff")   #newIdentNode()  #nskFunc  # genSym(nskLabel)#
  let eval_grad_f = genSym(nskProc, ident="eval_grad_ff")   #newIdentNode()  #nskFunc  # genSym(nskLabel)#
  let eval_jac_g = genSym(nskProc, ident="eval_jac_gg")
  let eval_h = genSym(nskProc, ident="eval_hh")  
  let g_L = newIdentNode("g_L")  
  let g_U = newIdentNode("g_U")

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
      for i in body[1]:   # Each part of the Hessian
        var tmp1:seq[seq[NimNode]] 
        #echo "==> ", repr i[1]
        for r in i[1]:
          var tmp2:seq[NimNode]
          for c in r[1]:
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
          for i in item[1][0]:
            tmp &= i
          
          g_grad &= tmp
        if eqIdent(item[0], "hess"):          
          for i in item[1]:   # Each part of the Hessian
            var tmp1:seq[seq[NimNode]] 
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
  let obj_hess = g_hess[0]
  for i in obj_hess:
    echo repr i
  #let n = x_U.len
  #echo ">> ", n
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
            let tmp = `obj_hess`[row][col] 
            #values[idx] = obj_factor * `obj_hess`[row][col]    #(2 * x[3])  # 0,0 
            #echo row, " ", col, "   ", repr `g_hess`[0][row][col]
            idx += 1
            #echo repr `obj_hess`[0][0]
            #echo row, " ", col

        #values[1] = obj_factor * (x[3])      # 1,0 
        #values[2] = 0                        # 1,1 

        #values[3] = obj_factor * (x[3])      # 2,0 
        #values[4] = 0                        # 2,1 
        #values[5] = 0                        # 2,2 

        #values[6] = obj_factor * (2 * x[0] + x[1] + x[2])   # 3,0 
        #values[7] = obj_factor * (x[0])                     # 3,1 
        #values[8] = obj_factor * (x[0])                     # 3,2 
        #values[9] = 0                                       # 3,3 

#[





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
]#














  #result.add quote do:
    #`body`
  #  assert( `x_L`.len == `x_U`.len, "Bounds limits should have the same length" )
  
  #echo repr result
#upper_limits [1.0,1.0,1.0,1.0] 
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
      function: x[0] * x[1] * x[2] * x[3] + data.g_offset[0]
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
      function: x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3] + data.g_offset[1]
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



#[
macro myMacro(body1,body2: untyped): untyped {.multiBodyMacro.} =
  echo body1.lispRepr # (StmtList (Command (Ident "echo") (IntLit 1)))
  echo body2.lispRepr # (StmtList (Command (Ident "echo") (IntLit 2)))

myMacro:
  body1:
    echo 1
  body2:
    echo 2
]#

#[
proc f(x:openArray[Number]):Number = x[0]*x[3]*(x[0]+x[1]+x[2]) + x[2]

proc gen_objective(f:proc(x:openArray[Number]):Number):Eval_F_CB =   #openArray[Number]
  return proc ( n:Index,
                x:ptr Number, 
                new_x:Bool,
                obj_value:ptr Number, 
                user_data:UserDataPtr ):Bool  {.cdecl,exportc.} =
            ## Callback function for evaluating objective function
            #assert(n == vec.len)
            let buf = cast[ptr UncheckedArray[Number]](x) 

            #let tmp = `tmp` #buf[0]*buf[3]*(buf[0]+buf[1]+buf[2]) + buf[2]
            obj_value[] = f(buf.toOpenArray(0,4))
            return TRUE
let eval_ff = gen_objective(f)
]#

#type
#  Limits = object
#[    
dumpAstGen:
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

nnkStmtList.newTree(
  nnkProcDef.newTree(
    newIdentNode("eval_f"),
    newEmptyNode(),
    newEmptyNode(),
    nnkFormalParams.newTree(
      newIdentNode("Bool"),
      nnkIdentDefs.newTree(
        newIdentNode("n"),
        newIdentNode("Index"),
        newEmptyNode()
      ),
      nnkIdentDefs.newTree(
        newIdentNode("x"),
        nnkPtrTy.newTree(
          newIdentNode("Number")
        ),
        newEmptyNode()
      ),
      nnkIdentDefs.newTree(
        newIdentNode("new_x"),
        newIdentNode("Bool"),
        newEmptyNode()
      ),
      nnkIdentDefs.newTree(
        newIdentNode("obj_value"),
        nnkPtrTy.newTree(
          newIdentNode("Number")
        ),
        newEmptyNode()
      ),
      nnkIdentDefs.newTree(
        newIdentNode("user_data"),
        newIdentNode("UserDataPtr"),
        newEmptyNode()
      )
    ),
    nnkPragma.newTree(
      newIdentNode("cdecl"),
      newIdentNode("exportc")
    ),
    newEmptyNode(),
    nnkStmtList.newTree(
      newCommentStmtNode("Callback function for evaluating objective function"),
      nnkCall.newTree(
        newIdentNode("assert"),
        nnkInfix.newTree(
          newIdentNode("=="),
          newIdentNode("n"),
          newLit(4)
        )
      ),
      nnkLetSection.newTree(
        nnkIdentDefs.newTree(
          newIdentNode("buf"),
          newEmptyNode(),
          nnkCast.newTree(
            nnkPtrTy.newTree(
              nnkBracketExpr.newTree(
                newIdentNode("UncheckedArray"),
                newIdentNode("Number")
              )
            ),
            newIdentNode("x")
          )
        )
      ),
      nnkLetSection.newTree(
        nnkIdentDefs.newTree(
          newIdentNode("tmp"),
          newEmptyNode(),
          nnkInfix.newTree(
            newIdentNode("+"),
            nnkInfix.newTree(
              newIdentNode("*"),
              nnkInfix.newTree(
                newIdentNode("*"),
                nnkBracketExpr.newTree(
                  newIdentNode("buf"),
                  newLit(0)
                ),
                nnkBracketExpr.newTree(
                  newIdentNode("buf"),
                  newLit(3)
                )
              ),
              nnkPar.newTree(
                nnkInfix.newTree(
                  newIdentNode("+"),
                  nnkInfix.newTree(
                    newIdentNode("+"),
                    nnkBracketExpr.newTree(
                      newIdentNode("buf"),
                      newLit(0)
                    ),
                    nnkBracketExpr.newTree(
                      newIdentNode("buf"),
                      newLit(1)
                    )
                  ),
                  nnkBracketExpr.newTree(
                    newIdentNode("buf"),
                    newLit(2)
                  )
                )
              )
            ),
            nnkBracketExpr.newTree(
              newIdentNode("buf"),
              newLit(2)
            )
          )
        )
      ),
      nnkAsgn.newTree(
        nnkBracketExpr.newTree(
          newIdentNode("obj_value")
        ),
        newIdentNode("tmp")
      ),
      nnkReturnStmt.newTree(
        newIdentNode("TRUE")
      )
    )
  )
)
]#


#[
dumpAstGen:
  let x_L = [1.0,2.0,3.0]

#[
nnkStmtList.newTree(
  nnkLetSection.newTree(
    nnkIdentDefs.newTree(
      newIdentNode("x_L"),
      newEmptyNode(),
      nnkBracket.newTree(
        newLit(1.0),
        newLit(2.0),
        newLit(3.0)
      )
    )
  )
)

]#
]#