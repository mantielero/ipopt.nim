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

macro problem(body:untyped):untyped =
  result = nnkStmtList.newTree()
  let x_L = newIdentNode("x_L")
  let x_U = newIdentNode("x_U")
  let eval_f = genSym(ident="eval_f")   #newIdentNode()  #nskFunc  # genSym(nskLabel)#

  for body in body:
    body.expectKind nnkCall
    if eqIdent(body[0], "upper_limits"):
      let tmp = body[1][0]
      result.add quote do:
        let `x_U` = `tmp`
    if eqIdent(body[0], "lower_limits"):
      let tmp = body[1][0]
      result.add quote do:
        let `x_L` = `tmp`        
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
            TRUE

  #echo repr body

  #result.add quote do:
    #`body`
  #  assert( `x_L`.len == `x_U`.len, "Bounds limits should have the same length" )
  
  #echo repr result
#upper_limits [1.0,1.0,1.0,1.0] 
#dumpAstGen:
expandMacros:
  problem:
    lower_limits: [1.0,1.0,1.0,1.0]
    upper_limits: [5.0,5.0,5.0,5.0]
    objective: x[0]*x[3]*(x[0]+x[1]+x[2]) + x[2]





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