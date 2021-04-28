
(******************************************************************************

    Building spinchain Hamiltonians : transverse-field anisotropic XY model

*******************************************************************************)


(* There are hopping terms and pairing terms - here are helper functions to build these. *)
hop = {
    Piecewise[{{Sigmap, j==1+i}, {Sigmam, j==i}, {Id, True}}],
    Piecewise[{{Sigmam, j==1+i}, {Sigmap, j==i}, {Id, True}}]    
};

pair = {
    Piecewise[{{Sigmap, j==1+i}, {Sigmap, j==i}, {Id, True}}],
    Piecewise[{{Sigmam, j==1+i}, {Sigmam, j==i}, {Id, True}}]    
};

(* Total spin operator for chain, length-n *)
Sz = Piecewise[{{Sigmaz, j==i}, {Id, True}}];
SzT[n_] := Total[Tensor @@@ Table[Sz, {i, 1, n}, {j, 1, n}]]


(* Now, the Hamiltonian *)

(* A utility function, toChain[] *)
(* Returns function which when given a hop or pair, etc returns 
   the corresponding sum: Sum[hopping terms] or Sum[Pairing terms] in the
   full composite-system Hilbert space. *)

toChain[n_] := Total[Tensor @@@ Flatten[Table[#, {i, 1, n - 1}, {j, 1, n}] & /@ #, 1]] &;

(* And then, the Hamiltonian: *)
H0[n_, g_] := toChain[n][hop] + g SzT[n];
V[n_, delta_] := delta toChain[n][pair];
HAnIsoXY[n_, g_, delta_] := H0[n, g] + V[n, delta];


(* ::Input::Initialization:: *)
PBC[n_,\[CapitalDelta]_]:=Block[
{
SigmapNSigmam1=Tensor@@Table[\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{"Sigmap", 
RowBox[{"j", "==", "n"}]},
{"Sigmam", 
RowBox[{"j", "==", "1"}]},
{"Id", "True"}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "ColumnsIndexed" -> {}, "Rows" -> {{1.}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "ColumnsIndexed" -> {}, "Rows" -> {{1.}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\),{j,n}],
Sigmap1SigmamN=Tensor@@Table[\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{"Sigmap", 
RowBox[{"j", "==", "1"}]},
{"Sigmam", 
RowBox[{"j", "==", "n"}]},
{"Id", "True"}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "ColumnsIndexed" -> {}, "Rows" -> {{1.}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "ColumnsIndexed" -> {}, "Rows" -> {{1.}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\),{j,n}],
SigmapNSigmap1=Tensor@@Table[\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{"Sigmap", 
RowBox[{"j", "==", "n"}]},
{"Sigmap", 
RowBox[{"j", "==", "1"}]},
{"Id", "True"}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "ColumnsIndexed" -> {}, "Rows" -> {{1.}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "ColumnsIndexed" -> {}, "Rows" -> {{1.}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\),{j,n}],
Sigmap1SigmapN=Tensor@@Table[\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{"Sigmap", 
RowBox[{"j", "==", "1"}]},
{"Sigmap", 
RowBox[{"j", "==", "n"}]},
{"Id", "True"}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "ColumnsIndexed" -> {}, "Rows" -> {{1.}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "ColumnsIndexed" -> {}, "Rows" -> {{1.}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\),{j,n}]
},
SigmapNSigmam1+Sigmap1SigmamN+\[CapitalDelta](SigmapNSigmap1+Sigmap1SigmapN)
]


(*
            Hamiltonian in the TOTAL SPIN BASIS.
*)

(* This function can act as U in U.p.Udag basis transformation *)
TotalSpinBasis[n_] := 
  Table[KroneckerDelta[i, j], {i, 2^n}, {j, 2^n}][[Ordering[SzT[n] // Normal // Diagonal]]] // SparseArray;
(* Produce the set of sub-spaces of the ISOTROPIC Hamiltonian *)
(* Helper: *)
f[m_, 0] := 1;
f[m_, r_] := f[m, r - 1] + Binomial[m, r - 1];

H0Blocks[n_, g_] := Module[
    {
        U = TotalSpinBasis[n],
        H0full = H0[n, g]
    },
    H0full = U . H0full . Transpose[U];
    Table[H0full[[f[n, m];;f[n, m + 1] - 1, f[n, m];;f[n, m + 1] - 1]], {m, 0, n}]
];

(*
            Spin lowering operators in the total spin basis
*)

(* Sigmam, length-m, at site i: *)
MinusBlocks[n_, i_] := Module[
    {
        U = TotalSpinBasis[n],
        SmT = Tensor@@Table[Piecewise[{{Sigmam, i==j}, {Id, True}}], {j, n}]
    },
    SmT = U . SmT . Transpose[U];
    Table[SmT[[f[n, m];;f[n, m + 1] - 1, f[n, m + 1];;f[n, m + 2] - 1]], {m, 0, n - 1}]
]

(* 
    Similarly, we can put the perturbative term into total spin basis, and extract
    only the terms which lower a state's spin number by two. 
    PlusBlocks and PlusPlusBlocks are obtained by mapping Transpose
*)
MinusMinusBlocks[n_, del_] := Module[
    {
        U = TotalSpinBasis[n],
        Vfull = V[n, del]
    },
    Vfull = U . Vfull . Transpose[U];
    Table[Vfull[[ f[n, m];;f[n, m+1]-1, f[n, m+2];;f[n, m+3]-1 ]], {m, 0, n-2}]
];


(* 
    Finally, of course we need to be able to invert this subspace projection, so 
    this builds us the full 2^N x 2^N matrix from the blocks. 
*)

Blockto2N[blockList_, mm_] := Module[
    {
        n = Length[blockList] - 1 + Abs[mm],
        out
    },
    out = SparseArray[{}, {2^n, 2^n}];
    Do[
        If[mm >= 0,
            out[[ f[n, i+mm];;f[n, i+mm+1]-1, f[n, i];;f[n, i+1]-1 ]] = blockList[[i + 1]],
            out[[ f[n, i];;f[n, i+1]-1, f[n, i-mm];;f[n, i-mm+1]-1 ]] = blockList[[i + 1]]
        ],
        
    {i, 0, Length[blockList] - 1}
    ];
    out
]




(* ::Input::Initialization:: *)
OnesBlock[n_,r_]:=Table[1,{i,Binomial[n,r]},{j,Binomial[n,r]}];


(* ::Input::Initialization:: *)
OnesFull[n_]:=SparseArray[Table[Band[{f[n,r],f[n,r]}]->OnesBlock[n,r],{r,0,n}]];


(* ::Input::Initialization:: *)
OnesFullk[n_,k_]:=
SparseArray[
Table[

If[k>=0,
Band[{f[n,r],f[n,r+k]}]->Table[1,{i,Binomial[n,r]},{j,Binomial[n,r+k]}],
Band[{f[n,r+Abs[k]],f[n,r]}]->Table[1,{i,Binomial[n,r+Abs[k]]},{j,Binomial[n,r]}]],

{r,0,n-Abs[k]}
],
{2^n,2^n}
]/;-n<=k<=n;


(* ::Input::Initialization:: *)
SuperBasisPositions[n_,k_]:=(SuperBasisPositions[n,k]=Position[OnesFullk[n,k]//Vectorise//Normal,1]//Flatten)/;-n<=k<=n;


(* ::Input::Initialization:: *)
SuperBlocksGet[\[ScriptCapitalL]_,n_,k_]:=(\[ScriptCapitalL][[SuperBasisPositions[n,k]]]\[Transpose][[SuperBasisPositions[n,k]]])\[Transpose];

(* ::Input::Initialization:: *)
SxSxCorrelations[\[Rho]_]:=Module[
{
SigmaxSite,
n=Log[2,\[Rho]//Length]
},
SigmaxSite=Table[Tensor@@Table[Piecewise[{{Sigmax,i==j},{Id,True}}],{i,n}],{j,n}];
Return[Table[Expect[SigmaxSite[[1]] . SigmaxSite[[i]],\[Rho]],{i,1,n}]]
]


(* ::Input::Initialization:: *)
SxSxDualityPlot[rhossDualityList_]:=Block[
{
n=Log[2,(rhossDualityList[[1,1,1]]//Dimensions)[[1]]],SxSxCorrList
},
SxSxCorrList=Map[SxSxCorrelations,rhossDualityList,{3}];

GraphicsGrid[Table[ListPlot[SxSxCorrList[[All,i,j]]//Re,Joined->True,PlotRange->{All,{-1,1}},AxesLabel->{"i","\!\(\*FormBox[\(\[LeftAngleBracket]\*SuperscriptBox[SubscriptBox[\(\[Sigma]\), \(1\)], \(x\)] \*SuperscriptBox[SubscriptBox[\(\[Sigma]\), \(i\)], \(x\)]\[RightAngleBracket]\),
TraditionalForm]\)"},Ticks->{Range[1,n],Automatic},PlotStyle->Thick,Mesh->All,PlotMarkers->Automatic],{i,2},{j,2}]
,ImageSize->Large]
];


(* ::Input::Initialization:: *)
EigOpBrute[{evals_,evecs_},J_]:=Block[
{
n=evals//Length,
SigmamSite,
EigOpSite,
bareEigOps
},
(* Initialise the Subscript[\[Sigma], i]^- operators for all i \[Element] [1, n]. *)
SigmamSite=Table[Tensor@@Table[Piecewise[{{Sigmam,i==j},{Id,True}}],{i,Log[2,n]}],{j,Log[2,n]}];

(* These are the projection operators of the eigenstates given. *)
bareEigOps=ParallelMap[SparseArray[Outer[Times,#,#]]&,evecs];

(* Enumeration of all operators \!\(
\*SubscriptBox[\(\[Sum]\), \(\[Epsilon], \[Nu]\)]\(J\((\[Epsilon] - \[Nu])\)
\*SubscriptBox[\(\[CapitalPi]\), \(\[Epsilon]\)]
\*SuperscriptBox[
SubscriptBox[\(\[Sigma]\), \(i\)], \(-\)]
\*SubscriptBox[\(\[CapitalPi]\), \(v\)]\ for\ each\ i\)\) \[Element] [1, n]. *)
EigOpSite=Total[Flatten[ParallelTable[J[evals[[i]]-evals[[j]]] bareEigOps[[i]] . # . bareEigOps[[j]],{i,n},{j,n}],1]]&/@SigmamSite;

Return[EigOpSite]
(* The super operator corresponding to the sum \!\(
\*SubscriptBox[\(\[Sum]\), \(i\)]\((
\*SuperscriptBox[
SubscriptBox[\(\[Sigma]\), \(i\)], \(+\)]
\*SuperscriptBox[\(A\), \(i\)]\[Rho]\  + \ \[Rho]\ 
\*SuperscriptBox[
SuperscriptBox[\(A\), \(i\)], \(\[Dagger]\)]
\*SuperscriptBox[
SubscriptBox[\(\[Sigma]\), \(i\)], \(-\)]\  - \ 
\*SuperscriptBox[\(A\), \(i\)]\[Rho]\ 
\*SuperscriptBox[
SubscriptBox[\(\[Sigma]\), \(i\)], \(+\)]\  - \ 
\*SuperscriptBox[
SubscriptBox[\(\[Sigma]\), \(i\)], \(-\)]\[Rho]\ 
\*SuperscriptBox[
SuperscriptBox[\(A\), \(i\)], \(\[Dagger]\)])\)\). (You would need to add the unitary part -\[ImaginaryI] [H, \[Rho]] to find a full Liouvillian.) *)
]


(* ::Input::Initialization:: *)
BlochRedfieldBrute[{evals_,evecs_},J_]:=Module[
{
n=evals//Length,
SigmamSite,
EigOpSite,
bareEigOps
},
(* Initialise the Subscript[\[Sigma], i]^- operators for all i \[Element] [1, n]. *)
SigmamSite=Table[Tensor@@Table[Piecewise[{{Sigmam,i==j},{Id,True}}],{i,Log[2,n]}],{j,Log[2,n]}];

(* These are the projection operators of the eigenstates given. *)
bareEigOps=ParallelMap[SparseArray[Outer[Times,#,#]]&,evecs];

(* Enumeration of all operators \!\(
\*SubscriptBox[\(\[Sum]\), \(\[Epsilon], \[Nu]\)]\(J\((\[Epsilon] - \[Nu])\)
\*SubscriptBox[\(\[CapitalPi]\), \(\[Epsilon]\)]
\*SuperscriptBox[
SubscriptBox[\(\[Sigma]\), \(i\)], \(-\)]
\*SubscriptBox[\(\[CapitalPi]\), \(v\)]\ for\ each\ i\)\) \[Element] [1, n]. *)
EigOpSite=Total[Flatten[ParallelTable[J[evals[[i]]-evals[[j]]] bareEigOps[[i]] . # . bareEigOps[[j]],{i,n},{j,n}],1]]&/@SigmamSite;

(* The super operator corresponding to the sum \!\(
\*SubscriptBox[\(\[Sum]\), \(i\)]\((
\*SuperscriptBox[
SubscriptBox[\(\[Sigma]\), \(i\)], \(+\)]
\*SuperscriptBox[\(A\), \(i\)]\[Rho]\  + \ \[Rho]\ 
\*SuperscriptBox[
SuperscriptBox[\(A\), \(i\)], \(\[Dagger]\)]
\*SuperscriptBox[
SubscriptBox[\(\[Sigma]\), \(i\)], \(-\)]\  - \ 
\*SuperscriptBox[\(A\), \(i\)]\[Rho]\ 
\*SuperscriptBox[
SubscriptBox[\(\[Sigma]\), \(i\)], \(+\)]\  - \ 
\*SuperscriptBox[
SubscriptBox[\(\[Sigma]\), \(i\)], \(-\)]\[Rho]\ 
\*SuperscriptBox[
SuperscriptBox[\(A\), \(i\)], \(\[Dagger]\)])\)\). (You would need to add the unitary part -\[ImaginaryI] [H, \[Rho]] to find a full Liouvillian.) *)

Total[ParallelTable[
SuperLR[EigOpSite[[i]],SigmamSite[[i]]\[Transpose]]
+SuperLR[SigmamSite[[i]],EigOpSite[[i]]\[HermitianConjugate]]
-SuperLR[IdentityMatrix[{n,n}],EigOpSite[[i]]\[HermitianConjugate] . SigmamSite[[i]]]
-SuperLR[SigmamSite[[i]]\[Transpose] . EigOpSite[[i]],IdentityMatrix[{n,n}]],
{i,Log[2,n]}]]
]


(* ::Input::Initialization:: *)
ZeroDeltaEigOp[{evalsList_, evecsList_}, J_] := 
  Module[
   {
    outerList = Map[Outer[Times, #, #] &, evecsList, {2}],
    outerPairs,
    EigOpBlock,
    n = Length[evalsList] - 1,
    SigmamSitesBlocks,
    SigmamSite,
    U,
    EigOpSite
    },
   SigmamSitesBlocks = Table[MinusBlocks[n, j], {j, n}];
   
   SigmamSite = Table[Tensor @@ Table[\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{"Sigmam", 
RowBox[{"i", "==", "j"}]},
{"Id", "True"}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "ColumnsIndexed" -> {}, "Rows" -> {{1.}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "ColumnsIndexed" -> {}, "Rows" -> {{1.}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False]\) , {i, n}], {j, n}];
   
   U = TotalSpinBasis[n];
   
   outerPairs=Table[{{evalsList[[r]],outerList[[r]]},{evalsList[[r+1]],outerList[[r+1]]}}, {r,n}];

EigOpBlock[oper_,pair_]:=Total[Flatten[
      ParallelTable[J[pair[[1,1,i]]-pair[[2,1,j]]]pair[[1,2,i]] . oper . pair[[2,2,j]],
{i,Length[pair[[1,1]]]},{j,Length[pair[[2,1]]]}],1]];
   
  EigOpSite=U\[Transpose] . # . U&/@Table[Blockto2N[EigOpBlock[#[[1]],#[[2]]]&/@Table[{SigmamSitesBlocks[[j]][[i]],outerPairs[[i]]},{i,n}], -1], {j,n}];
   
   EigOpSite
];


(* ::Input::Initialization:: *)
BlochRedfieldZero\[CapitalDelta][{evalsList_,evecsList_},J_]:=Module[
{
outerList=Map[Outer[Times,#,#]&,evecsList,{2}],
outerPairs,
EigOpBlock,
n=Length[evalsList]-1,
SigmamSitesBlocks,
SigmamSite,
U,
EigOpSite
},
SigmamSitesBlocks=Table[MinusBlocks[n,j],{j,n}];

SigmamSite=Table[Tensor@@Table[\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{"Sigmam", 
RowBox[{"i", "==", "j"}]},
{"Id", "True"}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "ColumnsIndexed" -> {}, "Rows" -> {{1.}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "ColumnsIndexed" -> {}, "Rows" -> {{1.}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False]\),{i,n}],{j,n}];

U=TotalSpinBasis[n];

outerPairs=Table[{{evalsList[[r]],outerList[[r]]},{evalsList[[r+1]],outerList[[r+1]]}},{r,n}];

EigOpBlock[oper_,pair_]:=Total[Flatten[ParallelTable[J[pair[[1,1,i]]-pair[[2,1,j]]] pair[[1,2,i]] . oper . pair[[2,2,j]],{i,Length[pair[[1,1]]]},{j,Length[pair[[2,1]]]}],1]];

EigOpSite=U\[Transpose] . # . U&/@Table[Blockto2N[EigOpBlock[#[[1]],#[[2]]]&/@Table[{SigmamSitesBlocks[[j]][[i]],outerPairs[[i]]},{i,n}],-1],{j,n}];

Total[ParallelTable[
SuperLR[EigOpSite[[i]],SigmamSite[[i]]\[Transpose]]
+SuperLR[SigmamSite[[i]],EigOpSite[[i]]\[HermitianConjugate]]
-SuperLR[IdentityMatrix[{2^n,2^n}]//SparseArray,EigOpSite[[i]]\[HermitianConjugate] . SigmamSite[[i]]]
-SuperLR[SigmamSite[[i]]\[Transpose] . EigOpSite[[i]],IdentityMatrix[{2^n,2^n}]//SparseArray],
{i,n}]]
];


(* ::Input::Initialization:: *)
EigVecKet[r_]:=Which[r-2<0,Subscript[\[Epsilon], r]+Subscript[P, r+2],r+2>n, Subscript[\[Epsilon], r]+Subscript[M, r-2],True,Subscript[\[Epsilon], r]+Subscript[M, r-2]+Subscript[P, r+2]]/;0<=r<=n;

(* & h.c. *)
EigVecBra[r_]:=Which[r-2<0,Subscript[\[Epsilon]T, r]+Subscript[PT, r+2],r+2>n, Subscript[\[Epsilon]T, r]+Subscript[MT, r-2],True,Subscript[\[Epsilon]T, r]+Subscript[MT, r-2]+Subscript[PT, r+2]]/;0<=r<=n;

(* The expansion of the outer products are encoded as the i,jth elements of the 'matrix of matrices' A\[Wedge]B *)
EigOpZero[r_,rp_]:=(Subscript[\[Epsilon], r]\[Wedge]Subscript[\[Epsilon]T, rp])/;(0<=r<=n&&0<=rp<=n);

EigOpFirst[r_,rp_]:=If[rp+2<=n,Subscript[\[Epsilon], r]\[Wedge]Subscript[PT, rp+2],0]+If[rp-2>=0,Subscript[\[Epsilon], r]\[Wedge]Subscript[MT, rp-2],0]+If[r+2<=n,Subscript[P, r+2]\[Wedge]Subscript[\[Epsilon]T, rp],0]+If[r-2>=0,Subscript[M, r-2]\[Wedge]Subscript[\[Epsilon]T, rp],0];


EigOpZero[r_,rp_]:=(Subscript[\[Epsilon], r]\[Wedge]Subscript[\[Epsilon]T, rp])/;(0<=r<=n&&0<=rp<=n);

EigOpFirst[r_,rp_]:=If[rp+2<=n,Subscript[\[Epsilon], r]\[Wedge]Subscript[PT, rp+2],0]+If[rp-2>=0,Subscript[\[Epsilon], r]\[Wedge]Subscript[MT, rp-2],0]+If[r+2<=n,Subscript[P, r+2]\[Wedge]Subscript[\[Epsilon]T, rp],0]+If[r-2>=0,Subscript[M, r-2]\[Wedge]Subscript[\[Epsilon]T, rp],0]/;(0<=r<=n&&0<=rp<=n);

LeadingOrderQ[A_,B_]:=(A===\[Epsilon]T&&B===\[Epsilon])||(A===\[Epsilon]T&&B===M)||(A===\[Epsilon]T&&B===P)||(A===MT&&B===\[Epsilon])||(A===PT&&B===\[Epsilon]);

FirstOrderQ[elem_]:=(elem/.{Subscript[PT, _]->0,Subscript[MT, _]->0,Subscript[P, _]->0,Subscript[M, _]->0}//TensorExpand)===0;

PertVects[{eEner_,eVec_}]:=Block[
{
n=Length[eEner]-1,
\[CapitalXi]iiP2,Vmm,\[CapitalXi]iiM2,Vpp,mi,pi,P2toN,M0toNm2,M,P,Pr,Mr,\[CapitalXi]iiP2denom,\[CapitalXi]iiM2denom
},

\[CapitalXi]iiP2=Table[{eVec[[i]],eVec[[i+2]]},{i,1,n-1}];

(* MINUS SIGN! *)
\[CapitalXi]iiP2denom=-Table[Table[1/(eEner[[i]][[i0]]-eEner[[i+2]][[i2]]),{i0,Binomial[n,i-1]},{i2,Binomial[n,i+1]}],{i,1,n-1}];

Vmm=Normal/@MinusMinusBlocks[n,1];

\[CapitalXi]iiM2=Table[{eVec[[i+2]],eVec[[i]]},{i,1,n-1}];

(* MINUS SIGN! *)
\[CapitalXi]iiM2denom=-Table[Table[1/(eEner[[i+2]][[i2]]-eEner[[i]][[i0]]),{i2,Binomial[n,i+1]},{i0,Binomial[n,i-1]}],{i,1,n-1}];

Vpp=Normal[#\[Transpose]]&/@MinusMinusBlocks[n,1];
(*
Print[Map[Dimensions,\[CapitalXi]iiP2,{2}]];
Print[Dimensions/@Vmm];
*)
(* Perturbative corrections *)
M=MapThread[#1[[1]] . #2 . #1[[2]]\[Transpose]&,{\[CapitalXi]iiP2,Vmm}];
P=MapThread[#1[[1]] . #2 . #1[[2]]\[Transpose]&,{\[CapitalXi]iiM2,Vpp}];

M=Table[M[[i]]*\[CapitalXi]iiP2denom[[i]],{i,Length[M]}];(* Why aren't these matched??*)
P=Table[P[[i]]*\[CapitalXi]iiM2denom[[i]],{i,Length[P]}];(* Yeah still confused *)
(* Rmemeber the demonimator !!!! *)

(* Corrective VECTORS *)
(* 
\[LeftAngleBracket]r-2,Subscript[j, r-2]|Subsuperscript[\[ScriptCapitalP], Subscript[i, r], (r+2)]\[RightAngleBracket] = \!\(
\*SubscriptBox[\(\[Sum]\), 
SubscriptBox[\(i\), \(r + 2\)]]
\*SubsuperscriptBox[\(P\), \(
\*SubscriptBox[\(i\), \(r\)]
\*SubscriptBox[\(i\), \(r + 2\)]\), \((r + 2)\)]\)|r+2,Subscript[i, r+2]\[RightAngleBracket] = Subscript[[(\[ScriptCapitalP]^(r+2))\[Transpose].E^(r+2)], Subscript[i, r]Subscript[j, r+2]]
*)
(*  AND  *)
(*
\[LeftAngleBracket]r-2,Subscript[j, r-2]|Subsuperscript[\[ScriptCapitalM], Subscript[i, r], (r-2)]\[RightAngleBracket] = \!\(
\*SubscriptBox[\(\[Sum]\), 
SubscriptBox[\(i\), \(r - 2\)]]
\*SubsuperscriptBox[\(M\), \(
\*SubscriptBox[\(i\), \(r\)]
\*SubscriptBox[\(i\), \(r - 2\)]\), \((r - 2)\)]\)|r-2,Subscript[i, r-2]\[RightAngleBracket] = Subscript[[(\[ScriptCapitalM]^(r-2))\[Transpose].E^(r-2)], Subscript[i, r]Subscript[j, r-2]]
*)
P2toN=Table[P[[i]]\[Transpose] . eVec[[i+2]],{i,n-1}];
Pr={0,0}~Join~P2toN;

M0toNm2=Table[M[[i]]\[Transpose] . eVec[[i]],{i,n-1}];
Mr=M0toNm2~Join~{0,0};

Return[{Pr,Mr}]
];


(* ::Input::Initialization:: *)
EigOpPerturbative[{eEner_,eVec_},\[CapitalDelta]_,J_]:=Block[
{
n=Length[eEner]-1,
SigmaMinusChain,
InnerMatrices,
EigOps,
\[Sigma]mblocks,
Pr,Mr,
\[Epsilon]v,\[Sigma]M,\[Sigma]P,InnerRules,
Operator,OperatorList,MatElem,ZeroOrderRules,FirstOrderRules,
CurrMatElem,OuterFirstOrderRules,OuterSum,
CurrMatElemZero,CurrMatElemFirst,U
},

\[Sigma]mblocks=Table[MinusBlocks[n,i],{i,1,n}];

InnerMatrices=Table[EigVecBra[r] . Subscript[\[Sigma]m, i] . EigVecKet[rp],{r,0,n},{rp,0,n}]//TensorExpand;
EigOps=Table[Which[
rp==r+1,EigOpFirst[r-1,rp-1]+EigOpZero[r-1,rp-1],
rp==r+3,EigOpZero[r-1,rp-1],
rp==r-1,EigOpZero[r-1,rp-1],
              True,0
	],
	{r,1,n+1},{rp,1,n+1}];

(* All those matrix elements that are zero to first order are removed. *)
InnerMatrices=InnerMatrices/.{Subscript[A_, q_] . Subscript[\[Sigma]m, i] . Subscript[B_, p_]:>If[LeadingOrderQ[A,B]&&p==q+1,Subscript[A, q] . Subscript[\[Sigma]m, i] . Subscript[B, p],0]};

{Pr,Mr}=PertVects[{eEner,eVec}];

OperatorList={};
Do[
ZeroOrderRules=Subscript[\[Epsilon]T, r_] . Subscript[\[Sigma]m, i] . Subscript[\[Epsilon], s_]:>eVec[[r+1]] . \[Sigma]mblock[[r+1]] . eVec[[s+1]]\[Transpose];
FirstOrderRules={
Subscript[\[Epsilon]T, r_] . Subscript[\[Sigma]m, i] . Subscript[P, s_]:>eVec[[r+1]] . \[Sigma]mblock[[r+1]] . Pr[[s+1]]\[Transpose],
Subscript[\[Epsilon]T, r_] . Subscript[\[Sigma]m, i] . Subscript[M, s_]:>eVec[[r+1]] . \[Sigma]mblock[[r+1]] . Mr[[s+1]]\[Transpose],
Subscript[MT, r_] . Subscript[\[Sigma]m, i] . Subscript[\[Epsilon], s_]:>Mr[[r+1]] . \[Sigma]mblock[[r+1]] . eVec[[s+1]]\[Transpose],
Subscript[PT, r_] . Subscript[\[Sigma]m, i] . Subscript[\[Epsilon], s_]:>Pr[[r+1]] . \[Sigma]mblock[[r+1]] . eVec[[s+1]]\[Transpose]
};

OuterSum[0]=0.;
OuterSum[Subscript[P, r_]\[Wedge]Subscript[\[Epsilon]T, s_],eEnR_,eEnS_,MatElem_]:=Sum[J[eEnR[[ir]]-eEnS[[is]]]*MatElem[[ir,is]]*SparseArray[{Band[{f[n,r],f[n,s]}]->Outer[Times,Pr[[r+1,ir]],eVec[[s+1,is]]]},{2^n,2^n}],{ir,1,Binomial[n,r-2]},{is,1,Binomial[n,s]}];

OuterSum[Subscript[M, r_]\[Wedge]Subscript[\[Epsilon]T, s_],eEnR_,eEnS_,MatElem_]:=Sum[J[eEnR[[ir]]-eEnS[[is]]]*MatElem[[ir,is]]*SparseArray[{Band[{f[n,r],f[n,s]}]->Outer[Times,Mr[[r+1,ir]],eVec[[s+1,is]]]},{2^n,2^n}],{ir,1,Binomial[n,r+2]},{is,1,Binomial[n,s]}];

OuterSum[Subscript[\[Epsilon], r_]\[Wedge]Subscript[PT, s_],eEnR_,eEnS_,MatElem_]:=Sum[J[eEnR[[ir]]-eEnS[[is]]]*MatElem[[ir,is]]*SparseArray[{Band[{f[n,r],f[n,s]}]->Outer[Times,eVec[[r+1,ir]],Pr[[s+1,is]]]},{2^n,2^n}],{ir,1,Binomial[n,r]},{is,1,Binomial[n,s-2]}];

OuterSum[Subscript[\[Epsilon], r_]\[Wedge]Subscript[MT, s_],eEnR_,eEnS_,MatElem_]:=Sum[J[eEnR[[ir]]-eEnS[[is]]]*MatElem[[ir,is]]*SparseArray[{Band[{f[n,r],f[n,s]}]->Outer[Times,eVec[[r+1,ir]],Mr[[s+1,is]]]},{2^n,2^n}],{ir,1,Binomial[n,r]},{is,1,Binomial[n,s+2]}];

Operator=SparseArray[{},{2^n,2^n}]; SetSharedVariable[Operator];

Do[
MatElem=InnerMatrices[[r+1,s+1]];

If[MatElem==0,Continue[]];

CurrMatElemZero=MatElem/.ZeroOrderRules;
CurrMatElemFirst=MatElem/.FirstOrderRules;

If[!FirstOrderQ[MatElem],
Operator+=Sum[J[eEner[[r+1,ir]]-eEner[[s+1,is]]]*
CurrMatElemZero[[ir,is]]*SparseArray[{Band[{f[n,r],f[n,s]}]->Outer[Times,eVec[[r+1,ir]],eVec[[s+1,is]]]},{2^n,2^n}],
{ir,1,Binomial[n,r]},{is,1,Binomial[n,s]}
];

Operator+=\[CapitalDelta] OuterSum[#,eEner[[r+1]],eEner[[s+1]],MatElem/.ZeroOrderRules]&/@(EigOps[[r+1,s+1]]/.{Subscript[\[Epsilon], r]\[Wedge]Subscript[\[Epsilon]T, s]->0});
Continue[]
];

Operator+=\[CapitalDelta] Sum[J[eEner[[r+1,ir]]-eEner[[s+1,is]]]*
CurrMatElemFirst[[ir,is]]*SparseArray[{Band[{f[n,r],f[n,s]}]->Outer[Times,eVec[[r+1,ir]],eVec[[s+1,is]]]},{2^n,2^n}],
{ir,1,Binomial[n,r]},{is,1,Binomial[n,s]}
];,

{r,0,n},{s,0,n}];

OperatorList=Append[OperatorList,Operator];,
{\[Sigma]mblock, \[Sigma]mblocks}];

U=TotalSpinBasis[n];
Return[U\[Transpose] . # . U&/@OperatorList]

];


(* ::Input::Initialization:: *)
BlochRedfieldPert[eSys_,Delta_,J_]:=Block[
{
n=(eSys[[1]]//Length)-1,SigmamSite,EigOpSite=EigOpPerturbative[eSys,Delta,J]
},

SigmamSite=Table[Tensor@@Table[Piecewise[{{Sigmam,i==j},{Id,True}}],{i,n}],{j,n}];

Total[
Table[
SuperLR[EigOpSite[[i]],SigmamSite[[i]]\[Transpose]]
+SuperLR[SigmamSite[[i]],EigOpSite[[i]]\[HermitianConjugate]]
-SuperLR[IdentityMatrix[{2^n,2^n}],EigOpSite[[i]]\[HermitianConjugate] . SigmamSite[[i]]]
-SuperLR[SigmamSite[[i]]\[Transpose] . EigOpSite[[i]],IdentityMatrix[{2^n,2^n}]],
{i,n}]
]
]


(* ::Input::Initialization:: *)
\[Rho]ss[nn_,g_,\[CapitalDelta]_,J_,option_]:=Module[
{myH,eigS,SPspec,SPbndwid,my\[ScriptCapitalD],my\[ScriptCapitalL],\[Rho]out,eval},

myH=HAnIsoXY[nn,g,\[CapitalDelta]];
eigS=Eigensystem[myH];

Which[
option=="BruteForce",
my\[ScriptCapitalD]=BlochRedfieldBrute[Eigensystem[myH],J],

option=="UnperturbedDissipator",
my\[ScriptCapitalD]=BlochRedfieldZero\[CapitalDelta][(Eigensystem/@H0Blocks[nn,g])\[Transpose],J],

option=="Perturbative",
my\[ScriptCapitalD]=BlochRedfieldPert[(Eigensystem/@H0Blocks[nn,g])\[Transpose],\[CapitalDelta],J] ,

option=="Lindblad",
my\[ScriptCapitalD]=J Total[Lindblad/@Table[Tensor@@Table[Piecewise[{{Sigmam,i==j},{Id,True}}],{i,nn}],{j,nn}]]
];


(* Overscript[\[Rho], .] = -\[ImaginaryI][H, \[Rho]] + \[ScriptCapitalD]\[Rho] ... *)
my\[ScriptCapitalL]=-I SuperComm[myH]+my\[ScriptCapitalD];

If[nn<=6,
{eval,\[Rho]out}=Eigensystem[my\[ScriptCapitalL]//Threshold,-1,Method->{"Arnoldi"}]\[Transpose][[1]],
{eval,\[Rho]out}=Eigensystem[my\[ScriptCapitalL]//Threshold,1,Method->{"Arnoldi","MaxIterations"->30000,"Criteria"->"RealPart"}]\[Transpose][[1]]
];

\[Rho]out=toHilbert[\[Rho]out];

(*Print["Checking that state is steady. Norm[\[ScriptCapitalL] \[Rho]] = ",(my\[ScriptCapitalL].Vectorise[\[Rho]out])//Norm," | Found Eigenvalue = ",eval];*)

Return[\[Rho]out/Tr[\[Rho]out]];
];
 
(* ::Input::Initialization:: *)
SteadyTimeEvo[n_,g_,\[CapitalDelta]_,{Jname_,J_},method_,tol_:10^-7]:=Module[
{
stateData,
\[Rho]init,\[Rho]and\[Rho]dot,
H,\[ScriptCapitalL],\[ScriptCapitalD],

fileNameBase=StringJoin[ToString/@({"[N=",n,"][g=",g,"][\[CapitalDelta]=",\[CapitalDelta],"][method=",method,"]"})]<>"["<>Jname<>"]",

optFlags,
CheckpointFlag=Checkpoint->True,
MaxNormFlag=MaxNorm->False,

processEq=NDSolve`ProcessEquations,
iterate=NDSolve`Iterate,
processSol=NDSolve`ProcessSolutions,
reInit=NDSolve`Reinitialize,

evolveState,tstep
},

tstep=0.1;
H=HAnIsoXY[n,g,\[CapitalDelta]];

(* Liouvillian *)
Print["Data/L_"<>fileNameBase];
If[
!FileExistsQ["Data/L_"<>fileNameBase<>".mtx"],

\[ScriptCapitalD]=Which[
method=="BruteForce",
BlochRedfieldBrute[Eigensystem[H],J],

method=="UnperturbedDissipator",
BlochRedfieldZero\[CapitalDelta][(Eigensystem/@H0Blocks[n,g])\[Transpose],J],

method=="PerturbationTheory",
BlochRedfieldPert[(Eigensystem/@H0Blocks[n,g])\[Transpose],\[CapitalDelta],J]
];

\[ScriptCapitalL]=-I SuperComm[H]+\[ScriptCapitalD];
Export["Data/L_"<>fileNameBase<>".mtx",\[ScriptCapitalL],"MTX"],

(* Else *)
\[ScriptCapitalL]=Import["Data/L_"<>fileNameBase<>".mtx","MTX"]
];

(* Initial density matrix *)
If[
!FileExistsQ["Data/rhoSS_"<>fileNameBase<>".mat"],

(* \[Rho] = \!\(
\*SubscriptBox[\(\[Sum]\), \(i\)]
\*SubscriptBox[\(\[Lambda]\), \(i\)]\)|Subscript[\[Psi], i]\[RightAngleBracket]\[LeftAngleBracket]Subscript[\[Psi], i]|,  {Subscript[\[Psi], i] \[Element] Configuration basis} *)
\[Rho]init=Flatten[2.^-n Subscript[\[ScriptCapitalI], 2^n]];Export["Data/rhoSS_"<>fileNameBase<>".mat",\[Rho]init],

(* Else *)
\[Rho]init=Import["Data/rhoSS_"<>fileNameBase<>".mat"][[1]];
Print[Norm[\[ScriptCapitalL] . SparseArray[\[Rho]init]]];
];

(* ------------------------------------------------------ *)

(* Steady state by time-evolution *)
stateData=processEq[{\[Rho]'[t]==\[ScriptCapitalL] . \[Rho][t],\[Rho][0]==\[Rho]init},\[Rho],{t,0,\[Infinity]},Method->"ExplicitRungeKutta"][[1]];

(* Basic evolution block *)
evolveState[t_]:=
(
iterate[stateData,t];
\[Rho]and\[Rho]dot=processSol[stateData,"Forward"];
stateData=reInit[stateData,Equal@@@Most[\[Rho]and\[Rho]dot]][[1]];
);

evolveState[tstep];
NestWhile[
(
evolveState[tstep*#];
#+1
)&,1,

Norm[\[ScriptCapitalL] . SparseArray[\[Rho]and\[Rho]dot[[1,2]]]]>tol&
];

Export["Data/rhoSS_"<>fileNameBase<>".mat",\[Rho]and\[Rho]dot[[1,2]]];
Return[\[Rho]and\[Rho]dot[[1,2]]];
]


(* ::Input::Initialization:: *)
SteadyStateTE[\[ScriptCapitalL]_,tol_:10^-7]:=Module[
{

n=Log[4,Length[\[ScriptCapitalL]]],
stateData,
\[Rho]init,\[Rho]and\[Rho]dot,

optFlags,
CheckpointFlag=Checkpoint->True,
MaxNormFlag=MaxNorm->False,

processEq=NDSolve`ProcessEquations,
iterate=NDSolve`Iterate,
processSol=NDSolve`ProcessSolutions,
reInit=NDSolve`Reinitialize,

evolveState,tstep
},

tstep=0.1;

(* Initial density matrix *)
(* \[Rho] = \!\(
\*SubscriptBox[\(\[Sum]\), \(i\)]
\*SubscriptBox[\(\[Lambda]\), \(i\)]\)|Subscript[\[Psi], i]\[RightAngleBracket]\[LeftAngleBracket]Subscript[\[Psi], i]|,  {Subscript[\[Psi], i] \[Element] Configuration basis} *)
\[Rho]init=Flatten[2.^-n Subscript[\[ScriptCapitalI], 2^n]];

(* ------------------------------------------------------ *)

(* Steady state by time-evolution *)
stateData=processEq[{\[Rho]'[t]==\[ScriptCapitalL] . \[Rho][t],\[Rho][0]==\[Rho]init},\[Rho],{t,0,\[Infinity]},Method->"ExplicitRungeKutta"][[1]];

(* Basic evolution block *)
evolveState[t_]:=
(
iterate[stateData,t];
\[Rho]and\[Rho]dot=processSol[stateData,"Forward"];
stateData=reInit[stateData,Equal@@@Most[\[Rho]and\[Rho]dot]][[1]];
);

evolveState[tstep];
NestWhile[
(
evolveState[tstep*#];
#+1
)&,1,

Norm[\[ScriptCapitalL] . SparseArray[\[Rho]and\[Rho]dot[[1,2]]]]>tol&
];

Return[\[Rho]and\[Rho]dot[[1,2]]];
]


(* ::Input::Initialization:: *)
Protect[Overwrite];

SaveDataCall[fTemplate_,ArgList_,flags___]:=Block[
{
FileName,
FunctionName,
ArgSymbs,
CallRules,
ArgsEq,
FlagsDefault,OverwriteFlag,
out
},

FlagsDefault={Overwrite->False};
{OverwriteFlag}={Overwrite}/.{flags}/.FlagsDefault;

If[Head[fTemplate]===HoldForm,

{FunctionName,ArgSymbs}={#[[1,0]],{ReleaseHold@@(Hold@@@#)}}&@fTemplate,

FunctionName=Head[fTemplate];ArgSymbs=List@@fTemplate
];

CallRules=Thread[ArgSymbs->ArgList];

ArgsEq=(ToString[#]<>"=")&/@ArgSymbs;
FileName=ToString[FunctionName]<>"("<>Riffle[ArgsEq,ToString/@ArgList]<>").mat";

If[
!FileExistsQ[FileName]||OverwriteFlag,
out=Reap[Export[FileName,Sow[fTemplate/.CallRules//ReleaseHold],"MTX"]][[2,1,1]],
out=Import[FileName,"MTX"]
];
out
];

SaveDataCallRange[f_,ArgsLists_,flags___]:=Block[
{
ArgSymbs=ArgsLists[[;;,1]], (* Labels defining each of the ranges*)
FileNamesFromArgs,
AllFileNames,
EvaluateWithArgs,
FlagsList,ArgsEq,
NonExistentFiles,
AllArgsValues,
TuplesToCalculate,
FlagsDefault,
OverwriteFlag
},

(* Flag handling. *)
FlagsDefault={Overwrite->False};
{OverwriteFlag}={Overwrite}/.{flags}/.FlagsDefault;

ArgsEq=(ToString[#]<>"=")&/@ArgSymbs;

FileNamesFromArgs[Args__]:=Table[ToString[f]<>"("<>Riffle[ArgsEq,ToString/@ArgSymbs]<>").mat",Args];
AllFileNames=FileNamesFromArgs@@ArgsLists;

Print[AllFileNames];

NonExistentFiles=Select[Flatten[AllFileNames],Not@*FileExistsQ];

(* Find the positions in the full list which have not been calculated and calculate the tuples of arguments that these positions correspond to. *)
AllArgsValues=Table[ArgSymbs,##]&@@ArgsLists;

TuplesToCalculate=Part[AllArgsValues,##]&@@@(Position[AllFileNames,#][[1]]&/@NonExistentFiles);

If[OverwriteFlag,
MapThread[(Export[#1,f@@#2,"MTX"])&,{AllFileNames,AllArgsValues},2],
MapThread[(Export[#1,f@@#2,"MTX"])&,{NonExistentFiles,TuplesToCalculate}]
];
];


(* ::Input::Initialization:: *)
GeneralFileName[n_,g_,\[CapitalDelta]_,method_,Jname_]:=StringJoin[ToString/@({"[N=",n,"][g=",g,"][\[CapitalDelta]=",\[CapitalDelta],"][method=",method,"]"})]<>"["<>Jname<>"]"


(* ::Input::Initialization:: *)
\[Rho]SSData[dat_,tmax_:0.,normL\[Rho]_:\[Infinity]]:=<|"dat"->dat,"tmax"->tmax,"normL\[Rho]"->normL\[Rho]|>;

\[Rho]SSDataFluc[dat_,tmax_:0.,normL\[Rho]_:\[Infinity],Fluc_:{{},{},{}}]:=<|"dat"->dat,"tmax"->tmax,"normL\[Rho]"->normL\[Rho],"FDT"->Fluc|>


SteadyStateDataInit[Method_,n_,g_,\[CapitalDelta]_,\[Rho]SS_:"NA",tmax_:0.,normL\[Rho]_:\[Infinity]]:=<|"Method"->Method,"n"->n,"g"->g,"\[CapitalDelta]"->\[CapitalDelta],"\[Rho]SS"->\[Rho]SSData[#, tmax, normL\[Rho]]&@If[\[Rho]SS==="NA",HAnIsoXY[n,g,\[CapitalDelta]],\[Rho]SS]|>


SteadyStateDataInitFluc[Method_,n_,g_,\[CapitalDelta]_,\[Rho]SS_:"NA",tmax_:0.,normL\[Rho]_:\[Infinity],Fluc_:{{},{},{}}]:=<|"Method"->Method,"n"->n,"g"->g,"\[CapitalDelta]"->\[CapitalDelta],"\[Rho]SS"->\[Rho]SSDataFluc[#, tmax, normL\[Rho],Fluc]&@If[\[Rho]SS==="NA",HAnIsoXY[n,g,\[CapitalDelta]],\[Rho]SS]|>


SteadyStateDataInitFluc["te", 3, 1., 1.]


ParamStr[{NameStr_, Val_}] := "[" <> NameStr <> "=" <> ToString[Val] <> "]";

ParamListString[ExperimentName_, ParamStrList_] := StringJoin[ExperimentName, ##, ".m"]& @@ (ParamStr /@ ParamStrList);


(* ::Input::Initialization:: *)
AdiabaticSweep[LiouvillianFunc_,{varMin_,varMax_,varStep_},flags___]:=
Module[
{
stateData,
\[Rho]init,\[Rho]and\[Rho]dot,
H,\[ScriptCapitalL],\[ScriptCapitalD],

FlagsDefault,
tol,OverwriteFlag,

optFlags,
CheckpointFlag=Checkpoint->True,
MaxNormFlag=MaxNorm->False,

processEq=NDSolve`ProcessEquations,
iterate=NDSolve`Iterate,
processSol=NDSolve`ProcessSolutions,
reInit=NDSolve`Reinitialize,

evolveState,tstep=0.1,res
},

(* FLAG HANDLING *)
FlagsDefault={MaxNorm->10^-7,Overwrite->False};
{tol,OverwriteFlag}={MaxNorm,Overwrite}/.{flags}/.FlagsDefault;


(* INITIALISE SOLVER *)
\[ScriptCapitalL]=LiouvillianFunc[varMin];\[Rho]init=Flatten[2.^-n Subscript[\[ScriptCapitalI], 2^n]];
stateData=processEq[{\[Rho]'[t]==\[ScriptCapitalL] . \[Rho][t],\[Rho][0]==\[Rho]init},\[Rho],{t,0,\[Infinity]},Method->"ExplicitRungeKutta"][[1]];


(* Basic evolution block *)
evolveState[t_]:=
(
iterate[stateData,t];
\[Rho]and\[Rho]dot=processSol[stateData,"Forward"];
stateData=reInit[stateData,Equal@@@Most[\[Rho]and\[Rho]dot]][[1]];
);

(* LOOP OVER SWEEP VARIABLES *)
res=Reap[Do[

evolveState[tstep];
NestWhile[
(
evolveState[tstep*#];
#+1
)&,1,

(#<15||Norm[\[ScriptCapitalL] . SparseArray[\[Rho]and\[Rho]dot[[1,2]]]]>tol)&
];

Sow[{var, \[Rho]and\[Rho]dot[[1,2]]//toHilbert//SparseArray}];

\[ScriptCapitalL]=LiouvillianFunc[var];
stateData=processEq[{\[Rho]'[t]==\[ScriptCapitalL] . \[Rho][t],\[Rho][0]==\[Rho]and\[Rho]dot[[1,2]]},\[Rho],{t,0,\[Infinity]},Method->"ExplicitRungeKutta"][[1]];


(* END DO *)
,{var,varMin,varMax,varStep}];
];
Return[res[[2,1]]];
]


(* ::Input::Initialization:: *)
TwoTimeCorrelator[\[Rho]ss_,\[ScriptCapitalL]_,X_,Y_,tstep_,tmax_]:=Module[
{
stateDataL,stateDataR,
\[Rho]init,\[Rho]and\[Rho]dotL,\[Rho]and\[Rho]dotR,
H,\[ScriptCapitalD],

FlagsDefault,
tol,OverwriteFlag,

optFlags,
CheckpointFlag=Checkpoint->True,
MaxNormFlag=MaxNorm->False,

processEq=NDSolve`ProcessEquations,
iterate=NDSolve`Iterate,
processSol=NDSolve`ProcessSolutions,
reInit=NDSolve`Reinitialize,
processreInit,
evolveState,res,tCurr=0.,i
},

stateDataL=processEq[{\[Rho]'[t]==\[ScriptCapitalL] . \[Rho][t],\[Rho][0]==SuperLR[Y,Subscript[\[ScriptCapitalI], Y//Length]] . Vectorise[\[Rho]ss]},\[Rho],{t,0,\[Infinity]},Method->"ExplicitRungeKutta"][[1]];
stateDataR=processEq[{\[Rho]'[t]==\[ScriptCapitalL] . \[Rho][t],\[Rho][0]==SuperLR[Subscript[\[ScriptCapitalI], Y//Length],Y] . Vectorise[\[Rho]ss]},\[Rho],{t,0,\[Infinity]},Method->"ExplicitRungeKutta"][[1]];

(* Basic evolution block *)
processreInit[]:=(
\[Rho]and\[Rho]dotL=processSol[stateDataL,"Forward"];
\[Rho]and\[Rho]dotR=processSol[stateDataR,"Forward"];
stateDataL=reInit[stateDataL,Equal@@@Most[\[Rho]and\[Rho]dotL]][[1]];
stateDataR=reInit[stateDataR,Equal@@@Most[\[Rho]and\[Rho]dotR]][[1]];);

evolveState[t_]:=
(
iterate[stateDataL,t];iterate[stateDataR,t];processreInit[]
);

i=1;
Reap[

processreInit[];
Sow[{tCurr,{Tr[X . toHilbert[\[Rho]and\[Rho]dotL[[1,2]]]],Tr[X . toHilbert[\[Rho]and\[Rho]dotL[[1,2]]]+X . toHilbert[\[Rho]and\[Rho]dotR[[1,2]]]]/2,
I Tr[X . toHilbert[\[Rho]and\[Rho]dotL[[1,2]]]-X . toHilbert[\[Rho]and\[Rho]dotR[[1,2]]]]}}];

While[tCurr<=tmax,
evolveState[tstep*i]; (*Evolve forward to t = tstep * i*)
i++;tCurr+=tstep;

(*
	S(t) = \[LeftAngleBracket]{X(t), X(0)}\[RightAngleBracket]/2 = Tr[X(t).Y(0).\[Rho] + Y(0).X(t).\[Rho]]/2 ;
	\[Chi](t) = \[ImaginaryI]\[LeftAngleBracket][X(t), X(0)]\[RightAngleBracket] = \[ImaginaryI] Tr[X(t).Y(0).\[Rho] - Y(0).X(t).\[Rho]]
*)

Sow[{tCurr,{Tr[X . toHilbert[\[Rho]and\[Rho]dotL[[1,2]]]],Tr[X . toHilbert[\[Rho]and\[Rho]dotL[[1,2]]]+X . toHilbert[\[Rho]and\[Rho]dotR[[1,2]]]]/2,
I Tr[X . toHilbert[\[Rho]and\[Rho]dotL[[1,2]]]-X . toHilbert[\[Rho]and\[Rho]dotR[[1,2]]]]}}];
]
][[2,1]]

]

