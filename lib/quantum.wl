(* ::Package:: *)

(*

    Conor Stevenson August 2019

    This mathematica package allows handling of quantum operators for 
    composite systems, primarily many-body problems, for example chains 
    of spins in 1D. It provides basic functionality to construct the operators 
    and store them as sparse arrays. It also provies routines to build 
    Hamiltonians for 1D chains of spins, in particular the transverse-field 
    anisotropic XY model Hamiltonian.

    It additionally provides functions to build dissipators for such 
    spin chains, in particular Lindblad dissipators, and Redfield 
    dissipators. The Redfield dissipators are constructed for the cases
    of structured environments at T=0. Such dissipators are formulated
    in the language of eigenoperators; a full enumeration of the
    eigenstates of the spin chain is thus required. This code can do
    this enumeration efficiently in the case of no anisotropy (unperturbed
    dissipator) because the Hamiltonian can be projected into it's 
    definite-spin sub-subspaces using the total spin operator S_z.
    In the case of anisotropy, one needs to fully diagonalise the
    problem.

    The work in progress is a method to extend the sub-space projection
    idea in the case of *weak* anisotropy using first-order perturbation
    theory. This will build on the unperturbed dissipator code, but needs
    to be designed at this stage.  

*)


(* Basic spin operator basis *)
{Id, Sigmax, Sigmaz, Sigmay, Sigmap, Sigmam} = 
    SparseArray /@ { {{1, 0}, {0, 1}}, {{0, 1}, {1, 0}}, {{1, 0}, {0, -1}}, 
                     {{0, -I}, {I, 0}}, {{0, 1}, {0, 0}}, {{0, 0}, {1, 0}} };



Subscript[\[Sigma], x] = Sigmax; Subscript[\[Sigma], y] = Sigmay;
Subscript[\[Sigma], z] = Sigmaz; 
SuperPlus[\[Sigma]] = Sigmap;
SuperMinus[\[Sigma]] = Sigmam;

Subscript[\[ScriptCapitalI], n_Integer] := IdentityMatrix[{n, n}];
Lower[n_Integer] := (Lower[n] = Table[If[j == i + 1, Sqrt[i], 0], {i, n}, {j, n}] // SparseArray);
Raise[n_Integer] := (Raise[n] = Lower[n]\[Transpose]);

Fock[nn_, m_] := UnitVector[nn, m + 1] // SparseArray;
NullVector[nn_] := ConstantArray[0, {nn}] // SparseArray

Coherent[nn_, \[Alpha]_] := Exp[-Abs[\[Alpha]]^2/2] MatrixExp[\[Alpha] Raise[nn]] . Fock[nn, 0];


\[CapitalGamma]Lorentzian[\[Omega]0_,a_] := a(a/((\[Omega]0-#)^2+a^2) - I (\[Omega]0-#)/((\[Omega]0-#)^2+a^2))&;
\[CapitalGamma]IMPROVED[\[Omega]0_,a_] := (I a)/(\[Omega]0-#+I a)&;


(******************************************************************************

    Building composite quantum systems' operators (tensor products)

*******************************************************************************)

(* Simply an alias for Mathematica's KroneckerProduct function a la QuTiP *)
Tensor[oper__] := KroneckerProduct[oper];
CircleTimes = Tensor;

(* Vectorisation of an operator. Maps density matrices into Liouville space. *)
Vectorise[mat_] := Flatten[Transpose[mat]];

(* Inverse of vectorisation *)
toHilbert[vect_] := ArrayReshape[vect, {Sqrt[Length[vect]], Sqrt[Length[vect]]}]//Transpose;

(* This returns a 'super'-operator, which acts on vectorised operators. 
   ApB = toHilbert[ SuperLR[A, B].Vectorise[p] ] *)
SuperLR[A_, B_] := Tensor[B // Transpose, A];

(* For short-hand, the super operator result of taking commutator. *)
SuperComm[oper_] :=  Module[{dims = Dimensions[oper]}, 
   SuperLR[oper, IdentityMatrix[dims]] - SuperLR[IdentityMatrix[dims], oper]
];

Lindblad[oper_] :=  2 SuperLR[oper, oper\[HermitianConjugate]] - 
                    ( SuperLR[oper\[HermitianConjugate] . oper, IdentityMatrix[oper // Dimensions]] + 
                      SuperLR[IdentityMatrix[oper // Dimensions], oper\[HermitianConjugate] . oper] );

(******************************************************************************

    Doing physics with these operators: Expectation values, 
    dynamics (dissipation), steady-states and time-evolution

*******************************************************************************)

Expect[oper_, rho_] := Tr[oper . rho];
(*hc = HermitianConjugate;
Lindblad[oper_] := 2 SuperLR[oper, oper//hc] - 
                (SuperLR[oper//hc].oper, IdentityMatrix[oper // Dimensions]] + 
                 SuperLR[IdentityMatrix[oper // Dimensions], oper//hc.oper]); *)


(* Now, we have the code which will calculate *)

MESolve[rhoinit_, \[ScriptCapitalL]_, {tmin_, tmax_}] := Module[
    {
        eqns, init, soln,
        state = Array[Subscript[rho, ##][t] &, {Length[rhoinit], Length[rhoinit]}] // Vectorise
    },
    eqns = Thread[D[state, t] == \[ScriptCapitalL] . state];
    init = Thread[(state /. {t -> tmin}) == Vectorise[rhoinit]];
    soln = toHilbert[state] /. NDSolve[{eqns, init}, state, {t, tmin, tmax}][[1]]
];