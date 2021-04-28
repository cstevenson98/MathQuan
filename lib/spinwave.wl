(* ::Package:: *)

Import["quantum.wl"];


(* ::Input::Initialization:: *)
CxxKilda::usage="CxxKilda[n, g, \[CapitalDelta], \[Kappa]] returns symbolic \!\(\*UnderscriptBox[\(\[Sum]\), \(k\)]\)\!\(\*SuperscriptBox[SubscriptBox[\(C\), \(k\)], \(xx\)]\)[\[Tau]] for an anisotropic XY-model with a flat bath J(\[Omega]) = \[Kappa].";
CxxKilda[n_,g_,\[CapitalDelta]_,\[Kappa]_]:=Module[
{
kList=N@Range[-\[Pi],\[Pi],2\[Pi]/(n-1)],
\[Epsilon]=2(g+Cos[#])&,
\[Eta]=2\[CapitalDelta] Cos[#]&,
\[Xi]
},
\[Xi]=Sqrt[\[Epsilon][#]^2-\[Eta][#]^2]&;

Exp[-\[Kappa] \[Tau]]Total[(Cos[\[Xi][#]\[Tau]]+I (\[Eta][#]-\[Epsilon][#]) Sin[\[Xi][#]\[Tau]]/\[Xi][#]+(\[Eta][#](\[Eta][#]-\[Epsilon][#]))/(\[Xi][#]^2+\[Kappa]^2) (Cos[\[Xi][#]\[Tau]]+\[Kappa] Sin[\[Xi][#]\[Tau]]/\[Xi][#]))&/@kList]
]


(* ::Input::Initialization:: *)
XYBogoliubovCompile=Compile[
{
{g,_Real},{\[CapitalDelta],_Real},
{\[Kappa],_Real},{\[Omega]0,_Real},{b,_Real},

{k,_Real},{#2,_Real}
},

Block[
{
\[Theta]k=ArcTanh[-\[CapitalDelta] Cos[k]/(g+Cos[k])],
\[Xi]k=2Sign[g+Cos[k]]Sqrt[(g+Cos[k])^2-\[CapitalDelta]^2 Cos[k]^2],

Jp,Sp,Jm,Sm,
Uk,Vk,
\[CapitalGamma],\[CapitalLambda],\[CapitalGamma]p,\[CapitalGamma]m,\[Mu]\[Mu],

Q,P,R,QR,PR,
x,y,xtilde,ytilde,
X,Y,\[Alpha],\[Beta],
Z,W,\[Gamma],\[Delta],

C1,C2,C3,\[Gamma]0,Qk,det
},

Uk=Cosh[\[Theta]k];Vk=Sinh[\[Theta]k];
(*Jp=(\[Kappa] b^2)/(b^2+(\[Omega]0-\[Xi]k)^2);
Sp=(\[Kappa] b (-\[Omega]0+\[Xi]k))/(b^2+(\[Omega]0-\[Xi]k)^2);Jm=(\[Kappa] b^2)/(b^2+(\[Omega]0+\[Xi]k)^2);Sm=(\[Kappa] b (-\[Omega]0-\[Xi]k))/(b^2+(\[Omega]0+\[Xi]k)^2);*)
Jp=(\[Kappa]*b^2)/((\[Omega]0-\[Xi]k)^2 + b^2);
Sp=(\[Kappa]*b*(\[Omega]0-\[Xi]k))/((\[Omega]0-\[Xi]k)^2 + b^2);
Jm=(\[Kappa]*b^2)/((\[Omega]0+\[Xi]k)^2 + b^2);
Sm=(\[Kappa]*b*(\[Omega]0+\[Xi]k))/((\[Omega]0+\[Xi]k)^2 + b^2);


\[CapitalGamma]p=Jp+I Sp;\[CapitalGamma]m=Jm+I Sm;
\[CapitalGamma]=(Uk^2 \[CapitalGamma]p-Vk^2 \[CapitalGamma]m\[Conjugate]);\[CapitalLambda]=Uk Vk(\[CapitalGamma]m\[Conjugate]-\[CapitalGamma]p);

\[Mu]\[Mu]=Sqrt[(Im[\[CapitalGamma]]+\[Xi]k)^2-Abs[\[CapitalLambda]]^2];

(* Spinwave-Redfield relevant terms *)

\[Gamma]0=-2Uk^2Re[\[CapitalGamma]p]+2Vk^2Re[\[CapitalGamma]m];
Qk=-2I \[Xi]k-2\[CapitalGamma];
det=-\[Gamma]0 Abs[Qk]^2-4Re[Qk]Abs[\[CapitalLambda]]^2;
C1=-1/det (2Vk^2Abs[Qk]^2Re[\[CapitalGamma]m]-2Uk Vk(\[CapitalLambda] Qk\[Conjugate] \[CapitalGamma]m+\[CapitalLambda]\[Conjugate]Qk \[CapitalGamma]m\[Conjugate]));
C2=-1/det (4 Vk^2 \[CapitalLambda]\[Conjugate] Qk\[Conjugate] Re[\[CapitalGamma]m]+2 Uk Vk(\[Gamma]0 Qk\[Conjugate] \[CapitalGamma]m+2Abs[\[CapitalLambda]]^2\[CapitalGamma]m-2(\[CapitalLambda]^2)\[Conjugate]\[CapitalGamma]m\[Conjugate]));
C3=-1/det (4 Vk^2 \[CapitalLambda] Qk Re[\[CapitalGamma]m]+2 Uk Vk(\[Gamma]0 Qk \[CapitalGamma]m\[Conjugate]+2Abs[\[CapitalLambda]]^2\[CapitalGamma]m\[Conjugate]-2\[CapitalLambda]^2\[CapitalGamma]m));

(* Langevin relevant terms *)

Q=1/2 (1-(\[Xi]k+Im[\[CapitalGamma]])/\[Mu]\[Mu]);P=1/2 (1+(\[Xi]k+Im[\[CapitalGamma]])/\[Mu]\[Mu]);R=I/2 \[CapitalLambda]/\[Mu]\[Mu];
QR=(Q+R);PR=(P-R);

x=(Jp Uk^2 QR PR\[Conjugate]+Jm Vk^2 QR PR\[Conjugate]-Uk Vk(Jp (PR\[Conjugate])^2+Jm QR^2));
y=(Jp Uk^2 PR PR\[Conjugate]+Jm Vk^2 QR QR\[Conjugate]-Uk Vk(Jp QR\[Conjugate]PR\[Conjugate]+Jm QR PR));

xtilde=(Jp Uk^2 PR QR\[Conjugate]+Jm Vk^2 PR QR\[Conjugate]-Uk Vk(Jp (QR\[Conjugate])^2+Jm PR^2));
ytilde=(Jp Uk^2 QR QR\[Conjugate]+Jm Vk^2 PR PR\[Conjugate]-Uk Vk(Jp QR\[Conjugate]PR\[Conjugate]+Jm QR PR));

X=2 QR PR\[Conjugate](Jp Uk^2+Jm Vk^2)-Uk Vk(QR^2+(PR\[Conjugate])^2)(Jp+Jm);
Y=(Abs[QR]^2+Abs[PR]^2)(Jp Uk^2+Jm Vk^2)-2Uk Vk Re[QR PR](Jp+Jm);

Z=Uk Vk(Jp-Jm)(QR^2-(PR\[Conjugate])^2);
W=(Jp Uk^2-Jm Vk^2)(Abs[PR]^2-Abs[QR]^2)+2I Uk Vk (Jp-Jm)Im[QR PR];

\[Gamma]=((Re[\[CapitalGamma]] Im[Z]-\[Mu]\[Mu] Re[Z])/(Re[\[CapitalGamma]]^2+\[Mu]\[Mu]^2)+Im[W]/Re[\[CapitalGamma]]);\[Delta]=((Re[\[CapitalGamma]] Re[Z]+\[Mu]\[Mu] Im[Z])/(Re[\[CapitalGamma]]^2+\[Mu]\[Mu]^2)+Re[W]/Re[\[CapitalGamma]]);
\[Alpha]=(Y/Re[\[CapitalGamma]]+(Re[\[CapitalGamma]] Re[X]-\[Mu]\[Mu] Im[X])/(Re[\[CapitalGamma]]^2+\[Mu]\[Mu]^2));\[Beta]=-((Re[\[CapitalGamma]] Im[X]+\[Mu]\[Mu] Re[X])/(Re[\[CapitalGamma]]^2+\[Mu]\[Mu]^2));
#1
],

CompilationTarget->"C",
RuntimeOptions->"EvaluateSymbolically"->False
]&;


(* ::Input::Initialization:: *)
(*\[InvisibleComma] \[LeftAngleBracket]\[Sigma]^x\[Sigma]^x\[RightAngleBracket](\[Omega]): *)
CxxLorentzianOmegaVsk::usage="CxxLorentzianOmegaVsk[g, \[CapitalDelta], \[Kappa], \[Omega]0, b, k, \[Omega]] returns numerical \!\(\*SuperscriptBox[SubscriptBox[\(C\), \(k\)], \(xx\)]\)[\[Omega]] for an anisotropic XY-model with a Lorentzian structured bath J(\[Omega]) = \!\(\*FractionBox[\(\[Kappa]\\\ \*SuperscriptBox[\(b\), \(2\)]\), \(\*SuperscriptBox[\(b\), \(2\)] + \*SuperscriptBox[\((\[Omega]0 - \[Omega])\), \(2\)]\)]\).";

CxxLorentzianOmegaVsk=XYBogoliubovCompile[
\[Pi] (Uk+Vk)^2 ((x/(Re[\[CapitalGamma]]-I \[Mu]\[Mu])+y/Re[\[CapitalGamma]]) (2Re[\[CapitalGamma]])/(Re[\[CapitalGamma]]^2+(\[Omega]+\[Mu]\[Mu])^2)
+(xtilde/(Re[\[CapitalGamma]]+I \[Mu]\[Mu])+ytilde/Re[\[CapitalGamma]]) (2Re[\[CapitalGamma]])/(Re[\[CapitalGamma]]^2+(\[Omega]-\[Mu]\[Mu])^2)),\[Omega]];

CxxLorentzianOmegaVskIMPROVED=XYBogoliubovCompile[(Uk+Vk)^2 /4((-QR/(-Re[\[CapitalGamma]]+I(\[Mu]\[Mu]-\[Omega]))+-PR/(-Re[\[CapitalGamma]]-I(\[Mu]\[Mu]+\[Omega]))+QR/(Re[\[CapitalGamma]]+I(\[Mu]\[Mu]-\[Omega]))+PR/(Re[\[CapitalGamma]]-I(\[Mu]\[Mu]+\[Omega])))(C1+C2+1)+(-QR\[Conjugate]/(-Re[\[CapitalGamma]]-I(\[Mu]\[Mu]+\[Omega]))+-PR\[Conjugate]/(-Re[\[CapitalGamma]]+I(\[Mu]\[Mu]-\[Omega]))+QR\[Conjugate]/(Re[\[CapitalGamma]]-I(\[Mu]\[Mu]+\[Omega]))+PR\[Conjugate]/(Re[\[CapitalGamma]]+I(\[Mu]\[Mu]-\[Omega])))(C1+C3)),\[Omega]];

CxxLorentzianOmegaIMPROVED[g_,\[CapitalDelta]_,\[Kappa]_,\[Omega]0_,b_,\[Omega]_]:=1/(2\[Pi])NIntegrate[CxxLorentzianOmegaVskIMPROVED[g,\[CapitalDelta],\[Kappa],\[Omega]0,b,k,\[Omega]],{k,-\[Pi],\[Pi]},Method->{Automatic,"SymbolicProcessing"->0}];

SxxLorentzianOmegaVskIMPROVED=XYBogoliubovCompile[(Uk+Vk)^2 /8((-QR/(-Re[\[CapitalGamma]]+I(\[Mu]\[Mu]-\[Omega]))+-PR/(-Re[\[CapitalGamma]]-I(\[Mu]\[Mu]+\[Omega]))+QR/(Re[\[CapitalGamma]]+I(\[Mu]\[Mu]-\[Omega]))+PR/(Re[\[CapitalGamma]]-I(\[Mu]\[Mu]+\[Omega])))(C1+C2+(C1+C3)\[Conjugate]+1)+(-QR\[Conjugate]/(-Re[\[CapitalGamma]]-I(\[Mu]\[Mu]+\[Omega]))+-PR\[Conjugate]/(-Re[\[CapitalGamma]]+I(\[Mu]\[Mu]-\[Omega]))+QR\[Conjugate]/(Re[\[CapitalGamma]]-I(\[Mu]\[Mu]+\[Omega]))+PR\[Conjugate]/(Re[\[CapitalGamma]]+I(\[Mu]\[Mu]-\[Omega])))(C1+C3+(C1+C2)\[Conjugate]+1)),\[Omega]];

SxxLorentzianOmegaIMPROVED[g_,\[CapitalDelta]_,\[Kappa]_,\[Omega]0_,b_,\[Omega]_]:=1/(2\[Pi])NIntegrate[SxxLorentzianOmegaVskIMPROVED[g,\[CapitalDelta],\[Kappa],\[Omega]0,b,k,\[Omega]],{k,-\[Pi],\[Pi]},Method->{Automatic,"SymbolicProcessing"->0}];

ChixxLorentzianOmegaVskIMPROVED=XYBogoliubovCompile[-I(Uk+Vk)^2/4 ((-QR/(-Re[\[CapitalGamma]]+I(\[Mu]\[Mu]-\[Omega]))+-PR/(-Re[\[CapitalGamma]]-I(\[Mu]\[Mu]+\[Omega])))(C1+C2-(C1+C3)\[Conjugate]+1)+(-QR\[Conjugate]/(-Re[\[CapitalGamma]]-I(\[Mu]\[Mu]+\[Omega]))+-PR\[Conjugate]/(-Re[\[CapitalGamma]]+I(\[Mu]\[Mu]-\[Omega])))(C1+C3-(C1+C2)\[Conjugate]-1)),\[Omega]];

ChixxLorentzianOmegaIMPROVED[g_,\[CapitalDelta]_,\[Kappa]_,\[Omega]0_,b_,\[Omega]_]:=1/(2\[Pi])NIntegrate[ChixxLorentzianOmegaVskIMPROVED[g,\[CapitalDelta],\[Kappa],\[Omega]0,b,k,\[Omega]],{k,-\[Pi],\[Pi]},Method->{Automatic,"SymbolicProcessing"->0}];

SxxSpatialLorentzianOmegaVskIMPROVED=XYBogoliubovCompile[1/4 (Exp[I m k](Uk^2 C1+Vk^2 (C1+1)+Uk Vk(C2+C3))+Exp[-I m k](Uk^2 C2+Vk^2 C3+Uk Vk(2C1+1))+Exp[-I m k](Uk^2 C1\[Conjugate]+Vk^2 (C1\[Conjugate]+1)+Uk Vk(C2\[Conjugate]+C3\[Conjugate]))+Exp[I m k](Uk^2 C2\[Conjugate]+Vk^2 C3\[Conjugate]+Uk Vk(2C1\[Conjugate]+1))),m];

SxxSpatialLorentzianOmegaIMPROVED[g_,\[CapitalDelta]_,\[Kappa]_,\[Omega]0_,b_,\[Omega]_]:=1/(2\[Pi])NIntegrate[SxxSpatialLorentzianOmegaVskIMPROVED[g,\[CapitalDelta],\[Kappa],\[Omega]0,b,k,\[Omega]],{k,-\[Pi],\[Pi]},Method->{Automatic,"SymbolicProcessing"->0}]



CxxLorentzianOmega::usage="CxxLorentzianOmega[g, \[CapitalDelta], \[Kappa], \[Omega]0, b, \[Omega]] returns numerical \!\(\*FractionBox[\(1\), \(2  \[Pi]\)]\)\!\(\*SubsuperscriptBox[\(\[Integral]\), \(-\[Pi]\), \(\[Pi]\)]\)\!\(\*SuperscriptBox[SubscriptBox[\(C\), \(k\)], \(xx\)]\)[\[Omega]]\[DifferentialD]k for an anisotropic XY-model with a Lorentzian structured bath J(\[Omega]) = \!\(\*FractionBox[\(\[Kappa]\\\ \*SuperscriptBox[\(b\), \(2\)]\), \(\*SuperscriptBox[\(b\), \(2\)] + \*SuperscriptBox[\((\[Omega]0 - \[Omega])\), \(2\)]\)]\).";

CxxLorentzianOmega[g_,\[CapitalDelta]_,\[Kappa]_,\[Omega]0_,b_,\[Omega]_]:=1/(2\[Pi]) NIntegrate[CxxLorentzianOmegaVsk[g,\[CapitalDelta],\[Kappa],\[Omega]0,b,k,\[Omega]],{k,-\[Pi],\[Pi]},Method->{Automatic,"SymbolicProcessing"->0}];

CxxLorentzianTimeVsk::usage="CxxLorentzianOmegaVsk[g, \[CapitalDelta], \[Kappa], \[Omega]0, b, k, \[Tau]] returns numerical \!\(\*SuperscriptBox[SubscriptBox[\(C\), \(k\)], \(xx\)]\)[\[Tau]] for an anisotropic XY-model with a Lorentzian structured bath J(\[Omega]) = \!\(\*FractionBox[\(\[Kappa]\\\ \*SuperscriptBox[\(b\), \(2\)]\), \(\*SuperscriptBox[\(b\), \(2\)] + \*SuperscriptBox[\((\[Omega]0 - \[Omega])\), \(2\)]\)]\).";

CxxLorentzianTimeVsk=XYBogoliubovCompile[\[Pi](Exp[-Re[\[CapitalGamma]]Abs[\[Tau]]+I \[Mu]\[Mu] \[Tau]](x/(Re[\[CapitalGamma]]-I \[Mu]\[Mu])+y/Re[\[CapitalGamma]])+Exp[-Re[\[CapitalGamma]] Abs[\[Tau]]-I \[Mu]\[Mu] \[Tau]](ytilde/Re[\[CapitalGamma]] +xtilde/(Re[\[CapitalGamma]]+I \[Mu]\[Mu]))),\[Tau]];
CxxLorentzianTime::usage="CxxLorentzianOmega[g, \[CapitalDelta], \[Kappa], \[Omega]0, b, \[Omega]] returns numerical \!\(\*FractionBox[\(1\), \(2  \[Pi]\)]\)\!\(\*SubsuperscriptBox[\(\[Integral]\), \(-\[Pi]\), \(\[Pi]\)]\)\!\(\*SuperscriptBox[SubscriptBox[\(C\), \(k\)], \(xx\)]\)[\[Tau]]\[DifferentialD]k for an anisotropic XY-model with a Lorentzian structured bath J(\[Tau]) = \!\(\*FractionBox[\(\[Kappa]\\\ \*SuperscriptBox[\(b\), \(2\)]\), \(\*SuperscriptBox[\(b\), \(2\)] + \*SuperscriptBox[\((\[Omega]0 - \[Omega])\), \(2\)]\)]\).";

CxxLorentzianTime[g_,\[CapitalDelta]_,\[Kappa]_,\[Omega]0_,b_,\[Omega]_]:=1/(2\[Pi]) NIntegrate[CxxLorentzianTimeVsk[g,\[CapitalDelta],\[Kappa],\[Omega]0,b,k,\[Omega]],{k,-\[Pi],\[Pi]},Method->{Automatic,"SymbolicProcessing"->0}];

SxxLorentzianOmegaVsk::usage="SxxLorentzianOmegaVsk[g, \[CapitalDelta], \[Kappa], \[Omega]0, b, k, \[Omega]] returns numerical \!\(\*SuperscriptBox[SubscriptBox[\(S\), \(k\)], \(xx\)]\)[\[Omega]] for an anisotropic XY-model with a Lorentzian structured bath J(\[Omega]) = \!\(\*FractionBox[\(\[Kappa]\\\ \*SuperscriptBox[\(b\), \(2\)]\), \(\*SuperscriptBox[\(b\), \(2\)] + \*SuperscriptBox[\((\[Omega]0 - \[Omega])\), \(2\)]\)]\).";

SxxLorentzianOmegaVsk=XYBogoliubovCompile[\[Pi] (Uk+Vk)^2 ((Re[\[CapitalGamma]](\[Alpha]-I \[Beta]))/(Re[\[CapitalGamma]]^2+(\[Mu]\[Mu]+\[Omega])^2)+(Re[\[CapitalGamma]](\[Alpha]+I \[Beta]))/(Re[\[CapitalGamma]]^2+(\[Mu]\[Mu]-\[Omega])^2)),\[Omega]];

SxxLorentzianOmega::usage="SxxLorentzianOmega[g, \[CapitalDelta], \[Kappa], \[Omega]0, b, \[Omega]] returns numerical \!\(\*FractionBox[\(1\), \(2  \[Pi]\)]\)\!\(\*SubsuperscriptBox[\(\[Integral]\), \(-\[Pi]\), \(\[Pi]\)]\)\!\(\*SuperscriptBox[SubscriptBox[\(S\), \(k\)], \(xx\)]\)[\[Omega]]\[DifferentialD]k for an anisotropic XY-model with a Lorentzian structured bath J(\[Omega]) = \!\(\*FractionBox[\(\[Kappa]\\\ \*SuperscriptBox[\(b\), \(2\)]\), \(\*SuperscriptBox[\(b\), \(2\)] + \*SuperscriptBox[\((\[Omega]0 - \[Omega])\), \(2\)]\)]\).";

SxxLorentzianOmega[g_,\[CapitalDelta]_,\[Kappa]_,\[Omega]0_,b_,\[Omega]_,Strategy_:Automatic]:=1/(2\[Pi]) NIntegrate[SxxLorentzianOmegaVsk[g,\[CapitalDelta],\[Kappa],\[Omega]0,b,k,\[Omega]],{k,-\[Pi],\[Pi]},Method->{Strategy,"SymbolicProcessing"->0}];

ChixxLorentzianOmegaVsk::usage="ChixxLorentzianOmegaVsk[g, \[CapitalDelta], \[Kappa], \[Omega]0, b, k, \[Omega]] returns numerical \!\(\*SuperscriptBox[SubscriptBox[\(\[Chi]\), \(k\)], \(xx\)]\)[\[Omega]] for an anisotropic XY-model with a Lorentzian structured bath J(\[Omega]) = \!\(\*FractionBox[\(\[Kappa]\\\ \*SuperscriptBox[\(b\), \(2\)]\), \(\*SuperscriptBox[\(b\), \(2\)] + \*SuperscriptBox[\((\[Omega]0 - \[Omega])\), \(2\)]\)]\).";

ChixxLorentzianOmegaVsk=XYBogoliubovCompile[\[Pi] (Uk+Vk)^2 ((\[Gamma]-I \[Delta])/(Re[\[CapitalGamma]]-I(\[Omega]+\[Mu]\[Mu]))+(\[Gamma]+I \[Delta])/(Re[\[CapitalGamma]]-I(\[Omega]-\[Mu]\[Mu]))),\[Omega]];

ChixxLorentzianOmega::usage="ChixxLorentzianOmega[g, \[CapitalDelta], \[Kappa], \[Omega]0, b, \[Omega]] returns numerical \!\(\*FractionBox[\(1\), \(2  \[Pi]\)]\)\!\(\*SubsuperscriptBox[\(\[Integral]\), \(-\[Pi]\), \(\[Pi]\)]\)\!\(\*SuperscriptBox[SubscriptBox[\(S\), \(k\)], \(xx\)]\)[\[Omega]]\[DifferentialD]k for an anisotropic XY-model with a Lorentzian structured bath J(\[Omega]) = \!\(\*FractionBox[\(\[Kappa]\\\ \*SuperscriptBox[\(b\), \(2\)]\), \(\*SuperscriptBox[\(b\), \(2\)] + \*SuperscriptBox[\((\[Omega]0 - \[Omega])\), \(2\)]\)]\).";

ChixxLorentzianOmega[g_,\[CapitalDelta]_,\[Kappa]_,\[Omega]0_,b_,\[Omega]_,Strategy_:Automatic]:=1/(2\[Pi]) NIntegrate[ChixxLorentzianOmegaVsk[g,\[CapitalDelta],\[Kappa],\[Omega]0,b,k,\[Omega]],{k,-\[Pi],\[Pi]},Method->{Strategy,"SymbolicProcessing"->0}];


(* ::Input::Initialization:: *)
TeffNLorentzian[n_,g_,\[CapitalDelta]_,\[Kappa]_,\[Omega]0_,b_,\[Omega]max_,NSamp_]:=
Block[
{
FinverseSamps=Table[{\[Omega],Re@Total@Table[SxxLorentzianOmegaVsk[g,\[CapitalDelta],\[Kappa],\[Omega]0,b,k,\[Omega]],{k,N@Range[-\[Pi],\[Pi],2\[Pi]/(n-1)]}]/Im@Total@Table[ChixxLorentzianOmegaVsk[g,\[CapitalDelta],\[Kappa],\[Omega]0,b,k,\[Omega]],{k,N@Range[-\[Pi],\[Pi],2\[Pi]/(n-1)]}]},{\[Omega],0.001,\[Omega]max,\[Omega]max/NSamp}]
},
myFitRule=FindFit[FinverseSamps,A Coth[B \[Omega]],{A,B},\[Omega]];
Return[A/(2B)/.myFitRule];
]


(* ::Input::Initialization:: *)
TeffInfiniteLorentzian[g_,\[CapitalDelta]_,\[Kappa]_,\[Omega]0_,b_,\[Omega]max_,NSamp_,Strategy_:Automatic]:=
Block[
{
FinverseSamps=Table[{\[Omega],Re@SxxLorentzianOmega[g,\[CapitalDelta],\[Kappa],\[Omega]0,b,\[Omega],Strategy]/Im@ChixxLorentzianOmega[g,\[CapitalDelta],\[Kappa],\[Omega]0,b,\[Omega],Strategy]},{\[Omega],0.1,\[Omega]max,\[Omega]max/NSamp}],
myFitRule
},
myFitRule=FindFit[FinverseSamps,A Coth[B \[Omega]],{A,B},\[Omega]];
Return[A/(2B)/.myFitRule];
];

TeffInfiniteLorentzianIMPROVED[g_,\[CapitalDelta]_,\[Kappa]_,\[Omega]0_,b_,\[Omega]max_,NSamp_]:=
Block[
{
FinverseSamps=Table[{\[Omega],Re@SxxLorentzianOmegaIMPROVED[g,\[CapitalDelta],\[Kappa],\[Omega]0,b,\[Omega]]/Im@ChixxLorentzianOmegaIMPROVED[g,\[CapitalDelta],\[Kappa],\[Omega]0,b,\[Omega]]},{\[Omega],0.1,\[Omega]max,\[Omega]max/NSamp}],
myFitRule
},
myFitRule=FindFit[FinverseSamps,A Coth[B \[Omega]],{A,B},\[Omega]];
Return[A/(2B)/.myFitRule];
]
