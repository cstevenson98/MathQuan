MyFourier[timeSeries_]:=Block[
    {
        \[CapitalDelta]t=(timeSeries[[-1,1]]-timeSeries[[1,1]])/(Length[timeSeries]-1),
        \[Omega]List
    },
    \[Omega]List=Table[(2\[Pi](s-1))/(Length[timeSeries]\[CapitalDelta]t),{s,1,Length[timeSeries]}];

    { \[Omega]List,Fourier[timeSeries[[All,2]]] } \[Transpose]
];

MyFourierSymmetric[timeSeries_]:=Block[
    {
        \[CapitalDelta]t=(timeSeries[[2,1]]-timeSeries[[1,1]]),
        n=Floor[Length[timeSeries]/2],
        \[Omega]List,
        rotatedTimeSeries
    },
    rotatedTimeSeries=RotateLeft[timeSeries,n];
    \[Omega]List=Table[(2\[Pi] s)/(Length[timeSeries]\[CapitalDelta]t),{s,-n,n}];

    { \[Omega]List,RotateRight[Fourier[rotatedTimeSeries[[All,2]]],n] } \[Transpose]
]
