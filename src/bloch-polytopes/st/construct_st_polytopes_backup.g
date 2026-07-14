Read("st.g");
gens:=ST(27);
G:= Group(gens[1]);

H:=Sum(G, g -> ComplexConjugate(TransposedMat(g))*g)/Order(G); 

#Only for diagonal H:
S:= List(H, r -> List(r,x -> Sqrt(x)));
iS:= Inverse(S);


#Conjugate Action:
ConjugateAct:= function(x,g)
	return S*g*iS*x*S*Inverse(g)*iS; 
end;

#Pick up a projection:
V:=[[1],[1],[1]];
M0:=V*TransposedMat(V);

OB:= Orbit(G,M0,ConjugateAct);

i:=E(4);

GellMann:=[];

GellMann[1]:= [[1, 0, 0],
              [0, 1, 0],
              [0, 0, 1]];
GellMann[1]:=GellMann[1]*Sqrt(2/3);

GellMann[2]:= [[0, 1, 0],
              [1, 0, 0],
              [0, 0, 0]];
GellMann[3]:= [[0, (-i), 0],
              [(+i), 0,  0],
              [0,  0,  0]];

GellMann[4]:= [[1, 0, 0],
              [0, (-1), 0],
              [0, 0, 0]];

GellMann[5]:= [[0, 0, 1],
              [0, 0, 0],
              [1, 0, 0]];

GellMann[6]:= [[0, 0, (-i)],
              [0, 0, 0],
              [(+i), 0, 0]];

GellMann[7]:= [[0, 0, 0],
              [0, 0, 1],
              [0, 1, 0]];

GellMann[8]:= [[0, 0, 0],
              [0, 0, (-i)],
              [0, (+i), 0]];

GellMann[9]:= [[1, 0, 0],
              [0, 1, 0],
              [0, 0, (-2)]];
GellMann[9]:=GellMann[9]*Sqrt(1/3);

ToBlochVector:= function(M)
	return List(GellMann, x -> Trace(x*M));
end;

L:= Length(OB);
VOB:= List(OB, M-> ToBlochVector(M));
FVOB:= Float(RealPart(VOB));

FILENAME:=Concatenation("polytope_st27_L",String(L),".txt");

FLOAT.VIEW_DIG=12
for V in FVOB do
	for x in V do
		AppendTo(FILENAME,x);
		AppendTo(FILENAME,"\t");
	od;
	AppendTo(FILENAME,"\n");
od;


