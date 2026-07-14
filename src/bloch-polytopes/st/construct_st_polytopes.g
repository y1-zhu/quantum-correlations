#Construction of bundle from complex reflection groups
#
#
#Currently work only for dimension 3;
Print("*************************************************************************\n");
Print("*The relevant polytopes from complex reflection groups:\n");
#Load character table lib:
LoadPackage( "ctbllib" );
#Load data from CHEVIE:
Read("st.g");

#Generate the generators of the group:
#ST number funs from 4 to 35 (more require more dimension) Test upto dim 6
#Not yet implemented for dim 7 and 8 (computation maybe too much)

STnum:= 27;
Print("*STnum: ",STnum,"\n");

Gen:=ST(STnum); #Unbind(ST);

MatGen:= Gen[1];
PermGen:= Gen[2];
Unbind(Gen);

PGroup:= GroupByGenerators(PermGen);;
OrderPGroup:= Size(PGroup);
Print("*Order: ", OrderPGroup,"\n");

MGroup:= GroupByGenerators(MatGen);;
d:=DimensionOfMatrixGroup(MGroup);;
Print("*Dimension: ", d,"\n");
#Conjugate action:
ConjugateAct:= function(x,g)
 return g*x*Inverse(g);
end;
#

#Mapping from permutation to matrices and back
PMHom:= GroupHomomorphismByImagesNC(PGroup, MGroup, PermGen, MatGen);;
MPHom:= GroupHomomorphismByImagesNC(MGroup, PGroup, MatGen, PermGen);;

#The Maske(?) inner product matrix
H:=Sum(MGroup, g -> ComplexConjugate(TransposedMat(g))*g)/Order(MGroup);

#Take the conjugacy classes of subgroups:
allConjugacyClassesSubgroups:= ConjugacyClassesSubgroups(PGroup);
Print("*Number of conjugacy classes of subgroups: ", Length(allConjugacyClassesSubgroups),"\n");

GoodSubgroups:=[];
GeneratingProjections:=[];
GeneratingAtomProjections:=[];
for theSubgroupClass in allConjugacyClassesSubgroups do
  #Take the relevant subgroup
  theSubgroup:= Representative(theSubgroupClass);

  #Order of the subgroup
  theOrder:= Size(theSubgroup);;

  #Compute the orbit length
  theOrbitLength:= OrderPGroup/theOrder;;

  #Compute the characters of the subgroup
  theCharacter:= List(Elements(theSubgroup), w -> Trace(Image(PMHom,w)));;

  #Check if there are only two irreps:
  #c:= Sum(theCharacter, x -> ComplexConjugate(x)*x)/Length(theCharacter);;

  #Look up the character table:
  theCharTable:= CharacterTable(theSubgroup);;
  theClassSizes:= SizesConjugacyClasses(theCharTable);;
  theConjugacyClasses:=ConjugacyClasses(theCharTable);;
  #Filter the 1D rep:
  theOneDIrreps:= Filtered(Irr(theCharTable),x -> x[1]=1);;

  #Recompute the character as class function
  theCharacterClass:= List(theConjugacyClasses, cw -> Trace(Image(PMHom,Representative(cw))));

  #Check whether there is 1D irrep
  haveOneDIrrep:=false;
  theProjection:=IdentityMat(d);
  for OneDIrrep in theOneDIrreps do
	#Compute the character overlap:
  	s:= Sum([1..Length(theClassSizes)], k -> theClassSizes[k]*ComplexConjugate(theCharacterClass[k])*OneDIrrep[k])/Length(theCharacter);
    #If there is a 1D representation, compute the projection on the isotropic component
	if (s>0) then #So there is 1D OneDIrrep in there
	  #Compute the projection on the the 1D
	  thePartialSums:=[];
      #Sum over each conjugacy class:
      for theClass in theConjugacyClasses do
         Add(thePartialSums,Sum(List(theClass),g->Image(PMHom,Inverse(g))));
      od;
      #Sum over all conjugacy class:
      theProjection:=Sum([1..Length(OneDIrrep)], k -> OneDIrrep[k]*thePartialSums[k])*OneDIrrep[1]/Length(theCharacter);
      haveOneDIrrep:=true;
      break;
    fi;
  od;

  #Add to the list of groups:
  rs := RandomSource(IsMersenneTwister); 
  if haveOneDIrrep then
	#Print("*Start testing...\n");
	theVector:=List([1..d],x -> Random(rs,-100,100)/Random(rs,-100,100));
	theVector:=TransposedMat([theProjection*theVector]);
	theAtomProjection:=theVector*ComplexConjugate(TransposedMat(theVector))*H; 
	theMSubgroup:=Image(PMHom,theSubgroup);
	#Print(List(theMSubgroup, g -> g*theAtomProjection*Inverse(g)-theAtomProjection,"\n"));
	#Print(theAtomProjection-theVector,"\n");
	Add(GeneratingProjections,theProjection);
	Add(GeneratingAtomProjections,theAtomProjection);
    Add(GoodSubgroups,theSubgroup);
  fi;
od;

Print("*Number of candidate generating projections: ", Length(GoodSubgroups),"\n");

#
#
#Reclassify orbit by types (i.e., by the conjugacy classes of the stabilizer group)
#
#

SClassIds:= [];
for k in [1..Length(GoodSubgroups)] do
  theSubgroup:= GoodSubgroups[k];
  theProjection:= GeneratingAtomProjections[k];
  allProjections:= Orbit(MGroup,theProjection,ConjugateAct);
  numProjections:= Length(allProjections);
  if (allProjections[1] <> theProjection) then
   Print("***Error: generating projections is not the first projection!");
  fi;
  #Action Homomorphism:
  ComHom:= CompositionMapping2(ActionHomomorphism(MGroup,allProjections,ConjugateAct),PMHom);
  ComAct:= function(x,g)
    return OnPoints(x,Image(ComHom,g));
  end;
  SGroup:= Stabilizer(PGroup,1,ComAct);
  Add(SClassIds,PositionProperty(allConjugacyClassesSubgroups, Class -> (SGroup in Class)));
  Print("*k=",k," L=",numProjections," L0=",OrbitLength(MGroup,GeneratingProjections[k],ConjugateAct),"\n",SClassIds,"\n");
od;

filteredGeneratingProjections:=[];
for id in Set(SClassIds) do
  Add(filteredGeneratingProjections,GeneratingAtomProjections[PositionProperty(SClassIds, k -> (k=id))]);
od;

Print("*Number of seclected generating projections: ", Length(filteredGeneratingProjections),"\n");
#
#
#
#


#
#
#Now we start to construct the orbits
#
#

AppendComplexMatrixForJuliaTo:= function(NAME,M)
  L1:=Length(M);
  L2:=Length(M[1]);	
  for q1 in [1..L1] do
	for q2 in [1..L2] do
		AppendTo(NAME,RealPart(Float(M[q1][q2])),"+",ImaginaryPart(Float(M[q1][q2])), "im\t");
	od;
	AppendTo(NAME,"\n");
  od;	
end;

LoadPackage("float");
SetFloats(MPC);
FLOAT.VIEW_DIG:=12;

I:= IdentityMat(d);
#Iterate over the selected subgroups:
Print("*The actual orbits:\n");
for k in [1..Length(filteredGeneratingProjections)] do
  theProjection:=filteredGeneratingProjections[k];
  Print(">>Generating projection: \n",theProjection,"\n");
  allProjections:= Orbit(MGroup,theProjection,ConjugateAct);
  numProjections:= Length(allProjections);
  Print(">Number of all projections: ", numProjections,"\n");
  #Check if the first projection is the original generating projection:
  if (allProjections[1] <> theProjection) then
   Print("***Error: generating projections is not the first projection! \n");
  fi;
  FILENAME:=Concatenation("polytope_st",String(STnum),"_d",String(d),"_L",String(numProjections),"_u",String(k),".txt");
  PrintTo(FILENAME, STnum, " #STnumber\n");
  AppendTo(FILENAME,0," #=Nonorthogonal\n");
  AppendTo(FILENAME,d," #Dimension\n");
  AppendTo(FILENAME, numProjections, " #OrbitLength\n");
  AppendComplexMatrixForJuliaTo(FILENAME,H);
  AppendTo(FILENAME,"\n");
  for eachProjection in allProjections do
	AppendComplexMatrixForJuliaTo(FILENAME,eachProjection);	
    AppendTo(FILENAME,"\n");
  od;
od;
Print("*Subprogram exited!");
