(* ::Package:: *)

(* ::Section:: *)
(*Begin Package*)


BeginPackage["QuantileNucleiSegmentation`"];


initialSegmentation::usage = "initialSegmentation[filename] takes in a path to an image stack and performs a crude segmentation";
mergeNeighbours::usage = "mergeNeighbours[segmentation] takes in the segmented imagestack and merge smaller blobs (false nuclei) with
their neighbours";
QuantileNuclei::usage = "QuantileNuclei[segmentation,ImageData,filename] takes in an image stack, its filename and the segmented stack
to bound nuclei using quantiles"
primaryStats::usage = "primaryStats[segmentation] takes in a segmented stack and outputs the basic stats";
segmentStack::usage = "segmentStack[filename] takes in a filename and starts the segmentation procedure";


Begin["`Private`"];


(* ::Subsection:: *)
(*Initial Segmentation*)


Options[initialSegmentation]:={"fillholes"->20,"smallcomponents"-> 100,"maxdetectthresh"-> 0.01,"gfilterthresh"-> 2,
"printStatus"-> False};
initialSegmentation[filename_?StringQ,OptionsPattern[]]:=Module[{stack,img3D,imgData,imgDim,roi,img3DBin,imgBinData,
imgDataDim,height, fillholes,distance,markers, seg, areas, background,centroids,fillh = OptionValue["fillholes"],
smallcomp=OptionValue["smallcomponents"]},

stack = Import@filename;
img3D = Image3D[stack];
If[OptionValue["printStatus"],Print@img3D];
imgDim = ImageDimensions@img3D;
roi = {ConstantArray[1,3],Reverse@imgDim}\[Transpose];
img3D = ImageTake[img3D,##]&@@roi;

img3DBin =Binarize[img3D];
imgBinData = ImageData[img3DBin]; (* data from binarized 3d image *)
imgData = ImageData[img3D]; (* data from non-binarized 3d image *)
imgDataDim =Dimensions[imgData];
height = First@imgDataDim;
(* now let us fill small-holes and gaps *)
fillholes[imagedata_,size_]:= Block[{img = imagedata, filledImage},
filledImage = (FillingTransform@Image[img])~DeleteSmallComponents~size;
ImageData[filledImage]
];
 Do[imgData[[i]] = imgData[[i]]~fillholes ~fillh;
imgBinData[[i]] = imgBinData[[i]]~fillholes~fillh,{i,height}
]; 
img3D = Image3D[imgData]//DeleteSmallComponents[#,smallcomp]&;
img3DBin = Image3D[imgBinData]//DeleteSmallComponents[#,smallcomp]&;
(* after filling holes we need to find seeds and segment the image using watershed *)
distance = ImageAdjust@DistanceTransform[img3DBin,Padding-> 0] ;(* distance transform of the image *)
markers = MaxDetect[distance,OptionValue["maxdetectthresh"]]; (* markers for segmentation *)
seg = WatershedComponents[GradientFilter[img3D,OptionValue["gfilterthresh"]],markers,Method->"Rainfall"];
(* watershed on non-binarized 3D image *)
(* removing background *)
areas = ComponentMeasurements[seg,"Area"]; 
background = MaximalBy[areas,Last][[1,1]];  
seg = ArrayComponents[seg,Length@areas,background-> 0];
centroids = ComponentMeasurements[seg,"Centroid"];
If[OptionValue["printStatus"],Print@Colorize[seg]];
seg
];


(* ::Subsection::Closed:: *)
(*merge small nuclei with neighbours*)


Options[mergeNeighbours]:={"offset"-> 10000,"maxdistance"-> 7.0,"areaval"-> 15000};
mergeNeighbours[seg_?ArrayQ,OptionsPattern[]]:=Module[{areasM, ncM,centroidsM,areaval=OptionValue@"areaval",
neighboursM,assoc,mergecandidates,closecentroids,nearest,nearestpairList,nearestneighboursList,mergeNeighboursHelper,
maxdistance=OptionValue@"maxdistance",offset = OptionValue["offset"],position,segT=seg},

{areasM,ncM,centroidsM,neighboursM} =Map[ComponentMeasurements[seg,#]&,{"Area","NeighborCount","Centroid","Neighbors"}];

mergeNeighboursHelper[lis_,pat:{_,_}]:= Module[{indices},
indices=#1&@@@Position[nearestpairList,pat[[#]]]&/@{1,2};
If[Length@indices>1,
Join[lis,{Union@Flatten@Extract[nearestpairList,Partition[indices,1]]}],
Join[lis,{pat}]]
];
assoc=Merge[{Association@areasM,Association@ncM},List@@#&]; 
mergecandidates = Keys@Select[assoc,(#[[1]]<areaval&&#[[2]]>0&)];(* cells with at least one neighbour and a size smaller than some
value *)
(* the part below finds the nearestneighbours of the individual cells in mergecandidates and sows the closest neighbour with the cell.
 any duplicate entries such as {cell1,cell2} and {cell2,cell1} are deleted. *)
nearestpairList=
Reap@If[mergecandidates!= {},
Table[closecentroids=Reverse/@Extract[centroidsM,Partition[neighboursM[[i]][[2]],1]];
nearest=Nearest[closecentroids,Last@@centroidsM[[i]],{1,maxdistance}];
If[nearest!={},Sow[{i,First@nearest}]],{i,mergecandidates}];
];
(* if one cell is also the nearest neighbour of another cell-pair then merge all the cells together *)
If[Last@nearestpairList != {},
(nearestpairList=DeleteDuplicates@*Map[Sort]@*First@*Last@nearestpairList;
nearestneighboursList= SortBy[First]@*Union@Fold[mergeNeighboursHelper,{},nearestpairList]; 
(* here we will reindex: merge the neighbouring cells together *)
Do[position=Position[segT,Alternatives@@i];
segT=ReplacePart[segT,Thread[position-> offset]];
offset=offset+1,{i,nearestneighboursList}])
];
(* now indexing them back *)
ArrayComponents[segT]
];


(* ::Subsection::Closed:: *)
(*primary stats*)


primaryStats[seg_]:= Module[{areasM,ncM,centroidsM,neighboursM,avgcellsize,areaval},
{areasM,ncM,centroidsM,neighboursM}= Map[ComponentMeasurements[seg,#]&,{"Area","NeighborCount","Centroid","Neighbors"}];
areaval = areasM[[All,2]];
Print["number of cells: ", Length@areasM];
Print[Style["average cell area in pixels: ",{Bold,Red,FontSize-> 14}],Style[avgcellsize= Mean@areaval,{Bold,FontSize->15}]];
Print@SmoothHistogram[areaval/avgcellsize,PlotStyle->{Thick,XYZColor[0,0,1,0.8]}];
Print@ListPlot[areaval/avgcellsize,PlotStyle->{XYZColor[0,0,1,0.6],PointSize[Large]}];
];


(* ::Subsection::Closed:: *)
(*ComponentMeasures/ImageProperties*)


imageStackData[file_,roi_]:= Module[{bc,bc3d,img,imgData},
bc=Import@file;
bc3d = Image3D[bc];
img = ImageTake[bc3d,Sequence@@roi];
{ImageData[img],ImageDimensions[img]}
];


componentMeasures[file_,seg_]:= Module[{imgdata,roi,masksM,boxesM,iW,iD,iH},
roi = Thread[{1,Dimensions@seg}];
{imgdata,{iW,iD,iH}} = imageStackData[file,roi];
(* these represent the masks for and bounding box surrounding each cell *)
masksM = ComponentMeasurements[seg,"Mask"]; (* yields a sparse array with dimensions same as the imagestack *)
boxesM = ComponentMeasurements[seg,"BoundingBox"];
{masksM,boxesM,imgdata,iW,iD,iH}
];


(* ::Subsection::Closed:: *)
(*Quantile Envelope*)


ClearAll[QuantileEnvelopeRegion];
QuantileEnvelopeRegion[points_?MatrixQ, quantile_?NumberQ, 
numberofDirections_?IntegerQ]:=Block[{nd = numberofDirections, dirs, rmats, qDirPoints,qRegion},
dirs = N@Flatten[
Table[{Cos[\[Theta]]Cos[\[Phi]], Sin[\[Theta]]Cos[\[Phi]],Sin[\[Phi]]},{\[Theta], 2 \[Pi]/(10 nd), 2 \[Pi],2 \[Pi]/nd},
{\[Phi],-\[Pi],\[Pi],2 \[Pi]/nd}],1];
rmats = RotationMatrix[{{1,0,0},#}]&/@dirs;
qDirPoints =
Flatten[Map[Function[{m},Quantile[(m.Transpose[points])[[3]],quantile]], rmats]];
qRegion = ImplicitRegion[MapThread[(#1.{x,y,z})[[3]] <= #2 &, {rmats, qDirPoints}],{x,y,z}];
qRegion
]/;Dimensions[points][[2]] == 3 && 0 < quantile <= 1;


(* ::Subsection:: *)
(*Nuclei Refinement*)


Options[QuantileNuclei]:={"quantile"->0.90,"directions"->20};
QuantileNuclei[segM_,file_,OptionsPattern[]]:= Module[{masksM,boxesM,iW,iH,iD,numcells,imgdata},

nucleiRefinement[seg_,id_]:=With[{cellID = id,directionOpt = OptionValue["directions"],quantileOpt = OptionValue["quantile"]},
Block[{b,singlecellData,singlecellImage,singlecellImageBin,pts,qRegion,qSurface,outlierpts,inlyingpts,cellpos,
segtemp,limits,unitcell},
Print["processing: ", cellID];
singlecellData=imgdata*masksM[[cellID,2]];
b = boxesM[[cellID,2]];
singlecellImage = ImageTake[Image3D@singlecellData,{iH-b[[2,3]],iH-b[[1,3]]+1},{iD-b[[2,2]],iD-b[[1,2]]+1},{b[[1,1]],b[[2,1]]+1}];

(* here we develop a quantile envelope for the binarized cell imagedata *)
singlecellImageBin = Binarize@singlecellImage;
pts= ImageValuePositions[singlecellImageBin,1];
qRegion=QuantileEnvelopeRegion[pts,quantileOpt,directionOpt];
qSurface = BoundaryDiscretizeRegion@qRegion;

(* make pixels outside qsurface zero *)
outlierpts=Select[pts,#\[NotElement]qRegion&];
singlecellImageBin=ReplaceImageValue[singlecellImageBin,outlierpts-> 0];

(* now we need to replace the nucleus in the segmentation structure: zero whatever is initially present and then
replace with inlyingpts *)
cellpos = Position[seg,cellID];
segtemp=ReplacePart[seg,cellpos-> 0];

(* lets pad the image to put cell in the right place *)
limits = {{b[[1,1]]-1,iW-b[[2,1]]-1},{b[[1,2]]-1,iD-b[[2,2]]-1},{b[[1,3]]-1,iH-b[[2,3]]}};
unitcell=MorphologicalComponents@ImagePad[singlecellImageBin,limits];
inlyingpts=Position[unitcell,1];
(* replace inlying pts with cell index *)
segtemp=ReplacePart[segtemp,inlyingpts-> cellID]]
];

{masksM,boxesM,imgdata,iW,iD,iH} = componentMeasures[file,segM];
numcells=Length@masksM;
Fold[nucleiRefinement[#1,#2]&,segM,Range@numcells]
];


(* ::Subsection::Closed:: *)
(*for visualization only*)


visualizeNuclei[file_,seg_,imgdata_,iH_,iD_,iW_]:=Module[{array,cellpos,pixelvals,cellnum},
cellnum = Length@ComponentMeasurements[seg,"Label"];
array=ConstantArray[0,{iH,iD,iW}]; 
(* initially empty to add nuclei *)
Do[
Print["processing: ", i];
cellpos=Position[seg,i];
pixelvals = Extract[imgdata,cellpos];
array=ReplacePart[array,Thread[cellpos->pixelvals]]
,{i,Range@cellnum}
];
Save[DirectoryName@file<>"nucleivisualization.res",array];
array
];


(* ::Subsection:: *)
(*Mains*)


Options[segmentStack]:= {"printseg"-> True,"printStats"->True,"mergeorder"-> 3,
"visualizeNuclei"-> True};
segmentStack[file_,OptionsPattern[]]:= Module[{seg,segM,imgdata,masksM,boxesM,iW,iH,iD,segNew},
seg=initialSegmentation[file];
segM=Nest[mergeNeighbours,seg,OptionValue@"mergeorder"];

If[OptionValue["printseg"]==True,Print@*Colorize@segM];
If[OptionValue["printStats"]==True,primaryStats[segM]];
Save[DirectoryName[file]<>"firstseg.res",segM];

segNew = QuantileNuclei[segM,file];

Save[DirectoryName[file]<>"postQuantileseg.res",segNew];

If[OptionValue["printseg"]==True,Print@*Colorize@segNew];
If[OptionValue["printStats"]==True,primaryStats[segNew]];

If[OptionValue["visualizeNuclei"]==True,
Print@*Image3D@visualizeNuclei[file,segNew,imgdata,iH,iD,iW]
];
segNew
];


End[];


(* ::Section:: *)
(*End Package*)


EndPackage[];
