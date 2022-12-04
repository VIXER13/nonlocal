(* Content-type: application/vnd.wolfram.mathematica *)

(*** Wolfram Notebook File ***)
(* http://www.wolfram.com/nb *)

(* CreatedBy='Mathematica 12.1' *)

(*CacheID: 234*)
(* Internal cache information:
NotebookFileLineBreakTest
NotebookFileLineBreakTest
NotebookDataPosition[       158,          7]
NotebookDataLength[     14120,        319]
NotebookOptionsPosition[     12850,        291]
NotebookOutlinePosition[     13340,        309]
CellTagsIndexPosition[     13297,        306]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{
Cell[BoxData[
 RowBox[{
  RowBox[{"SetDirectory", "[", 
   RowBox[{"NotebookDirectory", "[", "]"}], "]"}], " ", ";"}]], "Input",
 CellChangeTimes->{{3.8670625888441153`*^9, 3.86706260972563*^9}, {
  3.8670634857645283`*^9, 3.867063493536092*^9}, {3.8676513129588475`*^9, 
  3.8676513197106085`*^9}, {3.867651476677971*^9, 3.867651477766203*^9}, {
  3.867658910707138*^9, 3.8676589108683133`*^9}},
 CellLabel->"In[11]:=",ExpressionUUID->"54842709-66fe-49e0-897c-527e9e597c76"],

Cell[CellGroupData[{

Cell["README", "Subsection",
 CellChangeTimes->{{3.8737362369545727`*^9, 3.873736263996457*^9}, {
  3.879166723320084*^9, 
  3.8791667398465815`*^9}},ExpressionUUID->"a261da3e-841c-4a19-b9ce-\
ce5578667ea7"],

Cell[TextData[{
 "\:0414\:043e\:0431\:0440\:044b\:0439 \:0434\:0435\:043d\:044c! \:0414\:0430\
\:043d\:043d\:044b\:0439 \:0444\:0430\:0439\:043b Wolfram mathematica \:043f\
\:0440\:0435\:0434\:043d\:0430\:0437\:043d\:0430\:0447\:0435\:043d \:0434\
\:043b\:044f \:0432\:0438\:0437\:0443\:0430\:043b\:0438\:0437\:0430\:0446\
\:0438\:0438 \:0440\:0430\:0441\:0447\:0435\:0442\:043e\:0432 \:043f\:0430\
\:043a\:0435\:0442\:0430 \[OpenCurlyDoubleQuote]Thermal nonstationary 1d\
\[CloseCurlyDoubleQuote]. \n\:041d\:0430\:0437\:0432\:0430\:043d\:0438\:044f \
\:043f\:0430\:043f\:043e\:043a \:0441 \:0440\:0430\:0441\:0447\:0435\:0442\
\:0430\:043c\:0438 \:043d\:0430 \:0440\:0430\:0437\:043d\:044b\:0445 \:0432\
\:0440\:0435\:043c\:0435\:043d\:043d\:044b\:0445 \:0441\:043b\:043e\:044f\
\:0445 \:0438\:043c\:0435\:044e\:0442 \:0441\:043b\:0435\:0434\:0443\:044e\
\:0449\:0435\:044e \:0441\:0442\:0440\:0443\:043a\:0442\:0443\:0440\:0443:\n\n\
\[OpenCurlyDoubleQuote]test/<\:041d\:0430\:0437\:0432\:0430\:043d\:0438\:0435 \
\:0440\:0430\:0441\:0447\:0435\:0442\:0430>/Flux<\:043d\:043e\:043c\:0435\
\:0440 \:0432\:0440\:0435\:043c\:0435\:043d\:043d\:043e\:0433\:043e \:0441\
\:043b\:043e\:044f>.csv\[CloseCurlyDoubleQuote] \:0438\:043b\:0438 \
\[OpenCurlyDoubleQuote]test/<\:041d\:0430\:0437\:0432\:0430\:043d\:0438\:0435 \
\:0440\:0430\:0441\:0447\:0435\:0442\:0430>/T<\:043d\:043e\:043c\:0435\:0440 \
\:0432\:0440\:0435\:043c\:0435\:043d\:043d\:043e\:0433\:043e \:0441\:043b\
\:043e\:044f>.csv\[CloseCurlyDoubleQuote] ",
 StyleBox["\n", "Subsubsection"],
 "\n\:0412\:043d\:0443\:0442\:0440\:0438 \:043a\:0430\:0436\:0434\:043e\:0439 \
\:043f\:0430\:043f\:043a\:0438 \:043d\:0430\:0445\:043e\:0434\:044f\:0442\
\:0441\:044f \:0444\:0430\:0439\:043b\:044b \:0444\:043e\:0440\:043c\:0430\
\:0442\:0430 ",
 StyleBox["csv", "Subsubsection"],
 ", \:0432 \:043a\:0430\:0436\:0434\:043e\:043c \:0438\:0437 \:043a\:043e\
\:0442\:043e\:0440\:044b\:0445 \:043d\:0430\:0445\:043e\:0434\:0438\:0442\
\:0441\:044f \:043e\:0434\:0438\:043d \:0432\:0440\:0435\:043c\:0435\:043d\
\:043d\:043e\:0439 \:0441\:043b\:043e\:0439. \:041d\:0443\:043c\:0435\:0440\
\:0430\:0446\:0438\:044f \:043d\:0430\:0447\:0438\:043d\:0430\:0435\:0442\
\:0441\:044f \:0441 0.\n\:0414\:043b\:044f \:043f\:0440\:043e\:0441\:043c\
\:043e\:0442\:0440\:0430 \:0444\:0430\:0439\:043b\:043e\:0432 \:0432 \:043f\
\:0430\:043f\:043a\:0435 \:043c\:043e\:0436\:043d\:043e \:0438\:0441\:043f\
\:043e\:043b\:044c\:0437\:043e\:0432\:0430\:0442\:044c \:0444\:0443\:043d\
\:043a\:0446\:0438\:044e ",
 Cell[BoxData[
  RowBox[{"FileNames", "[", "]"}]], "Input",ExpressionUUID->
  "20e15c15-531b-4e15-8c29-5b1a46bf69fc"],
 "."
}], "Text",
 CellChangeTimes->{{3.873735464670876*^9, 3.873735779896126*^9}, {
  3.8737358514652853`*^9, 3.8737359942725997`*^9}, {3.8737360248390827`*^9, 
  3.873736222003023*^9}, {3.873736912137215*^9, 3.8737369391896467`*^9}, {
  3.8739017320039096`*^9, 3.873901747143067*^9}, {3.873901818152877*^9, 
  3.873901828544733*^9}, {3.8767262634424887`*^9, 3.8767262652084546`*^9}, {
  3.8791664319751444`*^9, 
  3.879166522735956*^9}},ExpressionUUID->"85b6a899-07a3-46b2-971b-\
925e5764243e"],

Cell[BoxData[{
 RowBox[{
  RowBox[{
   RowBox[{"WriteData", "[", 
    RowBox[{"foldername_", ",", " ", "timelayers_"}], "]"}], ":=", 
   RowBox[{"Module", "[", 
    RowBox[{
     RowBox[{"{", "data", "}"}], ",", "\[IndentingNewLine]", 
     RowBox[{
      RowBox[{"SetDirectory", "[", 
       RowBox[{
        RowBox[{"NotebookDirectory", "[", "]"}], "<>", "foldername"}], "]"}], 
      " ", ";", "\[IndentingNewLine]", 
      RowBox[{"(*", 
       RowBox[{
        RowBox[{"timelayers", "=", 
         RowBox[{"Max", "@", 
          RowBox[{"(", 
           RowBox[{
            RowBox[{
             RowBox[{"(", 
              RowBox[{"ToExpression", "[", 
               RowBox[{
                RowBox[{"StringSplit", "[", 
                 RowBox[{"#", ",", "\"\<.\>\""}], "]"}], 
                "\[LeftDoubleBracket]", "1", "\[RightDoubleBracket]"}], "]"}],
               ")"}], "&"}], "/@", 
            RowBox[{"FileNames", "[", "]"}]}], ")"}]}]}], ";", " ", 
        RowBox[{
        "\:041a\:043e\:043b\:0438\:0447\:0435\:0441\:0442\:0432\:043e", " ", 
         "\:0432\:0440\:0435\:043c\:0435\:043d\:043d\:044b\:0445", " ", 
         "\:0441\:043b\:043e\:0435\:0432"}]}], "*)"}], "\[IndentingNewLine]", 
      RowBox[{"dataFlux", " ", "=", 
       RowBox[{"Table", "[", 
        RowBox[{
         RowBox[{"Import", "[", 
          RowBox[{"\"\<Flux\>\"", "<>", 
           RowBox[{"ToString", "[", "i", "]"}], "<>", "\"\<.csv\>\""}], "]"}],
          " ", ",", 
         RowBox[{"{", 
          RowBox[{"i", ",", "0", ",", "timelayers", ",", "1"}], "}"}]}], 
        "]"}]}], ";", "\[IndentingNewLine]", 
      RowBox[{"dataT", " ", "=", 
       RowBox[{"Table", "[", 
        RowBox[{
         RowBox[{"Import", "[", 
          RowBox[{"\"\<T\>\"", "<>", 
           RowBox[{"ToString", "[", "i", "]"}], "<>", "\"\<.csv\>\""}], "]"}],
          " ", ",", 
         RowBox[{"{", 
          RowBox[{"i", ",", "0", ",", "timelayers", ",", "1"}], "}"}]}], 
        "]"}]}], ";", "\[IndentingNewLine]", 
      RowBox[{"SetDirectory", "[", 
       RowBox[{"NotebookDirectory", "[", "]"}], "]"}], " ", ";", 
      "\[IndentingNewLine]", 
      RowBox[{"Return", "[", 
       RowBox[{
        RowBox[{"{", "dataT", "}"}], "~", "Join", "~", 
        RowBox[{"{", "dataFlux", "}"}]}], "]"}], ";"}]}], 
    "\[IndentingNewLine]", "]"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"Visualize", "[", 
    RowBox[{
    "datalist_List", ",", "name_List", ",", " ", "start_", ",", " ", 
     "finish_"}], "]"}], ":=", 
   RowBox[{"Module", "[", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{
       RowBox[{"max", " ", "=", " ", 
        RowBox[{"Max", "[", 
         RowBox[{
          RowBox[{
           RowBox[{"(", 
            RowBox[{"#", "[", 
             RowBox[{"[", "2", "]"}], "]"}], ")"}], "&"}], "/@", 
          RowBox[{"Flatten", "[", 
           RowBox[{"datalist", ",", "2"}], "]"}]}], "]"}]}], ",", " ", 
       RowBox[{"min", " ", "=", " ", 
        RowBox[{"Min", "[", 
         RowBox[{
          RowBox[{
           RowBox[{"(", 
            RowBox[{"#", "[", 
             RowBox[{"[", "2", "]"}], "]"}], ")"}], "&"}], "/@", 
          RowBox[{"Flatten", "[", 
           RowBox[{"datalist", ",", "2"}], "]"}]}], "]"}]}], ",", " ", 
       RowBox[{"timelayers", " ", "=", " ", 
        RowBox[{"Min", "@", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{
            RowBox[{"(", 
             RowBox[{"Length", "@", "#"}], ")"}], "&"}], "/@", "datalist"}], 
          ")"}]}]}]}], "}"}], ",", "\[IndentingNewLine]", 
     RowBox[{
      RowBox[{"Print", "[", 
       RowBox[{
       "\"\<\:041a\:043e\:043b\:0438\:0447\:0435\:0441\:0442\:0432\:043e \
\:0432\:0440\:0435\:043c\:0435\:043d\:043d\:044b\:0445 \:0441\:043b\:043e\
\:0435\:0432: \>\"", "<>", 
        RowBox[{"ToString", "[", "timelayers", "]"}]}], "]"}], ";", 
      "\[IndentingNewLine]", 
      RowBox[{"Return", "[", 
       RowBox[{"Animate", "[", 
        RowBox[{
         RowBox[{"ListPlot", "[", 
          RowBox[{
           RowBox[{"Table", "[", 
            RowBox[{
             RowBox[{
              RowBox[{"datalist", "[", 
               RowBox[{"[", "j", "]"}], "]"}], "[", 
              RowBox[{"[", "i", "]"}], "]"}], ",", 
             RowBox[{"{", 
              RowBox[{"j", ",", "1", ",", 
               RowBox[{"Length", "@", "datalist"}]}], "}"}]}], "]"}], ",", 
           " ", 
           RowBox[{"PlotRange", "\[Rule]", 
            RowBox[{"{", 
             RowBox[{
              RowBox[{"{", 
               RowBox[{"start", ",", "finish"}], "}"}], ",", 
              RowBox[{"{", 
               RowBox[{"min", ",", "max"}], "}"}]}], "}"}]}], ",", 
           RowBox[{"PlotLegends", "\[Rule]", "name"}]}], "]"}], ",", 
         RowBox[{"{", 
          RowBox[{"i", ",", "1", ",", "timelayers", ",", "1"}], "}"}], ",", 
         RowBox[{"AnimationRunning", "\[Rule]", "False"}]}], "]"}], "]"}], 
      ";"}]}], "\[IndentingNewLine]", "]"}]}], ";"}]}], "Input",
 CellChangeTimes->{{3.8737376059524307`*^9, 3.8737377607858796`*^9}, {
  3.8737378043931446`*^9, 3.873737809385235*^9}, {3.873737894813581*^9, 
  3.873738120823843*^9}, {3.8737382057615194`*^9, 3.8737382103852444`*^9}, {
  3.8737382538579197`*^9, 3.873738287665844*^9}, {3.873738396737567*^9, 
  3.8737385261358547`*^9}, {3.8737385664122877`*^9, 3.873738576075731*^9}, {
  3.873771002402439*^9, 3.8737710129487143`*^9}, {3.8737710452277904`*^9, 
  3.873771055326606*^9}, {3.8738849577203517`*^9, 3.873884986844984*^9}, {
  3.873885061180679*^9, 3.873885063476758*^9}, {3.873902312355627*^9, 
  3.873902323068385*^9}, {3.8739025087525373`*^9, 3.873902574121496*^9}, {
  3.873902645624254*^9, 3.873902664709199*^9}, {3.8791665888518987`*^9, 
  3.879166624155263*^9}, {3.8791670587007713`*^9, 3.879167083232066*^9}},
 CellLabel->"In[12]:=",ExpressionUUID->"f571884d-41eb-4b7e-9af6-a63c00e2e2c0"]
}, Open  ]],

Cell[CellGroupData[{

Cell["\:041f\:0440\:0438\:043c\:0435\:0440 \:0441\:0447\:0438\:0442\:044b\
\:0432\:0430\:043d\:0438\:044f \:0444\:0430\:0439\:043b\:0430", "Subsection",
 CellChangeTimes->{{3.867651322624304*^9, 
  3.867651330370652*^9}},ExpressionUUID->"f591a71d-3929-4b3c-80f8-\
e27e942b3ac0"],

Cell["\:0414\:0435\:043c\:043e\:043d\:0441\:0442\:0440\:0430\:0446\:0438\:044f\
 \:0432\:043e\:0437\:043c\:043e\:0436\:043d\:043e\:0439 \:0432\:0438\:0437\
\:0443\:0430\:043b\:0438\:0437\:0430\:0446\:0438\:0438.", "Text",
 CellChangeTimes->{{3.873736305980817*^9, 
  3.873736336439495*^9}},ExpressionUUID->"63fe5c5b-1766-40e4-a4ca-\
e067d9ae79b3"],

Cell[BoxData[
 RowBox[{
  RowBox[{"datafortest", "=", 
   RowBox[{"WriteData", "[", 
    RowBox[{"\"\<tests\\left_flux_1_right_flux_-1\>\"", ",", " ", "100"}], 
    "]"}]}], ";"}]], "Input",
 CellChangeTimes->{{3.867651357818698*^9, 3.867651394272673*^9}, {
   3.873735326874516*^9, 3.87373533019947*^9}, {3.8737362898335705`*^9, 
   3.8737363033505836`*^9}, {3.8737363528852615`*^9, 3.873736493461133*^9}, {
   3.8737369453704586`*^9, 3.8737369486561575`*^9}, {3.8737373877848487`*^9, 
   3.8737374429448643`*^9}, 3.873737482825537*^9, {3.873737774931605*^9, 
   3.873737792924161*^9}, {3.873901771026435*^9, 3.873901787200121*^9}, {
   3.873901867610112*^9, 3.873901927877803*^9}, 3.8739019837599134`*^9, {
   3.8739020816906643`*^9, 3.873902091145079*^9}, {3.879166704381258*^9, 
   3.8791667166419067`*^9}, {3.8791667553785763`*^9, 3.879166801092543*^9}, {
   3.879167096617634*^9, 3.879167097539531*^9}, {3.8791671354142904`*^9, 
   3.8791671527609854`*^9}},
 EmphasizeSyntaxErrors->True,
 CellLabel->"In[20]:=",ExpressionUUID->"667ae529-4007-4a9f-b834-98b644bc3645"],

Cell[BoxData[
 RowBox[{"Length", "@", 
  RowBox[{"datafortest", "[", 
   RowBox[{"[", "2", "]"}], "]"}]}]], "Input",
 CellChangeTimes->{{3.8739024659340396`*^9, 3.873902471395752*^9}, {
  3.879167108245916*^9, 3.8791671114509563`*^9}},
 CellLabel->"In[17]:=",ExpressionUUID->"30343047-a86b-43a2-873d-3fee26e0d49e"],

Cell[BoxData[
 RowBox[{"ListPlot", "[", 
  RowBox[{"Table", "[", 
   RowBox[{
    RowBox[{"datafortest", "[", 
     RowBox[{"[", 
      RowBox[{"1", ",", "j"}], "]"}], "]"}], ",", 
    RowBox[{"{", 
     RowBox[{"j", ",", "1", ",", "100", ",", "5"}], "}"}]}], "]"}], 
  "]"}]], "Input",
 CellChangeTimes->{{3.8791674037602415`*^9, 3.8791674650587645`*^9}},
 CellLabel->"In[27]:=",ExpressionUUID->"84f39767-5854-4b2a-8565-ae2ca112f484"]
}, Open  ]]
},
WindowSize->{947.6307692307693, 465.7846153846154},
WindowMargins->{{
  23.26153846153846, Automatic}, {-3.876923076923049, Automatic}},
Magnification:>0.8 Inherited,
FrontEndVersion->"12.1 for Microsoft Windows (64-bit) (June 19, 2020)",
StyleDefinitions->"Default.nb",
ExpressionUUID->"ee5b16b9-5aaa-434e-b992-e61bbbe507c5"
]
(* End of Notebook Content *)

(* Internal cache information *)
(*CellTagsOutline
CellTagsIndex->{}
*)
(*CellTagsIndex
CellTagsIndex->{}
*)
(*NotebookFileOutline
Notebook[{
Cell[558, 20, 475, 8, 35, "Input",ExpressionUUID->"54842709-66fe-49e0-897c-527e9e597c76"],
Cell[CellGroupData[{
Cell[1058, 32, 207, 4, 44, "Subsection",ExpressionUUID->"a261da3e-841c-4a19-b9ce-ce5578667ea7"],
Cell[1268, 38, 3154, 49, 147, "Text",ExpressionUUID->"85b6a899-07a3-46b2-971b-925e5764243e"],
Cell[4425, 89, 5914, 142, 356, "Input",ExpressionUUID->"f571884d-41eb-4b7e-9af6-a63c00e2e2c0"]
}, Open  ]],
Cell[CellGroupData[{
Cell[10376, 236, 278, 4, 44, "Subsection",ExpressionUUID->"f591a71d-3929-4b3c-80f8-e27e942b3ac0"],
Cell[10657, 242, 347, 5, 29, "Text",ExpressionUUID->"63fe5c5b-1766-40e4-a4ca-e067d9ae79b3"],
Cell[11007, 249, 1072, 18, 23, "Input",ExpressionUUID->"667ae529-4007-4a9f-b834-98b644bc3645"],
Cell[12082, 269, 314, 6, 35, "Input",ExpressionUUID->"30343047-a86b-43a2-873d-3fee26e0d49e"],
Cell[12399, 277, 435, 11, 35, "Input",ExpressionUUID->"84f39767-5854-4b2a-8565-ae2ca112f484"]
}, Open  ]]
}
]
*)

