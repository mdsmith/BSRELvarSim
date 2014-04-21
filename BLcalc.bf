LoadFunctionLibrary("GrabBag");
LoadFunctionLibrary("dSdNTreeTools");
LoadFunctionLibrary("CF3x4");
LoadFunctionLibrary("BranchSiteTemplateMod.mdl");
LoadFunctionLibrary("TreeTools");

SetDialogPrompt ("Get branch lengths from this fit:");
ExecuteAFile    (PROMPT_FOR_FILE);

LoadFunctionLibrary("queryTree");
bNames = BranchName (givenTree, -1);

ChoiceList  (type,"What is the most complicated model in this tree?",1,NO_SKIP,"MG94", "one syn and nonsyn parameter per branch", "BSREL","t and one or more omega parameters", "BSRELsyn", "one or more syn and non syn parameters");

if (type < 0)
{
    return 1;
}

bls = {};

if (type == 0)
{
    bls = BranchLength (givenTree,-1);
}
if (type == 2)
{
    //bls = extractBranchLengthsFromBSREL_SRV ("mixtureTree");
    bls = extractBranchLengthsFromBSREL_SRV ("bsrel_tree");
}
if (type == 1)
{
    //bls = extractBranchLengthsFromBSREL ("mixtureTree");
    bls = extractBranchLengthsFromBSREL ("bsrel_tree");
}
if (type > 2)
{
    fprintf(stdout, "type not found: ");
    fprintf(stdout, type);
    fprintf(stdout, "\n");
}

fprintf(stdout, "bls: \n");
fprintf(stdout, bls);
fprintf(stdout, "\n");

newickString = Format(givenTree,1,0);
fprintf(stdout, "\n");
fprintf(stdout, newickString);
fprintf(stdout, "\n");
for (bI = 0; bI < Abs(bls); bI = bI + 1)
{
    newickString = newickString^{{bNames[bI] + ",", bNames[bI] + ":" + bls[bNames[bI]] + ","}};
    newickString = newickString^{{bNames[bI] + ")", bNames[bI] + ":" + bls[bNames[bI]] + ")"}};
}
fprintf(stdout, "\n");
fprintf(stdout, newickString);
fprintf(stdout, "\n");





// define aux models

/*
Model M1 = (MGMatrix1, codon3x4, 0);
Model M2 = (MGMatrix2, codon3x4, 0);
Model M3 = (MGMatrix3, codon3x4, 0);

GetString (M1L, M1, -1);
GetString (M2L, M2, -1);
GetString (M3L, M3, -1);

fprintf (stdout, "\nName,Length");

bnames = BranchName (mixtureTree, -1);

for (bid = 0; bid < Columns (bnames)-1; bid += 1)
{
    fprintf (stdout, "\n", bnames[bid], ",", computeBLByName ("mixtureTree", bnames[bid]));
}

fprintf (stdout, "\n");

//----------------------------------------------------------------------------------------

function computeBLByName (treeName, branchName)
{
    t = Eval ("`treeName`.`branchName`.t");
    t1 = Eval ("`treeName`.`branchName`.t1");
    t2 = Eval ("`treeName`.`branchName`.t2");
    t3 = Eval ("`treeName`.`branchName`.t3");
    omega1 = Eval ("`treeName`.`branchName`.omega1");
    omega2 = Eval ("`treeName`.`branchName`.omega2");
    omega3 = Eval ("`treeName`.`branchName`.omega3");
    Paux1 = Eval ("`treeName`.`branchName`.Paux1");
    Paux2 = Eval ("`treeName`.`branchName`.Paux2");

    return Eval ("(Paux1*(`M1L`)+(1-Paux1)*Paux2*(`M2L`)+(1-Paux1)*(1-Paux2)*(`M3L`))/3");
}
*/
