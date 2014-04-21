// XXX Omega input is treated as syn rate, not an actual omega. Results
// analysis will have to take this into account.
// test: /home/martin/Software/Simulation/alphaSimTest1/smallOut1.0
doSynRateVariation = 1;

SetDialogPrompt      ("Simulation settings:");
ExecuteAFile         (PROMPT_FOR_FILE);
//Tree bsrel_tree = (((1,2)5,(3,4)6)7);

//ExecuteAFile ("sim_settings.ibf");


//LoadFunctionLibrary("OmegaGenerator.bf");
//bsrel_settings = generate4TaxonSettings(2);

LoadFunctionLibrary ("chooseGeneticCode", {"0": "Universal"});
LoadFunctionLibrary ("BranchSiteTemplate");
LoadFunctionLibrary ("CF3x4");
LoadFunctionLibary  ("GrabBag");


nucCF						= CF3x4	(nuc3x4, GeneticCodeExclusions);
ModelMatrixDimension        = CountSenseCodons (_Genetic_Code);
codon3x4					= BuildCodonFrequencies (nucCF);

// get tree branch names

topology_branch_names = BranchName (bsrel_tree,-1);

matrices_defined          = 1;  
models_defined            = {};
branch_length_expressions = {};
solve_these_for_lengths   = {};

// ----- apply nucleoitide biases

nucleotide_bias_settings ["apply"][""];
function apply (key, value) {
    ^key = value;
    return 0;
}

// ----- DONE apply nucleoitide biases

for (branch_id = 0; branch_id < Columns(topology_branch_names)-1; branch_id += 1) {
    branch_name = topology_branch_names[branch_id];
    full_name = "bsrel_tree.`branch_name`";
    assert (Abs(bsrel_settings[branch_name]), "Missing simulation information for branch " + branch_name);
    assert (Abs((bsrel_settings[branch_name])["length"]) > 0.0, "Not a non-zero length for branch " + branch_name);
    
    omega_info = (bsrel_settings[branch_name])["omegas"];
    
    assert (Type (omega_info) == "Matrix", "Missing matrix with the omega distribution for " + branch_name);
    assert (Columns (omega_info) == 2, "Not a 2-column matrix with the omega distribution for " + branch_name);
    
    rate_classes   = Rows (omega_info);
    scaled_weights = omega_info[-1][1];
    scaled_weights = scaled_weights * (1/(+scaled_weights));
    
    if (models_defined [rate_classes] < 1) {
        for (; matrices_defined <= rate_classes; matrices_defined += 1) {
            //PopulateModelMatrix			  ("MGMatrix" + matrices_defined,  nucCF, "t", "omega" + matrices_defined, "");
            // XXX make the Modle matrix fit
            PopulateModelMatrix			  ("MGMatrix" + matrices_defined,  nucCF, "syn" + matrices_defined, "", "nonsyn" + matrices_defined);
            ExecuteCommands ("Model		MGlocal			= (MGMatrix" + matrices_defined +", codon3x4, 0)");
            GetString (bl_expression, MGlocal, -1);
            branch_length_expressions [matrices_defined] = bl_expression;
       }
        models_defined [rate_classes] = 1;
        
        if (rate_classes > 1) {        
            // XXX This produces essentially convolved branch length formulae
            // for the single length parameter. We need separate ones.
            freq_weights = generateFrequencyParameterization ("Paux", rate_classes);
            //expression = "";
            for (k = 0; k < rate_classes; k += 1) {
                expression = "";
                // XXX
                /*
                if (k) {
                    expression += "+";
                }
                */
                // XXX
                //expression += "(" + branch_length_expressions [k+1] + ")*" + freq_weights[k];
                expression += "(" + branch_length_expressions [k+1] + ")";
                // XXX?
                freq_weights[k] =  "Exp(MGMatrix" + (k+1) + ")*" + freq_weights[k];
            }
            ExecuteCommands ("Model 		MG" + rate_classes + "=(\"" + Join("+",freq_weights) + "\",codon3x4,EXPLICIT_FORM_MATRIX_EXPONENTIAL);");  
            solve_these_for_lengths [rate_classes] = expression; 
        } else {
            Model		MG1				= (MGMatrix1, codon3x4, 0);
            solve_these_for_lengths [rate_classes] = branch_length_expressions[1];
        }
    }
    
    
    ExecuteCommands ("SetParameter (`full_name`, MODEL, MG" + rate_classes +");");

    fprintf(stdout, "length formulae:\n", solve_these_for_lengths, "\n"); // XXX
    left_over_weight = 1;
    for (rate_class = 1; rate_class <= rate_classes; rate_class += 1) {
        // XXX no longer one omega, plus convert omega to nonsyn
        //*(full_name + ".omega" + rate_class) = omega_info[rate_class-1][0];
        //*("omega" + rate_class) = omega_info[rate_class-1][0];
        *(full_name + ".nonsyn" + rate_class) = omega_info[rate_class-1][0];
        *("nonsyn" + rate_class) = omega_info[rate_class-1][0];
        if (rate_class < rate_classes) {
            *(full_name + ".Paux" + rate_class) = scaled_weights[rate_class-1]/left_over_weight;
            *("Paux" + rate_class) = scaled_weights[rate_class-1]/left_over_weight;
            left_over_weight = left_over_weight - scaled_weights[rate_class-1];
        }

        // XXX Move into the rate class loop above
        desired_length = (bsrel_settings[branch_name])["length"];
        // XXX change length term
        //ExecuteCommands ("FindRoot (length_term, " + "((" + solve_these_for_lengths[rate_classes] + ")/3) - desired_length, t, 0, 10000)"); // XXX
        /*
        fprintf(stdout, "\n");
        ExecuteCommands ("fprintf(stdout, nonsyn1)"); // XXX
        fprintf(stdout, "\n");
        */
        fprintf(stdout, "FindRoot (length_term, " + "((" + solve_these_for_lengths[rate_class] + ")/3) - (desired_length * omega_info[rate_class-1][1], syn" + rate_class + ", 0, 10000)"); // XXX
        fprintf(stdout, "\n");
        fprintf(stdout, "DL total: " + desired_length); // XXX
        fprintf(stdout, "\n");
        fprintf(stdout, "\n");
        fprintf(stdout, "GT= " + GT); // XXX
        fprintf(stdout, "\n");
        fprintf(stdout, "CT= " + CT); // XXX
        fprintf(stdout, "\n");
        fprintf(stdout, "CG= " + CG); // XXX
        fprintf(stdout, "\n");
        fprintf(stdout, "nonsyn2 = " + nonsyn2); // XXX
        fprintf(stdout, "\n");
        fprintf(stdout, "AC= " + AC); // XXX
        fprintf(stdout, "\n");
        fprintf(stdout, "AT= " + AT); // XXX
        fprintf(stdout, "\n");
        ExecuteCommands ("fprintf(stdout, \"DL this class: \" + (desired_length * omega_info[rate_class-1][1]))"); // XXX
        fprintf(stdout, "\n");
        ExecuteCommands ("FindRoot (length_term, " + "((" + solve_these_for_lengths[rate_class] + ")/3) - (desired_length * omega_info[rate_class-1][1]), syn" + rate_class + ", 0, 10000)"); // XXX
        //ExecuteCommands("FindRoot (length_term, ((0.1191346538733348*syn1+0.2524983209976207*nonsyn1+0.0588460548502598*GT*syn1+0.2916039026420982*GT*nonsyn1+0.1354306490726096*CT*syn1+0.2150199128563157*CT*nonsyn1+0.05884610403133145*CG*syn1+0.2916065281183093*CG*nonsyn1+0.07176604092115288*AT*syn1+0.278683880588994*AT*nonsyn1+0.0876239628864896*AC*syn1+0.2628286332845085*AC*nonsyn1)/3) - desired_length, syn1, 0, 10000)");
        //ExecuteCommands("FindRoot (length_term, ((0.1191346538733348*syn1+0.2524983209976207*nonsyn1+0.0588460548502598*GT*syn1+0.2916039026420982*GT*nonsyn1+0.1354306490726096*CT*syn1+0.2150199128563157*CT*nonsyn1+0.05884610403133145*CG*syn1+0.2916065281183093*CG*nonsyn1+0.07176604092115288*AT*syn1+0.278683880588994*AT*nonsyn1+0.0876239628864896*AC*syn1+0.2628286332845085*AC*nonsyn1)/3) - desired_length, syn1, 0, 10000)");
        // XXX change length term
        //*(full_name + ".t") = length_term;
        fprintf(stdout, "\nsolved length: " + length_term + "\n\n");
        *(full_name + ".syn" + rate_class) = length_term;
    }
    printNodeDesc (full_name, rate_classes);
}

// ----- apply nucleoitide biases
/*

nucleotide_bias_settings ["apply"][""];
function apply (key, value) {
    ^key = value;
    return 0;
}
*/

// ----- DONE apply nucleoitide biases

codonCharacters = {{"A","C","G","T"}
			  			   {"3",GeneticCodeExclusions,"",""}};
			  			   
SetDialogPrompt ("Save simulated data to:");
fprintf 		(PROMPT_FOR_FILE, CLEAR_FILE);
save_sims_to     = LAST_FILE_PATH;
IS_TREE_PRESENT_IN_DATA = 1;
DATAFILE_TREE = Format (bsrel_tree,1,1);

fprintf (stdout, "\n");

for (it = 0; it < replicates; it += 1) {		  
    DataSet sim   = Simulate 	(bsrel_tree,codon3x4,codonCharacters,codons,0);
    DataSetFilter simFilter = CreateFilter (sim,1);	
	fName = save_sims_to + "." + it;
	fprintf 		(fName, CLEAR_FILE,simFilter);
    SetParameter (STATUS_BAR_STATUS_STRING, "Replicate " + (it+1) + " of " + replicates + " generated", 0);
    lfOut = save_sims_to + ".fit";
    DataSetFilter newFilter = CreateFilter(sim, 3, "", "", GeneticCodeExclusions);
    LikelihoodFunction sim_LF = (newFilter, bsrel_tree);
    LIKELIHOOD_FUNCTION_OUTPUT = 7;
    fprintf(lfOut, CLEAR_FILE, sim_LF);
    LIKELIHOOD_FUNCTION_OUTPUT = 2;
}		
//------------------------------------------------------------------------------------------------------------------------

lfunction generateFrequencyParameterization (prefix, count) {
    result = {1, count};
    if (count > 1) {
        result[0] = "`prefix`1";
        buffer = "(1-`prefix`1)";
        for (i = 1; i < count - 1; i+=1) {
            result[i] = buffer + "*" + "`prefix`" + (i+1);
            buffer   += "*(1-`prefix`" + (i+1) + ")";
        } 
        result[count-1] = buffer;
    } else {
        result [0] = "1";
    }
    return result;
}

//------------------------------------------------------------------------------------------------------------------------

lfunction printNodeDesc (ref, rate_classes) {
    do_srv = ^("doSynRateVariation");
    wts = generateFrequencyParameterization ("Paux", rate_classes);
 
    fprintf (stdout, "Node: ", ref);
    if (! do_srv) {
        fprintf (stdout,"\n\tLength parameter = ", ^(ref+".t"));
    }
             
    for (k = 1; k < rate_classes; k+=1) {
        Eval ("Paux" + k + "= ^(ref+\".Paux\" + k)");
    }

    for (k = 0; k < rate_classes; k+=1) {
        if (do_srv) {
            syn = ^(ref+".syn" + (k+1));
            nonsyn = ^(ref+".nonsyn" + (k+1));
            nananana = 1.0/3.0;
            omega = -1.0;
             if (syn != 0)
            {
                    _ratio = nonsyn/syn;
                    if (_ratio < 10000)
                    {
                            omega = Format (_ratio,10,4);
                    }
            }
            else
            {
                    if (nonsyn == 0)
                    {
                            omega = " Undefined";
                    }
                    else omega = "  Infinite";
            }

            fprintf (stdout, "\n\tClass ", k+1,
                     "\n\t\tsyn     : ", Format(^(ref+".syn" + (k+1)),10,4), 
                     "\n\t\tnon-syn : ", Format(^(ref+".nonsyn" + (k+1)),10,4), 
                     "\n\t\tomega   : ", omega,
    /*
                     "\n\t\tomega   : ", _standardizeRatio (^(ref+".nonsyn" + (k+1)),^(ref+".syn" + (k+1))),
        */
                     "\n\t\tweight  : ", Format(Eval(wts[k]),10,4));            
        } else {
            fprintf (stdout, "\n\tClass ", k+1,
                    "\n\t\tomega = ", Format(^(ref+".omega" + (k+1)),5,3), 
                    "\n\t\tweight = ", Format(Eval(wts[k]),5,3));
        }
    }
              
    fprintf (stdout, "\n"); 
    return 0;
}
