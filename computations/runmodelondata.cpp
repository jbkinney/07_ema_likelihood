/* runmodelondata.cpp - computes scores of all regions and 
 * locates all sites, as well as their respective energies
 *
 *  MATLAB usage: 
 *  [mi, regionStats, siteStats] = 
 *         runmodelondata(model, dataset, data, <maxNumSites>)
 *
 *
 */

#include "mex.h"
#include "mi.h"

#define MODEL           prhs[0]
#define DATASET         prhs[1]
#define DATA            prhs[2]
#define MAXNUMSITES     prhs[3]
#define POS1            prhs[4]
#define POS2            prhs[5]

#define MI              plhs[0]
#define REGIONSTATS     plhs[1]
#define SITESTATS       plhs[2]

// Main function
void mexFunction(   int nlhs,       mxArray *plhs[],
                    int nrhs, const mxArray *prhs[])
{
    int maxNumSites;
    bool useAllRegions;
    bool recordRegionScoresByBase = false;
    int pos1 = 0;
    int pos2 = 1;
    
    if ((nrhs<3) || (nrhs>6)) {
        mexErrMsgTxt("Funciton takes 3-6 arguments.");
    }
    if (nrhs==6) {
        pos1 = ((int) mxGetScalar(POS1)) - 1;
        pos2 = ((int) mxGetScalar(POS2)) - 1;
    }
    if (nrhs<4)
        maxNumSites = 0;
    else
        maxNumSites = (int) mxGetScalar(MAXNUMSITES);
    
    //////////////////////////////////////////////////////
    // Delclare variables
    int i, j, b, c, n, x, r, b1, b2;
    
    // model.emat related variables
    mxArray *EMAT;
    int L;
    double *emat;
    double *fwemat;
    double *rcemat;
    double e, maxe, mine, maxsitee, minsitee;
    double fwenergy, rcenergy, score;
    double scoreByBase[16];
    double maxScoreByBase[16];
    double *maxeByPos;
    int fwBaseIndex, rcBaseIndex;
    
    // model.cutoff related variables
    mxArray *CUTOFF;
    double *cutoffArray;
    double cutoff;
    
    // dataSet realted variables
    int dataSetN;
    int dataSetM;
    int numDataSets;
    int dataSetNum;
    
    // data.regions related variables
    int numRegions, regionIndex;
    mxArray *REGIONS;
    mxArray *REGIONINDICES;
    double *regionIndicesDbl;
    mxArray *REGIONTYPES;
    double *regionTypesDbl;
    mxArray *NUMTYPES;
    int numTypes;
    
    // data.regions.nseq related variables
    int seqLength, numForwardSites;
    mxArray *SEQ;
    double *seq; //int *seq;
    
    // regionStats related varibles
    const char *regionStatsFields[] = {"hits", "scores", "scoresByBase", "czx"};
    mxArray *REGIONHITS;
    double *regionHits;
    mxArray *REGIONSCORES;
    double *regionScores;
    mxArray *REGIONSCORESBYBASE;
    double *regionScoresByBase;
    mxArray *CZX;
    double *czx;
    //mxArray *REGIONSCORESBYBASESITENUM;
    //double *regionScoresByBaseSiteNum;
    int fwb1, fwb2, rcb1, rcb2;
    int fwSiteNum;
    int rcSiteNum = 0;
    
    // siteStats related variables
    int siteNum, numSites;
    const char *siteStatsFields[] = 
    {"orientations", "regions", "positions", "energies"};
    mxArray *SITEREGIONS;
    double *siteRegions = (double*) mxMalloc(maxNumSites*sizeof(double)); //new double[maxNumSites];
    double *srptr;
    mxArray *SITEPOSITIONS;
    double *sitePositions = (double*) mxMalloc(maxNumSites*sizeof(double)); //new double[maxNumSites];
    double *spptr;
    mxArray *SITEENERGIES;
    double *siteEnergies = (double*) mxMalloc(maxNumSites*sizeof(double)); //new double[maxNumSites];
    double *septr;
    mxArray *SITEORIENTATIONS;
    double *siteOrientations  = (double*) mxMalloc(maxNumSites*sizeof(double)); //new double[maxNumSites];
    double *soptr;
    
    // mi related variables 
    double *miPtr;

    // Get number of data sets
    dataSetN = (int) mxGetN(DATASET);
    dataSetM = (int) mxGetM(DATASET);
    numDataSets = dataSetN*dataSetM;
    const int dataSetDims[2] = {dataSetM, dataSetN};
    
    // Get model.emat
    EMAT = mxGetField(MODEL, 0, "emat");
    if (!EMAT)
        mexErrMsgTxt("Input 'model' does not have field 'emat'");
    if (mxGetN(EMAT) != 4)
        mexErrMsgTxt("Input 'model.emat' does not have 4 columns.");
    L = mxGetM(EMAT);
    emat = (double *) mxGetPr(EMAT);
    
    mxArray* FWEMAT = mxCreateDoubleMatrix(L,5,mxREAL);
    fwemat = (double *) mxGetPr(FWEMAT);
    mxArray* RCEMAT = mxCreateDoubleMatrix(L,5,mxREAL);
    rcemat = (double *) mxGetPr(RCEMAT);
    mxArray* MAXEBYPOS = mxCreateDoubleMatrix(L,1,mxREAL);
    maxeByPos = (double *) mxGetPr(MAXEBYPOS);
    
    //fwemat = new double[5*L];
    //rcemat = new double[5*L];
    //maxeByPos = new double[L];
    
    // Make sure positions are between 0 and L-1
    if ((pos1 < 0) || (pos1 >= L) || (pos2 < 0) || (pos2 >= L))
        mexErrMsgTxt("ERROR: Invalid pos1 or pos2");
    
    // Determine maximum and minimum site energies
    minsitee = 0.0;
    maxsitee = 0.0;
    for (i=0; i<L; i++) {
        maxe = emat[L*0 + i];
        mine = emat[L*0 + i];
        for (b=1; b<4; b++) {
            e = emat[L*b + i];
            if (e > maxe)
                maxe = e;
            if (e < mine)
                mine = e;
        }
        maxeByPos[i] = maxe;
        maxsitee += maxe;
        minsitee += mine;
    }
    
    // Determine maximum energy for site by base
    for (b1=0;b1<4;b1++) {
        for (b2=0;b2<4;b2++) {
            n = 4*b1+b2;
            maxScoreByBase[n] = maxsitee 
                - maxeByPos[pos1] + emat[L*b1+pos1] 
                - maxeByPos[pos2] + emat[L*b1+pos2];
        }
    }
    
    // Make two 5xL emats for forward and backward searching
    for (i=0; i<L; i++) {
        j=L-i-1;
        fwemat[L*0+i] = maxsitee - minsitee;    // Ensure bad site  score isn't counted
        rcemat[L*0+j] = maxsitee - minsitee;    // Ensure bad site  score isn't counted
        for (b=1; b<=4; b++) {
            c=5-b;
            fwemat[L*b + i] = emat[L*(b-1) + i];
            rcemat[L*c + j] = emat[L*(b-1) + i];
        }
    }
    
    // Get model.cutoff. Make sure it is the right size
    CUTOFF = mxGetField(MODEL, 0, "cutoff");
    if (!CUTOFF)
        mexErrMsgTxt("Input 'model' does not have field 'cutoff'");
    if ((mxGetM(CUTOFF) != dataSetM) || (mxGetN(CUTOFF) != dataSetN))
        mexErrMsgTxt("model.cutoff does not have the same size as dataSet");
    cutoffArray = (double *) mxGetPr(CUTOFF); 
    
    // Make sure data is the right size
    if ((mxGetM(DATA) != dataSetM) || (mxGetN(DATA) != dataSetN))
        mexErrMsgTxt("data does not have the same size as dataSet");
    
    // Create REGIONSTATS and SITESTATS
    MI = mxCreateDoubleMatrix(dataSetM, dataSetN, mxREAL);
    miPtr = (double *) mxGetPr(MI);
    REGIONSTATS = mxCreateStructArray(2, dataSetDims, 4, regionStatsFields);
    SITESTATS = mxCreateStructArray(2, dataSetDims, 4, siteStatsFields);
    
    for (dataSetNum=0; dataSetNum<numDataSets; dataSetNum++) {
    
        // Get dataSet.regionIndices
        REGIONINDICES = mxGetField(DATASET, dataSetNum, "regionIndices");
        if (!REGIONINDICES)
            mexErrMsgTxt("Input 'dataSet' does not have field 'regionIndices'");
        regionIndicesDbl = (double *) mxGetPr(REGIONINDICES);
    
        // Get dataSet.regionTypes
        REGIONTYPES = mxGetField(DATASET, dataSetNum, "regionTypes");
        if (!REGIONTYPES)
            mexErrMsgTxt("Input 'dataSet' does not have field 'regionTypes'");
        regionTypesDbl = (double *) mxGetPr(REGIONTYPES);
        
        // Get dataSet.types
        NUMTYPES = mxGetField(DATASET, dataSetNum, "numTypes");
        if (!NUMTYPES)
            mexErrMsgTxt("Input 'dataSet' does not have field 'numTypes'");
        numTypes = (int) mxGetScalar(NUMTYPES);
        
        // Make sure REGIONINDICES has the same number of elements as
        // REGIONTYPES
        if (mxGetNumberOfElements(REGIONINDICES) != 
                mxGetNumberOfElements(REGIONTYPES))
            mexErrMsgTxt("dataSet.regionIndices has a different size than dataSet.regionTypes");
            
    
        // Get data.regions
        REGIONS = mxGetField(DATA, dataSetNum, "regions");
        if (!REGIONS)
            mexErrMsgTxt("Input 'data' does not have field 'regions'");
    
        // Get numRegions
        numRegions = mxGetNumberOfElements(REGIONINDICES);
        
        // Declare hits
        REGIONHITS = mxCreateDoubleMatrix(numRegions,1,mxREAL);
        regionHits = (double *) mxGetPr(REGIONHITS);
    
        // Declare scores
        REGIONSCORES = mxCreateDoubleMatrix(numRegions,1,mxREAL);
        regionScores = (double *) mxGetPr(REGIONSCORES);
        
        // Declare scoresByBase
        REGIONSCORESBYBASE = mxCreateDoubleMatrix(numRegions,16,mxREAL);
        regionScoresByBase = (double *) mxGetPr(REGIONSCORESBYBASE);
        
        // Declare scoresByBaseSiteNum
        //REGIONSCORESBYBASESITENUM = mxCreateDoubleMatrix(numRegions,16,mxREAL);
        //regionScoresByBaseSiteNum = (double *) mxGetPr(REGIONSCORESBYBASESITENUM);
        
        // Get cutoff for this data set
        cutoff = cutoffArray[dataSetNum];
        
        // Compute scores of all regions
        siteNum = 0;
        for (r=0;r<numRegions;r++) {
            regionIndex = ((int) regionIndicesDbl[r])-1;
            
            SEQ = mxGetField(REGIONS,regionIndex,"nseq");
            if (!SEQ) {
                mexPrintf("Error in regions(model.regionIndex(%d)):\n", r+1);
                mexErrMsgTxt("Quiting at location 0.");
            }
            seq = (double *) mxGetPr(SEQ);
            seqLength = mxGetNumberOfElements(SEQ);
            
            // Make sure seq entries are kosher. 
            for (n=0;n<seqLength;n++)   {
                x = (int) seq[n];
                if ((x < 0) || (x > 4)) {
                    mexPrintf("ERROR: data(%d).region(%d).nseq(%d) is invalid.", dataSetNum, regionIndex, n);
                    mexErrMsgTxt("Quitting runmodelondata");
                }
            }    
        
            // Deterimine number of possible sites in the sequence
            numForwardSites = seqLength-L+1;
            if (numForwardSites < 0)
                numForwardSites = 0;
    
            // Set scores to maximum value possible
            score = maxsitee;
            for (n=0;n<16;n++) {
                scoreByBase[n] = maxScoreByBase[n];
            }
            
            // Compute energies of all valid sites
            for (n=0;n<numForwardSites;n++) {
                fwenergy = 0.0;
                rcenergy = 0.0;
            
                // Compute energy matrix contributions
                for (i=0;i<L;i++)   {
                    x = L*((int)seq[n+i])+i;
                    fwenergy += fwemat[x];
                    rcenergy += rcemat[x];
                } 
                
                fwSiteNum = 0;
                rcSiteNum = 0;
                
                // Record site if energy < cutoff
                if ((fwenergy < cutoff) && (siteNum < maxNumSites)) {
                    siteRegions[siteNum] = (double) (regionIndex+1);
                    sitePositions[siteNum] =(double) (n+1);
                    siteEnergies[siteNum] = fwenergy;
                    siteOrientations[siteNum] = 1.0;
                    fwSiteNum = siteNum+1;
                    siteNum++;
                }
            
                // Record site if energy < cutoff
                if ((rcenergy < cutoff) && (siteNum < maxNumSites)) {
                    siteRegions[siteNum] = (double) (regionIndex+1);
                    sitePositions[siteNum] =(double) (n+1);
                    siteEnergies[siteNum] = rcenergy;
                    siteOrientations[siteNum] = -1.0;
                    rcSiteNum = siteNum+1;
                    siteNum++;
                }
                

                // Replace scoreByBase (or score) if new energy is lower 
                // than current score. Check fwBaseIndex and rcBaseIndex too.
                fwBaseIndex = 4*(((int)seq[n+pos1])-1)+(((int)seq[n+pos2])-1);
                rcBaseIndex = 4*(4-((int)seq[n+L-(pos1+1)])) + (4-((int)seq[n+L-(pos2+1)]));
                if ((0 <= fwBaseIndex) && (fwBaseIndex < 16)) {
                    if (fwenergy < scoreByBase[fwBaseIndex]) {
                        scoreByBase[fwBaseIndex] = fwenergy;
                        //scoreByBaseSiteNum[fwBaseIndex] = (double) fwSiteNum;
                        if (fwenergy < score)
                            score = fwenergy;
                    }
                }
                else {
                    //mexPrintf("fwBaseIndex == %d\n", fwBaseIndex);
                    //mexPrintf("seq[n+pos1] == %d\n", seq[n+pos1]);
                    //mexPrintf("seq[n+pos2] == %d\n", seq[n+pos2]);
                    //mexErrMsgTxt("ERROR: INVALID fwBaseIndex");
                }
      
                if ((0 <= rcBaseIndex) && (rcBaseIndex < 16)) {
                    if (rcenergy < scoreByBase[rcBaseIndex]) {
                        scoreByBase[rcBaseIndex] = rcenergy;
                        //scoreByBaseSiteNum[rcBaseIndex] = (double) rcSiteNum;
                        if (rcenergy < score)
                            score = rcenergy;
                    }
                }
                else {
                    //mexPrintf("rcBaseIndex == %d\n", rcBaseIndex);
                    //mexPrintf("seq[n+L-(pos1+1)] == %d\n", seq[n+L-(pos1+1)]);
                    //mexPrintf("seq[n+L-(pos2+1)] == %d\n", seq[n+L-(pos2+1)]);
                    //mexErrMsgTxt("ERROR: INVALID rcBaseIndex");
                }
            }  
            
            // Record scores
            regionHits[r] = (score < cutoff) ? 1.0 : 0.0;
            regionScores[r] = score;
            for (n=0;n<16;n++) {
                regionScoresByBase[n*numRegions+r] = scoreByBase[n];
                //regionScoresByBaseSiteNum[n*numRegions+r] = scoreByBaseSiteNum[n];
            }
        }
    
        // Set site arrays to return to matlab. Make them the right length
        numSites = siteNum;
        SITEREGIONS = mxCreateDoubleMatrix(numSites,1,mxREAL);
        srptr = (double *) mxGetPr(SITEREGIONS);
        SITEPOSITIONS = mxCreateDoubleMatrix(numSites,1,mxREAL);
        spptr = (double *) mxGetPr(SITEPOSITIONS);
        SITEENERGIES = mxCreateDoubleMatrix(numSites,1,mxREAL);
        septr = (double *) mxGetPr(SITEENERGIES);
        SITEORIENTATIONS = mxCreateDoubleMatrix(numSites,1,mxREAL);
        soptr = (double *) mxGetPr(SITEORIENTATIONS);
        for (i=0; i<numSites; i++) {
            // Valid sites will all have regions > 0
            srptr[i] = siteRegions[i];
            spptr[i] = sitePositions[i];
            septr[i] = siteEnergies[i];
            soptr[i] = siteOrientations[i];
        }

        // Compute MI
        int x,z,k;
        CZX = mxCreateDoubleMatrix(numTypes,2,mxREAL);
        czx = (double *) mxGetPr(CZX);
        for (k=0;k<2*numTypes;k++) {
            czx[k] = 0.;
        }
        for (r=0;r<numRegions;r++) {
            x = int(regionHits[r]);
            z = int(regionTypesDbl[r]);
            czx[z+numTypes*x] += 1.0;
        }  
        miPtr[dataSetNum] = normalLogLikelihood(czx, numTypes, 2);
        
        // Copy back to REGIONSTATS and SITESTATS mxArray objects
        mxSetField(REGIONSTATS, dataSetNum, "scores", REGIONSCORES);
        mxSetField(REGIONSTATS, dataSetNum, "hits", REGIONHITS);
        mxSetField(REGIONSTATS, dataSetNum, "scoresByBase", REGIONSCORESBYBASE);
        mxSetField(REGIONSTATS, dataSetNum, "czx", CZX);
        //mxSetField(REGIONSTATS, dataSetNum, "scoresByBaseSiteNum", REGIONSCORESBYBASESITENUM);
        mxSetField(SITESTATS, dataSetNum, "regions", SITEREGIONS);
        mxSetField(SITESTATS, dataSetNum, "positions", SITEPOSITIONS);
        mxSetField(SITESTATS, dataSetNum, "energies", SITEENERGIES);
        mxSetField(SITESTATS, dataSetNum, "orientations", SITEORIENTATIONS); 
    }
    
    // Delete allocated arrays
    //delete[] fwemat;
    //delete[] rcemat;
    mxFree(siteRegions); //delete[] siteRegions;
    mxFree(sitePositions); //delete[] sitePositions;
    mxFree(siteEnergies); //delete[] siteEnergies;
    mxFree(siteOrientations); //delete[] siteOrientations;
    //delete[] maxeByPos;
}
