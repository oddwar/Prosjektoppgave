/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fishmodel;

import fishmodel.enkf.AssimSettings;
import fishmodel.enkf.EnsembleKF;
import fishmodel.enkf.Util;
import fishmodel.hydraulics.Balanced3DHydraulics;
import fishmodel.hydraulics.SimpleTankHydraulics;
import fishmodel.pellets.*;
import fishmodel.enkf.MpiHandler;

import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.text.SimpleDateFormat;
import java.util.*;

import fishmodel.sim.CurrentMagicFields;
import fishmodel.sim.InputDataNetcdf;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import save.SaveNetCDF;
import ucar.nc2.NetcdfFileWriteable;


/**
 *
 * @author malv
 */

/**
 * Setup for submergable seacage Scale AQ - Specialization project TTK 4551
 */
public class RunSimStein {

    public static final double HYPOXIA_THRESHOLD = 5.0976; //threshold value for low oxygen values - using Scale AQs values for reference 59% at 13 degrees

    public static LinearInterpolator interpol = new LinearInterpolator();


    public static void main(String[] args) {

        // Save files:
        String saveDir = "./";
        String simNamePrefix = "Projectassignment-test13-alternating saturation with 8 hours?-"; //"assim6_o2pert_lbeta_nopar_dropout";
        String simNamePostfix = "";

        boolean doMPI = false; // Will be set to true if we are running is EnKF mode using MPI
        AssimSettings as = new AssimSettings(); // Settings related to the EnKF are gathered in AssimSettings
        if (as.dryRun)
            simNamePrefix += "_dr";

        int rank = 0, N=1;
        MpiHandler mpi = null;
        EnsembleKF enKF = null;
        try {
            mpi = new MpiHandler(args);
            System.out.println("rank="+mpi.getRank()+", N="+mpi.getN());
            rank = mpi.getRank();
            N = mpi.getN();
            doMPI = true;
            // TEST TEST TEST:
            AdvectPellets.disableMultiprocessing(); // No internal parallelization when running MPI
            //
            if (as.usePerturbations && ((rank < N-1) || !as.useTwin || as.perturbTwin))
                as.perturbThisMember = true;
        } catch (Throwable ex) {
            ex.printStackTrace();
            System.out.println("Not doing MPI");
        }
        boolean isRoot = (rank==0); // For convenience, isRoot tells us if this is the rank 0 process.

/**
 * Changing this during testing.
 *maskO2Whensaving = false
 * All values outside grid = NaN
 */

        // Simulation settings:
        boolean maskO2WhenSaving = false;
        boolean varyAmbient = false; // Reduction in ambient values towards the rest of the farm
        //double addRedMult = 0.65*0.015; // Scale factor for reduction in ambient values

        boolean decreasingCurrentFactor = true; // Model gradual decrease in current factor due to
                                                // increasing cage net biofouling

        boolean useCurrentMagic = false; // Use spatially variable current flow field
        CurrentMagicFields cmf = null;
        if (useCurrentMagic) {
            cmf = new CurrentMagicFields("C:/Users/stein/Documents/01 Industriel kybernetikk/5. semester/Prosjektoppgave/Filer 03 aug/bjoroya_data.nc");
        }
        //C:\Users\stein\Documents\01 Industriel kybernetikk\5. semester\Prosjektoppgave\Filer 03 aug

        boolean use3dBalancedCurrent = false; // Use Balanced3DHydraulics

        boolean useVerticalDist = true; // Use non-uniform vertical distribution (defined further down)
                                        // for non-feeding fish

        boolean useInstantaneousAmbientVals = true; // true to use ambient value time series, false to use daily averages

        //boolean includeHypoxiaAvoidance = true;     int checkAvoidanceInterval = 30, checkAvoidanceCount = 0;

        // Simulation start time: Bjørøya
        int initYear = 2022, initMonth = Calendar.JUNE, initDate = 22, initHour = 0, initMin = 0, initSec = 0;
        //int initYear = , initMonth = Calendar.JUNE, initDate = 22, initHour = 0, initMin = 0, initSec = 0;
/**
 * Changing this during testing.
 * 17 hours first
 * lower resolution during testing
 *
 */
        double simHours = 20.;
        double t_end = simHours*3600;//1*24*3600; // Duration of simulation //
        int nSim = 1; // Number of days to simulate (separate sims)
        int startAt = 0; // Set to >0 to skip one of more simulations, but count them in the sim numbering

        // Common settings:
        double rad = 25;
        double depth = 50;
        double topOfCage = 12.5;
        double cylDepth = 32.5;
        double totDepth = 47.5; // Cage size (m)
        double dxy = 2, dz = dxy; // Model resolution (m) originally set to 2
        double dt = 0.5 * dxy; // Time step (s) originally set to 0.5 seconds
        int storeIntervalFeed = 60, storeIntervalInfo = 60;
        double fishMaxDepth = 47.5; // The maximum depth of the fish under non-feeding conditions

        //constant feeding added
        double nominalFeedingRate = 12771.2*1000/(simHours*3600); // updated to 9 hours // Approximate feeding over 10 hours based on calculated data 1.3 % of fishweight and max biomass 982 400 kg, oppdatert fra 0.67% etter verdier fra Scale AQ
                                                            // angis i antall gram


        double currentReductionFactor = (useCurrentMagic ? 1.0 : 0.8); // Multiplier for inside current
                                                                    // as function of outside


        // Environmental conditions:
        double currentSpeedInit = 0.20; // External current speed (m/s) // = 0,02 knop, no current. The current has upper bounds in this model
        double T_w = 13;    // water temp (celcius),
        //Todo make a formula depending on temperature and percentage
        double avO2 = 7.776; //  oxygen  (mg / l)  should bee 90% saturation determined by temperature
        double currentDirection = 25.; // direction of the current in degrees (0 - 360)
        boolean daysToSimulate = false; //is the simulation running over several days

        // Oxygen diffusion constant (values updated further down)
        double diffKappaO2 = 0.1, diffKappaO2Z = 0.1;

        double[] currentOffset = new double[] {0,0,0}; // Global current vector
        double[] currentOffset_r = new double[] {0,0,0}; // Perturbed global current vector

        // Fish setup (N, mean weight and std.dev weight):
        double nFishInit = 196480; // Biomass divided by 5 kg each salmon
        double meanWeight = 5000.0; // Estimated maximum weight at the harvest
        double[] wFish = new double[] {meanWeight, 0.2*meanWeight};
        
        // Wind speed (x, y components in m/s) affecting feed spreader:
        double[] windSpeed = new double[] {0, 0};
        // Pellet setup:
        double[] sizes = new double[] {3, 6, 9, 12};
        double[] speeds = new double[] {0.0773, 0.0815, 0.1284, 0.1421};
        int di = 2; // Index of chosen pellet size
        double pelletWeight=0.2; // Pellet weight (g)
        double sinkingSpeed = speeds[di];
        double kappa_ref = 0.00012;
        double kappa_add = 0.2;
        double refSize = 9;
        double diffKappa = kappa_ref*(kappa_add + Math.pow(sizes[di]/refSize,2));
        double kappa_z_mult = 25;
        double diffKappaZ = diffKappa*kappa_z_mult;

        // Set up cage dimensions and cage grid:
        double modelDim = 2*(rad+8*dxy);
        //double modelDim = 2*(rad+2.5*rad); // TEST TEST TEST extra padding
        int[] cageDims = new int[3];
        cageDims[0] = (int)Math.ceil(modelDim/dxy);
        cageDims[1] = cageDims[0];
        cageDims[2] = (int)Math.ceil(depth/dz)+1;
        boolean[][][] mask = null;
        //mask = CageMasking.circularMasking(cageDims, dxy, rad, false); // null

        mask = CageMasking.cylindroConicalMaskingSubmerged(cageDims, dxy, dz, rad, cylDepth, totDepth, topOfCage);
        //                  cylindroConicalMaskingSubmerged(int[] dims, double dxy, double dz, double radius, double cylDepth, double totDepth, double topDepth) {
        //
        boolean useWalls = false;
        System.out.println("Domain dimensions: ("+cageDims[0]+", "+cageDims[1]+", "+cageDims[2]+")");

        //todo adjust feeding position to submerged

        // Feeding setup:
        int[][] feedingPos = new int[][] {{cageDims[0]/2, cageDims[1]/2}};
        // Feeding periods (start/end in s):
        // Fra Eskil (Bjørøya): måltidene varte fra ca. kl. 07:30-17:30, i gjennomsnitt. Dette tilsvarer 7 dager med mating.

        /**
        int[][] feedingPeriods = new int[][] {{27000, 63000}, {86400+27000, 86400+63000}, {2*86400+27000, 2*86400+63000},
                {3*86400+27000, 3*86400+63000}, {4*86400+27000, 4*86400+63000}, {5*86400+27000, 5*86400+63000},
                {6*86400+27000, 6*86400+63000}, {7*86400+27000, 7*86400+63000}};
        */
        int[][] feedingPeriods = new int[][] {{0, 72000},
                {3*86400+27000, 3*86400+63000}, {4*86400+27000, 4*86400+63000}, {5*86400+27000, 5*86400+63000},
                {6*86400+27000, 6*86400+63000}, {7*86400+27000, 7*86400+63000}};

        int nPeriods = feedingPeriods.length;
        for (int i=0; i<nPeriods; i++) {
            System.out.println("Feeding period "+(i+1)+": "+feedingPeriods[i][0]+" to "+feedingPeriods[i][1]);
        }
        Object sourceTerm = null;

        double[] ambientValueFeed = new double[cageDims[2]]; // Outside feed concentrations are set to 0
        for (int i = 0; i < ambientValueFeed.length; i++) {
            ambientValueFeed[i] = 0;
        }

        // Set up array for ambient temperature:
        double[] ambientTemp = new double[cageDims[2]];
        for (int i = 0; i < ambientTemp.length; i++) {
            ambientTemp[i] = T_w;
        }

        // Determine number of fish, and feeding rate:
        double nFish = nFishInit;
        System.out.println("N fish = "+nFish);

        double feedingRateMult = 0; // Set each timestep
        System.out.println("Feeding rate = "+nominalFeedingRate);

        // Oxygen sensor positions:
        Measurements.MeasurementSet ms = Measurements.setupSensorPositionsBjoroya(cageDims, dxy, dz, rad);
        String[] o2Names = ms.names;
        int[][] o2Pos = ms.pos;
        // List O2 sensor positions and which are masked / not masked:
        for (int i=0; i<o2Names.length; i++) {
            System.out.println(o2Names[i]+": "+o2Pos[i][0]+" , "+o2Pos[i][1]+" , "+o2Pos[i][2]+", mask="+mask[o2Pos[i][0]][o2Pos[i][1]][o2Pos[i][2]]);
        }

        // Feed affinity (determines how "interested" the fish is in feeding in each model grid cell):
        double[][][] affinity = new double[cageDims[0]][cageDims[1]][cageDims[2]];
        // O2 affinity (determines the typical vertical distribution of the fish):
        double[][][] o2Affinity = new double[cageDims[0]][cageDims[1]][cageDims[2]];

        //todo Add depth adjustment to distribution of fish
        double depthAdjustment = 1.;
        // Define vertical distribution based on telemetry data (from 8 individuals):
        double[] affProfile_orig = new double[] {0.0110 , 0.0913 , 0.8601, 2.1406, 2.7774, 2.6903, 2.5195, 2.2987, 2.0137,
                1.7448, 1.5883, 1.3667, 1.2348, 1.0724, 0.9379, 0.7764, 0.7104, 0.5895, 0.5607, 0.4668, 0.3933,
                0.4009, 0.2935, 0.1801, 0.1260, 0.0787, 0.0457, 0.0304};

        double[] affProfile_flat = new double[] {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        // Choose the flat or variable profile according to the "useVerticalDist" setting:
        double[] affProfile = useVerticalDist ? affProfile_orig : affProfile_flat;
        // Depths values for the values defined above:
        // distributing the fish

        double depthDiff = topOfCage;
        double[] affDepths = new double[] {0.5000+depthDiff, 1.5000+depthDiff, 2.5000+depthDiff, 3.5000+depthDiff, 4.5000+depthDiff, 5.5000+depthDiff,
                6.5000+depthDiff, 7.5000+depthDiff, 8.5000+depthDiff, 9.5000+depthDiff, 10.5000+depthDiff, 11.5000+depthDiff, 12.5000+depthDiff,
                13.5000+depthDiff, 14.5000+depthDiff, 15.5000+depthDiff, 16.5000+depthDiff, 17.5000+depthDiff, 18.5000+depthDiff, 19.5000+depthDiff,
                20.5000+depthDiff, 21.5000+depthDiff, 22.5000+depthDiff, 23.5000+depthDiff, 24.5000+depthDiff, 25.5000+depthDiff, 26.5000+depthDiff,
                27.5000+depthDiff};
        double[] affinityProfile = new double[cageDims[2]];
        interpolateVertical(affinityProfile, affDepths, affProfile, cageDims[2], dz);

        // Set feeding and O2 affinity based on the chosen vertical profile and cage mask:
        double o2AffSum = setO2AffinityWithVerticalProfile(cageDims, dz, fishMaxDepth, mask, affinityProfile, o2Affinity, affinity);
        //double o2AffSum = setO2AffinityWithVerticalProfileAndEdgeDecrease(cageDims, rad, dxy, dz, fishMaxDepth, mask, affinityProfile, o2Affinity, affinity);
        int availableCellsForO2Uptake = countAvailableCellsForOxygenUptake(cageDims, dz, fishMaxDepth, mask);

        // Oxygen
        double[] ambientValueO2 = new double[cageDims[2]]; // Ambient value of O2
        double[] ambientValueO2_r = new double[cageDims[2]]; // Possibly perturbed ambient value of O2
        for (int i = 0; i < ambientValueO2.length; i++) {
            ambientValueO2[i] = avO2;
        }

        // If we are running in EnKF mode, let rank 0 initialize the EnKF class:
        if (doMPI && (rank==0)) {
            enKF = new EnsembleKF(simNamePrefix, cageDims, as.nPar, dxy, ms);
        }

        // Set up initial perturbation and parameter values:
        double ambientO2_perturb = 0;
        double[] current_perturb = new double[2];
        double o2Cons_perturb = 0; // Relative perturbation to total oxygen consumption to be updated per time step
        double o2Cons_perturb_r = 0; // The perturbation to apply at this particular time step - may be set
            // equal to o2Cons_perturb, or to the sum of o2Cons_perturb and an estimated consumption parameter.
        double[] parVal = new double[] {0, 0, 0, 0}; //[nPar];

        // Initialize number formatter:
        NumberFormat nf1 = NumberFormat.getNumberInstance(Locale.ENGLISH);
        nf1.setMaximumFractionDigits(1); nf1.setMinimumFractionDigits(0);
        // Initialize random number generator:
        Random rnd = new Random();
        // Other initialization:
        double lastMeanFeedDepth = -1;


        // --------------------------------------------------------------------------
        // Simulations to run
        // --------------------------------------------------------------------------

        for (int sim=0; sim<nSim; sim++) {
            if (sim<startAt)
                continue;
            int feedingPeriodPiv = 0;
            boolean isFeeding = false;

            Calendar c = Calendar.getInstance();
            c.set(initYear, initMonth, initDate, initHour, initMin, initSec);
            c.add(Calendar.DAY_OF_MONTH, sim);
            Date startTime = c.getTime();


            if(sim==2 && sim !=1 && sim!=0 && daysToSimulate == true){
                T_w = 8.0;
                avO2 = 9.52; // 80% saturation at 8.0 celsius
                // Set up array for ambient temperature:
                for (int i = 0; i < ambientTemp.length; i++) {
                    ambientTemp[i] = T_w;
                }

                for (int i = 0; i < ambientValueO2.length; i++) {
                    ambientValueO2[i] = avO2;
                }
            }

            // Current field
            double[][][][] hydro;
            double[][] currentProfile = new double[cageDims[2] + 1][3];

            //**
            // Current profile from the Aquaexcel project
            double[] initialLowCurrent = {0.041, 0.03, 0.021, 0.015, 0.0145, 0.013, 0.012};
            double[] currentDepths  = {2.5, 5, 7.5, 12.5, 17.5, 27.5, 32.5};
            double[] lowCurrent = new double[initialLowCurrent.length];

            double currentMult = currentSpeedInit/0.041;
            double currentMult_2 = 0.8 * currentMult; // 16 cm/s
            double currentMult_3 = 0.6 * currentMult; // 12 cm/s
            double currentMult_4 = 0.4 * currentMult; // 8 cm/s
            double currentMult_5 = 0.2 * currentMult; // 4 cm/s



            for (int i = 0; i < initialLowCurrent.length; i++){
                lowCurrent[i] = initialLowCurrent[i]*currentMult;
            }



            // Here you set up the current profile (3D current vector per depth layer):
            // currentProfile[0] = x, currentProfile[1] = y currentProfile[2] = z

            double [] tempCurrentProfile = new double[cageDims[2]+1];
            interpolateVertical(tempCurrentProfile, currentDepths, lowCurrent, cageDims[2], dz);

            // Due to numerical errors in the model it is best not to have current from a straight direction
            double theta = currentDirection*Math.PI/180.;

            for (int i = 0; i< tempCurrentProfile.length;i++) {
                currentProfile[i][0] = tempCurrentProfile[i]*Math.cos(theta);
                currentProfile[i][1] = tempCurrentProfile[i]*Math.sin(theta);
            }

            //*/

            /** //original currentProfile
            for (int i=0; i<cageDims[2]+1; i++) {

                //updating currentprofile from 0 to CurrentSpeedInit
                currentProfile[i][0] = currentSpeedInit;
            }
            */


            //Hydro is a 4-d matrix with current for x,y,z direction.
            hydro = SimpleTankHydraulics.getProfileHydraulicField(cageDims, currentProfile);

            //Todo: Do i need this for two days of simulation?
            if (decreasingCurrentFactor)
                currentReductionFactor = 0.8 + 0.05 - ((double)sim)*(0.2/*0.25*//8.0);

            AdvectPellets ap = new AdvectPellets();
            AdvectPellets apOx = new AdvectPellets();
            if (varyAmbient) {
                apOx.setVaryAmbient(true, affinityProfile);

            }

            // Initialize environmental input data:
            //String inDataFile = "C:/Users/alver/OneDrive - NTNU/prosjekt/O2_Bjørøya/bjoroya_data.nc";
            String inDataFile = "C:/Users/stein/Documents/01 Industriel kybernetikk/5. semester/Prosjektoppgave/1. Simuleringer og plotting/bjoroya_data.nc";
            if (!(new File(inDataFile)).exists())
                inDataFile = "bjoroya_data.nc";
            InputDataNetcdf inData = new InputDataNetcdf(inDataFile, useInstantaneousAmbientVals);
            inData.setStartTime(startTime);


            // Format a unit string for the time variable to save to NetCDF giving the initial time:
            SimpleDateFormat formatter = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
            String unitString = "seconds since "+formatter.format(startTime);
            System.out.println("Unit string: "+unitString);

            NumberFormat nf = NumberFormat.getNumberInstance(Locale.US);
            nf.setMaximumFractionDigits(2);
            nf.setMinimumFractionDigits(2);

            // Initialize states:
            double[][][] fc = new double[cageDims[0]][cageDims[1]][cageDims[2]];
            double[][][] o2 = new double[cageDims[0]][cageDims[1]][cageDims[2]];
            double[][][] ingDist = new double[cageDims[0]][cageDims[1]][cageDims[2]];
            double[][][] o2consDist = new double[cageDims[0]][cageDims[1]][cageDims[2]];
            AdvectPellets.initField(o2, avO2);
            AdvectPellets.initField(ingDist, 0);
            double outFlow = 0., outFlow_net = 0.;

/**
 * Need to give initial values manually in order to test the limits
  */

            // Initialize O2 field based on first ambient values:
            double[] ambVal = new double[] { avO2,avO2 ,avO2 };
            //double[] ambVal = new double[] {inData.getO2Ambient5(), inData.getO2Ambient10(), inData.getO2Ambient15()};
            interpolateVertical(ambientValueO2, new double[] {5, 10, 15}, ambVal, cageDims[2], dz);
            for (int i=0; i<cageDims[0]; i++)
                for (int j=0; j<cageDims[1]; j++)
                    for (int k=0; k<cageDims[2]; k++) {
                        o2[i][j][k] = ambientValueO2[k];
                    }

            // Setup of surface feeding:
            double[][][] feedingRate = new double[cageDims[0]][cageDims[1]][cageDims[2]];
            double[][] fTemp = new double[cageDims[0]][cageDims[1]];

            double[][] surfFeed = new double[cageDims[0]][cageDims[1]];
            for (int i=0; i<feedingPos.length; i++) {
                PelletSpreaderModel.setPelletDist(fTemp, feedingPos[i][0], feedingPos[i][1], dxy, 0, 35, 0, true, windSpeed, null);
                // Add to total distribution:
                for (int ii=0; ii<surfFeed.length; ii++)
                    for (int jj = 0; jj < surfFeed[ii].length; jj++) {
                        surfFeed[ii][jj] = surfFeed[ii][jj] + fTemp[ii][jj];
                    }
            }

            for (int i=0; i<fTemp.length; i++)
                for (int j=0; j<fTemp[i].length; j++) {
                    feedingRate[i][j][0] = fTemp[i][j];
                }
            sourceTerm = feedingRate;

            // Initialize simple (grouped) fish model:
            SimpleFish fish = new SimpleFish(nFish, wFish[0], wFish[1]);
            double[][][] fishTmp = new double[fish.getNGroups()][1][1];

            SimpleDateFormat filenameForm = new SimpleDateFormat("dd_MM");
            String filePrefix = simNamePrefix+filenameForm.format(startTime);

            // Establish file names to write data to:
            NetcdfFileWriteable ncfile = null;
            NetcdfFileWriteable fishfile = null;
            // gave new name
            String ncfilePath = saveDir + filePrefix + simNamePostfix + (doMPI ? "_"+String.format("%02d", rank) : "")+".nc";
            String fishfilePath = saveDir + filePrefix + simNamePostfix + (doMPI ? "_"+String.format("%02d", rank) : "")+"_fish.nc";
            boolean firstStore3d = true, firstStoreScalars = true;

            double totFeedAdded = 0;

            //double currentDirection = 90;
            double currentSpeed = currentSpeedInit;
            double f_currentDir = Math.exp(-0.01*dt); // For simulating Gauss-Markov process
            double sigma_currentDir = 10;


            double t = 0;
            int n_steps = (int) (t_end / dt);
            long stime = System.currentTimeMillis();
            boolean currentReduction_1 = true; // true - enables first current reduction
            boolean currentReduction_2 = true;
            boolean currentReduction_3 = true;
            boolean currentReduction_4 = true;



//--------- start of loop -------------------------------------------------------------------------------------------------------------

            for (int i = 0; i < n_steps; i++) {

                //System.out.println("t = "+t);
                double tMin = t / 60;

  /**
   //* Using manually entered data
                //------------- vv Current magic vv --------------------------------------------------------------------------------------
                //                  Noe er moffens her

                if (inData.advance(t) || (i==0)) {
                    //double[] tempVal = {inData.getTemperature5(), inData.getTemperature10(), inData.getTemperature15()};
                    double[] tempVal = {T_w, T_w, T_w};
                    //tempVal[0] = 18; tempVal[1] = 10; tempVal[2] = 2;
                    interpolateVertical(ambientTemp, new double[] {5, 10, 15}, tempVal, cageDims[2], dz);
                    //ambVal = new double[] {inData.getO2Ambient5(), inData.getO2Ambient10(), inData.getO2Ambient15()};
                    ambVal = new double[] {avO2 , avO2, avO2};
                    interpolateVertical(ambientValueO2, new double[] {5, 10, 15}, ambVal, cageDims[2], dz);

                    //currentSpeed = inData.getExtCurrentSpeed();
                    //currentDirection = inData.getExtCurrentDir();

                    double profileValCur = avO2;
                    double profileValDir = currentDirection;
                    int profileMult = (int) totDepth ;

                    double[] obsCurrentProfile = inData.getExtCurrentSpeedProfile();
                    //double[] obsCurrentProfile = makeProfile(profileMult, profileValCur);
                    double[] obsCurrentDirProfile = inData.getExtCurrentDirProfile();
                    //double[] obsCurrentDirProfile = makeProfile(profileMult,profileValDir);

                    double[] obsCurrentComp1 = new  double[obsCurrentProfile.length],
                             obsCurrentComp2 = new  double[obsCurrentProfile.length];
                    // Current directions are given as the direction the current flows towards, with
                    // 0 degrees being north and 90 degrees being east. Verified by comparing histograms
                    // with the textual descriptions in the report by Aqua Kompetanse.
                    // x component: speed*sin(direction)
                    // y component: speed*cos(direction)


                    for (int j = 0; j < obsCurrentComp1.length; j++) {
                        double currentReductionFactorHere = currentReductionFactor;
                        // If we are using 3D balanced currents, use a factor of 1 here,
                        // so we keep the original ambient currents. The reduction factor is instead
                        // taken into account when the 3D field is calculated.
                        if (use3dBalancedCurrent)
                            currentReductionFactorHere = 1.0;

                        //currentReductionFactorHere *= Math.min(1.0, (1.+12.*(obsCurrentProfile[j]-0.06))); // TEST
                        obsCurrentComp1[j] = currentReductionFactorHere*
                                obsCurrentProfile[j]*Math.sin(obsCurrentDirProfile[j]*Math.PI/180.);
                        obsCurrentComp2[j] = currentReductionFactorHere*
                                obsCurrentProfile[j]*Math.cos(obsCurrentDirProfile[j]*Math.PI/180.);
                    }

                    double[] obsCurrentDepths = inData.getCurrentDepths();
                    // since the current is uniform the depth is irrelevant, it is relevant that it covers the full depth

                    //double[] obsCurrentDepths = makeDepthProfile(totDepth);
                    double[] interpProfile1 = new double[cageDims[2]],
                             interpProfile2 = new double[cageDims[2]];
                    interpolateVertical(interpProfile1, obsCurrentDepths, obsCurrentComp1, cageDims[2], dz);
                    interpolateVertical(interpProfile2, obsCurrentDepths, obsCurrentComp2, cageDims[2], dz);




                    if (!useCurrentMagic) {
                        for (int j = 0; j < interpProfile1.length; j++) {
                            //double speedHereNow = Math.sqrt(interpProfile1[j]*interpProfile1[j] + interpProfile2[j]*interpProfile2[j]);
                            //currentProfile[j][0] = 0;//0.1317*speedHereNow;
                            //currentProfile[j][1] = 0;//0.9913*speedHereNow;
                            currentProfile[j][0] = interpProfile1[j];
                            currentProfile[j][1] = interpProfile2[j];
                            currentProfile[j][2] = 0.;
                        }
                        // Update current field using the new profile:
                        if (!use3dBalancedCurrent)
                            // Create uniform current field for each layer:
                            SimpleTankHydraulics.getProfileHydraulicField(hydro, cageDims, currentProfile);
                        else
                            // Create a balanced 3D current field based on the current profile:
                            //Balanced3DHydraulics.getBjoroyaHydraulicField(cageDims, dxy, rad,
                            //        currentReductionFactor, currentProfile, hydro);
                            Balanced3DHydraulics.getTurbulentHydraulicField(cageDims, dxy, rad,
                                    currentReductionFactor, currentProfile, hydro);
                            //Balanced3DHydraulics.stats(hydro);
                    }
                    else {
                        double[] lDirections = new double[cageDims[2]],
                                lSpeeds = new double[cageDims[2]];
                        for (int j = 0; j < cageDims[2]; j++) {
                            lSpeeds[j] = Math.sqrt(interpProfile1[j]*interpProfile1[j] + interpProfile2[j]*interpProfile2[j]);
                            lDirections[j] = Math.atan2(interpProfile1[j], interpProfile2[j])*180./Math.PI;

                        }
                        cmf.setCurrentField(hydro, lSpeeds, lDirections);

                    }
                }
  //*/ //------------- ^^ Current magic ^^ --------------------------------------------------------------------------------------

                // Update variable current speed offset:
                currentOffset[0] = 0.;//currentReductionFactor*currentSpeed*Math.cos(currentDirection*Math.PI/180.);
                currentOffset[1] = 0.;//currentReductionFactor*currentSpeed*Math.sin(currentDirection*Math.PI/180.);

                //diffKappaO2 = 0.2*Math.min(0.5, 10*Math.pow(currentReductionFactor*0.06,2)); // Math.min(0.5, 10*Math.pow(currentReductionFactor*0.04,2));
                //diffKappaO2Z = 5.0*0.1*Math.min(0.5, 10*Math.pow(currentReductionFactor*0.06,2)); // Math.min(0.5, 10*Math.pow(currentReductionFactor*0.04,2));
                diffKappaO2 = Math.min(0.5, 10*Math.pow(currentReductionFactor*0.06,2)); // Math.min(0.5, 10*Math.pow(currentReductionFactor*0.04,2));
                diffKappaO2Z = 5.0*0.1*Math.min(0.5, 10*Math.pow(currentReductionFactor*0.06,2)); // Math.min(0.5, 10*Math.pow(currentReductionFactor*0.04,2));
                //System.out.println("DiffKappa O2: "+diffKappaO2);

                // Update feeding rate depending on preset feeding periods:
                if (!isFeeding) { // Not already feeding. Check if we should start:
                    if ((feedingPeriodPiv < feedingPeriods.length) && (t >= feedingPeriods[feedingPeriodPiv][0])) {
                        isFeeding = true;
                        feedingRateMult = nominalFeedingRate;
                        // Reset gut content at start of feeding period:
                        fish.resetAllV();
                    } else
                        feedingRateMult = 0;
                } else { // Already feeding. Check if we should stop:
                    if (t >= feedingPeriods[feedingPeriodPiv][1]) {
                        isFeeding = false;
                        feedingRateMult = 0;
                        feedingPeriodPiv++; // Update so we start looking for next period
                    } else
                        feedingRateMult = nominalFeedingRate;
                }

                /*if (includeHypoxiaAvoidance) {
                    checkAvoidanceCount++;
                    if (checkAvoidanceCount == checkAvoidanceInterval) {
                        checkAvoidanceCount = 0;
                        o2AffSum = setO2AffinityWithAvoidance(cageDims, dz, fishMaxDepth, mask, o2, affinity, o2Affinity);
                        //System.out.println("Updating O2 affinity. sum="+o2AffSum);

                    }
                }*/

                //** Alternating the current after 4 hours = 14 400 sec
                if (currentReduction_1 && t/14400 == 1 ){

                    //Run new current profile set up with reduction.
                    for (int o = 0; o < initialLowCurrent.length; o++){
                        lowCurrent[o] = initialLowCurrent[o]*currentMult_2; // 11cm/S
                    }
                    interpolateVertical(tempCurrentProfile, currentDepths, lowCurrent, cageDims[2], dz);
                    for (int k = 0; k< tempCurrentProfile.length;k++) {
                        currentProfile[k][0] = tempCurrentProfile[k]*Math.cos(theta);
                        currentProfile[k][1] = tempCurrentProfile[k]*Math.sin(theta);
                    }

                    hydro = SimpleTankHydraulics.getProfileHydraulicField(cageDims, currentProfile);
                    currentReduction_1 = false;
                }

                // Alternating the current after 8 hours = 28 800 sec
                if (currentReduction_2 && t/28800 == 1 ){

                    //Run new current profile set up with reduction.
                    for (int o = 0; o < initialLowCurrent.length; o++){
                        lowCurrent[o] = initialLowCurrent[o]*currentMult_3; // 9,5cm/S
                    }
                    interpolateVertical(tempCurrentProfile, currentDepths, lowCurrent, cageDims[2], dz);
                    for (int k = 0; k< tempCurrentProfile.length;k++) {
                        currentProfile[k][0] = tempCurrentProfile[k]*Math.cos(theta);
                        currentProfile[k][1] = tempCurrentProfile[k]*Math.sin(theta);
                    }

                    hydro = SimpleTankHydraulics.getProfileHydraulicField(cageDims, currentProfile);
                    currentReduction_2 = false;
                }

                // Alternating the current after 12 hours = 43 200 sec
                if (currentReduction_3 && t/43200 == 1 ) {

                    //Run new current profile set up with reduction.
                    for (int o = 0; o < initialLowCurrent.length; o++) {
                        lowCurrent[o] = initialLowCurrent[o] * currentMult_4; // 8 cm/s
                    }
                    interpolateVertical(tempCurrentProfile, currentDepths, lowCurrent, cageDims[2], dz);
                    for (int k = 0; k < tempCurrentProfile.length; k++) {
                        currentProfile[k][0] = tempCurrentProfile[k] * Math.cos(theta);
                        currentProfile[k][1] = tempCurrentProfile[k] * Math.sin(theta);
                    }


                    hydro = SimpleTankHydraulics.getProfileHydraulicField(cageDims, currentProfile);
                    currentReduction_3 = false;
                }

                // Alternating the current after 16 hours = 57 600 sec
                if (currentReduction_4 && t/57600 == 1 ) {

                    //Run new current profile set up with reduction.
                    for (int o = 0; o < initialLowCurrent.length; o++) {
                        lowCurrent[o] = initialLowCurrent[o] * currentMult_5; // 6.5 cm/s
                    }
                    interpolateVertical(tempCurrentProfile, currentDepths, lowCurrent, cageDims[2], dz);
                    for (int k = 0; k < tempCurrentProfile.length; k++) {
                        currentProfile[k][0] = tempCurrentProfile[k] * Math.cos(theta);
                        currentProfile[k][1] = tempCurrentProfile[k] * Math.sin(theta);
                    }

                    hydro = SimpleTankHydraulics.getProfileHydraulicField(cageDims, currentProfile);
                    currentReduction_4 = false;
                }

                    //*/



                totFeedAdded += dt * feedingRateMult;

                o2Cons_perturb_r = 0.; //initially zero

                // Perturb if we are using MPI, except if we are using a twin, and this is the twin, and the
                // twin is not to be perturbed.
                if (doMPI && as.perturbThisMember) {
                    ambientO2_perturb = Util.updateGaussMarkov(ambientO2_perturb, as.ambientO2Beta, as.ambientO2Std, dt, rnd);
                    for (int j = 0; j < ambientValueO2.length; j++) {
                        ambientValueO2_r[j] = ambientValueO2[j] + ambientO2_perturb;
                    }
                    for (int j = 0; j < 2; j++) {
                        current_perturb[j] = Util.updateGaussMarkov(current_perturb[j], as.currentBeta, as.currentStd, dt, rnd);
                        currentOffset_r[j] = currentOffset[j] + current_perturb[j];
                    }
                    o2Cons_perturb = Util.updateGaussMarkov(o2Cons_perturb, as.o2ConsBeta, as.o2ConsStd, dt, rnd);
                    o2Cons_perturb_r = o2Cons_perturb;
                } else if (as.useTwin && (rank == N-1)) {
                    // This is the twin, introduce possible model error here.
                    for (int j = 0; j < ambientValueO2.length; j++) {
                        ambientValueO2_r[j] = ambientValueO2[j] + 0.25;

                    }
                    System.arraycopy(currentOffset, 0, currentOffset_r, 0, currentOffset.length);
                } else {
                    // Copy ambientValueO2 and currentOffset without perturbations:
                    System.arraycopy(ambientValueO2, 0, ambientValueO2_r, 0, ambientValueO2.length);
                    System.arraycopy(currentOffset, 0, currentOffset_r, 0, currentOffset.length);
                }

                // Perturb parameters according to their std. setting if we have any:
                if (doMPI && as.perturbThisMember && as.nPar > 0) {

                    for (int j=0; j<as.nPar; j++) {
                        parVal[j] += as.parStd[j]*dt*rnd.nextGaussian();
                    }

                    // The first three parameters are perturbations to ambient O2 at 5, 10 and 15 m. We need to
                    // calculate a linear interpolation of these to all model depths before applying it:
                    double[] ambO2Par = new double[] {parVal[0], parVal[1], parVal[2]}; // Make array of amb O2 related parameters
                    double[] interpAmbO2Par = new double[cageDims[2]];
                    interpolateVertical(interpAmbO2Par, as.parDepths, ambO2Par, cageDims[2], dz);

                    // Apply parameter values to model:
                    // Param 0: offset to ambient O2 values:
                    for (int j=0; j<cageDims[2]; j++)
                        ambientValueO2_r[j] += interpAmbO2Par[j];
                    /*for (int j = 0; j < ambientValueO2.length; j++) {
                        ambientValueO2_r[j] += parVal[0];
                    }*/

                    // Parameter number 4 is additional perturbation to total o2 consumption:
                    o2Cons_perturb_r += parVal[3];
                }

                double[] r = ap.step(dt, fc, dxy, dz, useWalls, mask, sinkingSpeed, diffKappa, diffKappaZ,
                        hydro, currentOffset_r, sourceTerm, feedingRateMult, ambientValueFeed);
                outFlow = r[0]; // Feed lost from grid (not used)
                outFlow_net = r[1]; // Feed lost from the unmasked part of the grid (feed lost through side)

                double[] o2OutFlow = apOx.step(dt, o2, dxy, dz, useWalls, mask, 0, diffKappaO2, diffKappaO2Z,
                        hydro, currentOffset_r, feedingRate, 0, ambientValueO2_r);



                double[] res = IngestionAndO2Tempprofile.calculateIngestion(dt, fc, o2, affinity, o2Affinity, o2AffSum,
                        availableCellsForO2Uptake, ingDist, o2consDist, dxy, dz, mask, pelletWeight, ambientTemp, fish,
                        o2Cons_perturb_r);
                double totalIntake = res[0], rho = res[1], o2ConsumptionRate = res[2];

                t = t + dt;

                // Check if we are running with EnKF. If so, check if it is time to do model correction:
                if (doMPI && i>0 && ((t/((double)as.enKFInterval) - Math.floor(t/(double)as.enKFInterval)) < 1e-5)) {
                    double[][] X = mpi.gatherStateToRank0(o2, parVal);
                    if (isRoot) {
                        System.out.println("Calling EnKF");
                        long tic = System.currentTimeMillis();
                        double[][] X_a = enKF.doAnalysis(t, X, as.useTwin, as.locDist/dxy, as.locZMultiplier,
                                as.ensembleInflation, as.ensembleInflationFactor, inData);
                        long duration = System.currentTimeMillis() - tic;
                        if (duration > 1000L)
                            System.out.println("Analysis took "+(duration/1000L)+" seconds.");
                        else
                            System.out.println("Analysis took "+duration+" ms.");
                        if (!as.dryRun && !as.isDropOutActive(t)) { // We only apply corrections if we are not doing a dry run and not in a dropout interval
                            mpi.distributeAnalysisFromRank0(X_a, o2, parVal, cageDims, as.nPar);
                        }
                    } else if (!as.dryRun && !as.isDropOutActive(t)) { // We only apply corrections if we are not doing a dry run and not in a dropout interval
                        if (!as.useTwin || (rank < N-1)) {
                            mpi.receiveAnalysisFromRank0(o2, parVal, cageDims, as.nPar);
                        }
                    }
                }

                // Check if it is time to store 3D fields of feed and O2:
                if (i>0 && ((t/((double)storeIntervalFeed) - Math.floor(t/(double)storeIntervalFeed)) < 1e-5)) {
                    if (isRoot) {
                        double elapsed = (double) ((System.currentTimeMillis() - stime)) / 60000.;
                        double fractionCompleted = ((double) i) / ((double) n_steps);
                        double remaining = (elapsed / fractionCompleted) - elapsed;
                        System.out.println("t = " + nf1.format(t) + " - Estimated time to complete: " + nf1.format(remaining) + " minutes");
                    }

                    if (firstStore3d) {
                        firstStore3d = false;
                        ncfile = SaveNetCDF.initializeFile(ncfilePath, cageDims, 1, 1, unitString, ms);
                        SaveNetCDF.createCageVariables(ncfile, "feed", "ingDist", "o2", "o2consDist");
                    }
                    else {
                        try {
                            ncfile = NetcdfFileWriteable.openExisting(ncfilePath);
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                    }

                    SaveNetCDF.saveCageVariable(ncfile, t, "feed", fc, mask, true);
                    SaveNetCDF.saveCageVariable(ncfile, t, "ingDist", ingDist, mask, false);
                    SaveNetCDF.saveCageVariable(ncfile, t, "o2", o2, (maskO2WhenSaving ? mask : null), false);
                    SaveNetCDF.saveCageVariable(ncfile, t, "o2consDist", o2consDist, mask, false);

                    try {
                        ncfile.close();
                    } catch (IOException e) {
                        e.printStackTrace();
                    }

                }

                // Check if it is time to store scalar output values:
                if (i>0 && ((t/((double)storeIntervalInfo) - Math.floor(t/(double)storeIntervalInfo)) < 1e-5)) {

                    if (firstStoreScalars) {
                        firstStoreScalars = false;
                        fishfile = SaveNetCDF.initializeFile(fishfilePath, new int[]{fish.getNGroups(), 1, cageDims[2]}, 1, 1, unitString, ms);
                        SaveNetCDF.createProfileVariable(fishfile, "appetite", 0);
                        SaveNetCDF.createProfileVariable(fishfile, "ingested", 0);
                        SaveNetCDF.createProfileVariable(fishfile, "V", 0);
                        SaveNetCDF.createScalarVariables(fishfile, "rho", "feedingRate", "o2ConsumptionRate",
                                "min_O2", "mean_O2", "frac_hypoxia",
                                "meanFeedDepth", "d_meanFeedDepth", "totIngRate", "totIngested", "totFeed",
                                "waste", "waste_net");
                        SaveNetCDF.createProfileVariable(fishfile, "ext_O2", 2); // dim=2 means along z dim
                        SaveNetCDF.createProfileVariable(fishfile, "temperature", 2);
                        SaveNetCDF.createProfileVariable(fishfile, "ext_currentU", 2);
                        SaveNetCDF.createProfileVariable(fishfile, "ext_currentV", 2);
                        SaveNetCDF.createScalarVariables(fishfile, o2Names);
                    }
                    else {
                        try {
                            fishfile = NetcdfFileWriteable.openExisting(fishfilePath);
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                    }

                    double[] groupArray = new double[fish.getNGroups()];

                    for (int j = 0; j < fishTmp.length; j++)
                        groupArray[j] = fish.getV(j);
                    SaveNetCDF.saveProfileVariable(fishfile, t, "V", 0, groupArray, true);
                    SaveNetCDF.saveScalarVariable(fishfile, t, "waste", outFlow, false);
                    SaveNetCDF.saveScalarVariable(fishfile, t, "waste_net", outFlow_net, false);
                    
                    double totI = 0, totIngRate = 0;
                    for (int j = 0; j < fishTmp.length; j++) {
                        groupArray[j] = fish.getIngested(j);
                        totI += fish.getN(j) * fishTmp[j][0][0];
                        totIngRate += fish.getN(j) * fish.getIngRate(j);
                    }
                    SaveNetCDF.saveProfileVariable(fishfile, t, "ingested", 0, groupArray, false);


                    for (int j = 0; j < groupArray.length; j++) {
                        groupArray[j] = fish.getAppetite(j);
                    }
                    SaveNetCDF.saveProfileVariable(fishfile, t, "appetite", 0, groupArray, false);

                    SaveNetCDF.saveScalarVariable(fishfile, t, "totIngested", totI, false);

                    SaveNetCDF.saveScalarVariable(fishfile, t, "totIngRate", totIngRate, false);

                    double totalFeed = 0;
                    for (int ii = 0; ii < fc.length; ii++)
                        for (int j = 0; j < fc[0].length; j++)
                            for (int k = 0; k < fc[0][0].length; k++)
                                totalFeed += fc[ii][j][k];
                    SaveNetCDF.saveScalarVariable(fishfile, t, "totFeed", totalFeed, false);

                    double meanFeedDepth = AdvectPellets.getMeanFeedDepth(fc, dz, mask);
                    SaveNetCDF.saveScalarVariable(fishfile, t, "meanFeedDepth", meanFeedDepth, false);

                    double dMeanFeedDepth = (lastMeanFeedDepth > 0) ? (meanFeedDepth - lastMeanFeedDepth) / dt : 0;
                    SaveNetCDF.saveScalarVariable(fishfile, t, "d_meanFeedDepth", dMeanFeedDepth, false);

                    SaveNetCDF.saveScalarVariable(fishfile, t, "rho", rho, false);


                    SaveNetCDF.saveScalarVariable(fishfile, t, "o2ConsumptionRate", o2ConsumptionRate, false);

                    // Save temperature (input value):
                    SaveNetCDF.saveProfileVariable(fishfile, t, "temperature", 2, ambientTemp, false);

                    // Save external O2 (input value):
                    SaveNetCDF.saveProfileVariable(fishfile, t, "ext_O2", 2, ambientValueO2, false);


                    // Save external current speed and direction(input values):
                    double[] currentComp = new double[cageDims[2]];
                    for (int j = 0; j < currentComp.length; j++) {
                        currentComp[j] = currentProfile[j][0];
                    }
                    SaveNetCDF.saveProfileVariable(fishfile, t, "ext_currentU", 2, currentComp, false);
                    for (int j = 0; j < currentComp.length; j++) {
                        currentComp[j] = currentProfile[j][1];
                    }
                    SaveNetCDF.saveProfileVariable(fishfile, t, "ext_currentV", 2, currentComp, false);

                    // Save feeding rage (input value):
                    SaveNetCDF.saveScalarVariable(fishfile, t, "feedingRate", feedingRateMult, false);

                    // Save minimum O2 value:
                    double[] values = minValueMeanAndFracHypoxia(o2, mask);
                    SaveNetCDF.saveScalarVariable(fishfile, t, "min_O2", values[0], false);
                    SaveNetCDF.saveScalarVariable(fishfile, t, "mean_O2", values[1], false);
                    SaveNetCDF.saveScalarVariable(fishfile, t, "frac_hypoxia", values[2], false);

                    // Save o2 level at sensor positions:
                    for (int j=0; j<o2Names.length; j++) {
                        // Store value at sensor position:
                        SaveNetCDF.saveScalarVariable(fishfile, t, o2Names[j], o2[o2Pos[j][0]][o2Pos[j][1]][o2Pos[j][2]], false);
                        /*// Get values in a neighbourhood of the sensor to calculate spatial variability:
                        double nearStd = getStdAround(o2, o2Pos[j]);
                        SaveNetCDF.saveScalarVariable(fishfile, t, o2Names[j]+"_std", nearStd, false);*/
                    }

                    lastMeanFeedDepth = meanFeedDepth;

                    try {
                        fishfile.close();
                    } catch (IOException e) {
                        e.printStackTrace();
                    }

                }


            }


            double totI = 0;
            for (int j = 0; j < fishTmp.length; j++) {
                totI += fish.getN(j) * fish.getIngested(j);
            }

            System.out.println("totFeedAdded = " + totFeedAdded);
            System.out.println("totI = " + totI);
            System.out.println("Feed wastage: " + nf.format(100 * (totFeedAdded - totI) / totFeedAdded) + " %");

        }
    }


    private static void setAllVal(double[][][][] h, double[] v) {
        for (int i=0; i<h.length; i++)
            for (int j=0; j<h[i].length; j++)
                for (int k=0; k<h[i][j].length; k++) {
                    for (int z=0; z<v.length; z++)
                        h[i][j][k][z] = v[z];
                }
    }

    public static double standardDev(double[] data) {
        // The mean average
        double mean = 0.0;
        for (int i=0; i<data.length; i++) {
            mean += data[i];
        }
        mean /= data.length;
        // The variance
        double variance = 0.0;
        for (int i=0; i<data.length; i++) {
            variance += Math.pow(data[i]-mean, 2);
        }
        variance /= (data.length-1);
        // Standard Deviation
        return Math.sqrt(variance);
    }
    private static double[] makeProfile (int n, double val) {
        double[] newProfile = new double[n];
        for (int i = 0; i < newProfile.length; i++) {
            newProfile[i] = val;                    // replacing each element is 0
        }
        return newProfile;
    }

    private static double[] makeDepthProfile ( double val) {
        int n = (int) val;
        double[] newProfile = new double[n];
        for (int i = 0; i < newProfile.length; i++) {
            newProfile[i] = i;                    // replacing each element is 0
        }
        return newProfile;
    }



    private static double[] minValueMeanAndFracHypoxia(double[][][] h, boolean[][][] mask) {
        double minval = Double.MAX_VALUE;
        double meanVal = 0;
        int hypoCells = 0, nonHypoCells = 0, cageCells = 0;
        for (int i=0; i<h.length; i++)
            for (int j=0; j<h[i].length; j++)
                for (int k=0; k<h[i][j].length; k++) {
                    if (mask[i][j][k]) {
                        cageCells++;
                        meanVal += h[i][j][k];
                        if (h[i][j][k] < minval)
                            minval = h[i][j][k];
                        if (h[i][j][k] < HYPOXIA_THRESHOLD)
                            hypoCells++;
                        else
                            nonHypoCells++;
                    }
                }
        return new double[] {minval, meanVal/((double)cageCells),
                ((double)hypoCells)/((double)(nonHypoCells+hypoCells))};
    }

    private static int countAvailableCellsForOxygenUptake(int[] cageDims, double dz, double fishMaxDepth,
                                                          boolean[][][] mask) {
        int res = 0;
        for (int i=0; i<cageDims[0]; i++)
            for (int j=0; j<cageDims[1]; j++)
                for (int k=0; k<cageDims[2]; k++) {
                    double lDepth = ((double) k + 0.5) * dz;
                    if ((mask == null || mask[i][j][k]) && (lDepth < fishMaxDepth)) {
                        res++;
                    }
                }

        return res;
    }

    private static double setO2AffinityWithVerticalProfileAndEdgeDecrease(int[] cageDims, double rad, double dxy, double dz, double fishMaxDepth, boolean[][][] mask,
                                                           double[] affinityProfile, double[][][] affinity,
                                                           double[][][] o2Affinity) {
        double o2AffSum = 0;


        for (int i=0; i<cageDims[0]; i++)
            for (int j=0; j<cageDims[1]; j++) {

                double xDist = dxy*(double)(i - cageDims[0]/2);
                double yDist = dxy*(double)(j - cageDims[1]/2);
                double distFromCenter = Math.sqrt(xDist*xDist + yDist*yDist);

                for (int k = 0; k < cageDims[2]; k++) {
                    double lDepth = ((double) k + 0.5) * dz;
                    if ((mask == null || mask[i][j][k]) && (lDepth < fishMaxDepth)) {
                        // Check oxygen avoidance criterion:
                        affinity[i][j][k] = affinityProfile[k];
                        o2Affinity[i][j][k] = affinityProfile[k];

                    } else {
                        affinity[i][j][k] = 0;
                        o2Affinity[i][j][k] = 0;

                    }

                    // TEST TEST TEST
                    // Decrease affinity towards edge:
                    double edgeFactor = 1.;
                    if (distFromCenter/rad > 0.5) {
                        edgeFactor -= Math.pow((distFromCenter-0.5*rad)/(0.5*rad), 4);
                        affinity[i][j][k] *= edgeFactor;
                        o2Affinity[i][j][k] *= edgeFactor;
                    }

                    o2AffSum += o2Affinity[i][j][k];
                }
            }
        return o2AffSum;
    }

    private static double setO2AffinityWithVerticalProfile(int[] cageDims, double dz, double fishMaxDepth, boolean[][][] mask,
                                                           double[] affinityProfile, double[][][] affinity,
                                                           double[][][] o2Affinity) {
        double o2AffSum = 0;


        for (int i=0; i<cageDims[0]; i++)
            for (int j=0; j<cageDims[1]; j++)
                for (int k = 0; k < cageDims[2]; k++) {
                    double lDepth = ((double) k + 0.5) * dz;
                    if ((mask == null || mask[i][j][k]) && (lDepth < fishMaxDepth)) {
                        // Check oxygen avoidance criterion:
                        affinity[i][j][k] = affinityProfile[k];
                        o2Affinity[i][j][k] = affinityProfile[k];
                        //affinity[i][j][k] = 1.;//affinityProfile[k];
                        //o2Affinity[i][j][k] = 1.;//affinityProfile[k];

                    } else {
                        affinity[i][j][k] = 0;
                        o2Affinity[i][j][k] = 0;

                    }

                    o2AffSum += o2Affinity[i][j][k];
                }

        return o2AffSum;
    }

    private static double setO2AffinityWithVerticalProfileAndDirection(int[] cageDims, double dz, double fishMaxDepth, double[] dirVector, boolean[][][] mask,
                                                           double[] affinityProfile, double[][][] affinity,
                                                           double[][][] o2Affinity) {
        double o2AffSum = 0;
        double centerX = ((double)cageDims[0])/2.0,
                centerY = ((double)cageDims[1])/2.0;
        double dirVectorLength = Math.sqrt(dirVector[0]*dirVector[0] + dirVector[1]*dirVector[1]);

        for (int i=0; i<cageDims[0]; i++)
            for (int j=0; j<cageDims[1]; j++) {
                // Set up unit vector in the direction from center to this horizontal position:
                double distX = ((double)i) - centerX,
                        distY = ((double)j) - centerY;
                double distTot = Math.sqrt(distX*distX+distY*distY);
                if (distTot > 0) {
                    distX = distX/distTot;
                    distY = distY/distTot;
                }
                // Calculate dot product between direction vector and the input dir vector. This gives a value
                // that is equal to the length of the input dir vector multiplied by the cosine of the angle between
                // the two vectors:
                double dotProduct = distX*dirVector[0] + distY*dirVector[1];
                //dotProduct = (dotProduct+1)*(dotProduct+1)-1;
                System.out.println("distX="+distX+", distY="+distY+", dotProd="+dotProduct);


                for (int k = 0; k < cageDims[2]; k++) {
                    double lDepth = ((double) k + 0.5) * dz;
                    if ((mask == null || mask[i][j][k]) && (lDepth < fishMaxDepth)) {
                        // Check oxygen avoidance criterion:
                        //affinity[i][j][k] = Math.max(0., affinityProfile[k] + dotProduct);
                        //o2Affinity[i][j][k] = Math.max(0., affinityProfile[k] + dotProduct);
                        if (dotProduct/dirVectorLength > 0*0.2588) {
                            affinity[i][j][k] = affinityProfile[k];
                            o2Affinity[i][j][k] = affinityProfile[k];
                        } else {
                            affinity[i][j][k] = affinityProfile[k];
                            o2Affinity[i][j][k] = 0;
                        }
                        //affinity[i][j][k] = affinityProfile[k];
                        //o2Affinity[i][j][k] = affinityProfile[k];

                    } else {
                        affinity[i][j][k] = 0;
                        o2Affinity[i][j][k] = 0;

                    }

                    o2AffSum += o2Affinity[i][j][k];
                }
            }
        return o2AffSum;
    }

    private static double setO2AffinityWithAvoidance(int[] cageDims, double dz, double fishMaxDepth, boolean[][][] mask, double[][][] o2,
                                                     double[][][] affinity, double[][][] o2Affinity) {
        double o2AffSum = 0;

        double surfaceAvoidanceDepth = 4; // Range of (gradual) surface avoidance, m
        double avoidanceThresh = 4;

        for (int i=0; i<cageDims[0]; i++)
            for (int j=0; j<cageDims[1]; j++)
                for (int k=0; k<cageDims[2]; k++) {
                    double lDepth = ((double)k+0.5)*dz;
                    if ((mask==null || mask[i][j][k]) && (lDepth < fishMaxDepth)) {
                        // Check oxygen avoidance criterion:
                        if (o2[i][j][k] > avoidanceThresh) {
                            affinity[i][j][k] = 1;
                            o2Affinity[i][j][k] = 1;
                        }
                        else {
                            affinity[i][j][k] = o2[i][j][k]/avoidanceThresh; // Linear decrease below threshold
                            o2Affinity[i][j][k] = o2[i][j][k]/avoidanceThresh; // Linear decrease below threshold
                        }

                        // Check surface avoidance criterion:
                        /*if (lDepth < surfaceAvoidanceDepth) {
                            affinity[i][j][k] *= lDepth/surfaceAvoidanceDepth;
                            o2Affinity[i][j][k] *= lDepth/surfaceAvoidanceDepth;
                        }*/

                    }
                    o2AffSum += o2Affinity[i][j][k];
                }
        return o2AffSum;
    }


    private static double getStdAround(double[][][] field, int[] pos) {
        double[] values = new double[9];
        values[0] = field[pos[0]-1][pos[1]-1][pos[2]];
        values[1] = field[pos[0]-1][pos[1]][pos[2]];
        values[2] = field[pos[0]-1][pos[1]+1][pos[2]];
        values[3] = field[pos[0]][pos[1]-1][pos[2]];
        values[4] = field[pos[0]][pos[1]][pos[2]];
        values[5] = field[pos[0]][pos[1]+1][pos[2]];
        values[6] = field[pos[0]+1][pos[1]-1][pos[2]];
        values[7] = field[pos[0]+1][pos[1]][pos[2]];
        values[8] = field[pos[0]+1][pos[1]+1][pos[2]];

        return standardDev(values);
    }

    /**
     * Set up vertical profile for model grid based on values at set depths. Extrapolate beyond end values
     * Interpolates using linear interpolator from Apache Commons Math.
     * @param res The array to put interpolated values into
     * @param depths Depths at which values are given (increasing values)
     * @param values Values at given depths.
     * @param kmax Number of vertical layers
     * @param dz Vertical resolution
     * @return Interpolated/extrapolated profile
     */
    public static void interpolateVertical(double[] res, double[] depths, double[] values, int kmax, double dz) {

        double minDepth = depths[0], maxDepth = depths[depths.length-1],
                topValue = values[0], bottomValue = values[values.length-1];
        PolynomialSplineFunction interp = interpol.interpolate(depths, values);
        for (int i=0; i<res.length; i++) {
            double currDepth = ((double)i + 0.5)*dz;
            if (currDepth < minDepth)
                res[i] = topValue;
            else if (currDepth > maxDepth)
                res[i] = bottomValue;
            else
                res[i] = interp.value(currDepth);
        }


    }
}
