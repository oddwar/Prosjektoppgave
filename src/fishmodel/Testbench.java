/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package fishmodel;
import java.time.LocalDate;
import java.time.Month;
import java.time.DayOfWeek;
import fishmodel.CageMasking;
import fishmodel.Measurements;
import fishmodel.enkf.AssimSettings;
import fishmodel.enkf.EnsembleKF;
import fishmodel.enkf.MpiHandler;
import fishmodel.enkf.Util;
import fishmodel.hydraulics.Balanced3DHydraulics;
import fishmodel.hydraulics.SimpleTankHydraulics;
import fishmodel.pellets.AdvectPellets;
import fishmodel.pellets.IngestionAndO2Tempprofile;
import fishmodel.pellets.PelletSpreaderModel;
import fishmodel.pellets.SimpleFish;
import fishmodel.sim.CurrentMagicFields;
import fishmodel.sim.InputDataNetcdf;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import save.SaveNetCDF;
import ucar.nc2.NetcdfFileWriteable;

import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Date;
import java.util.Locale;
import java.util.Random;


/**
 *
 * @author Oddwar
 * This file is made in order to test part of the code and different theories
 * first up is the currentprofile
 */

import fishmodel.hydraulics.SimpleTankHydraulics;

public class Testbench {
    public static LinearInterpolator interpol = new LinearInterpolator();

    public static void main(String[] args) {



        double[][][][] hydro;
        double rad = 25;
        double depth = 40;
        double topOfCage = 0;
        double cylDepth = 20;
        double totDepth = 35; // Cage size (m)
        double dxy = 2, dz = dxy; // Model resolution (m) originally set to 2
        double dt = 0.5 * dxy; // Time step (s) originally set to 0.5
        int storeIntervalFeed = 60, storeIntervalInfo = 60;
        double fishMaxDepth = 35; // The maximum depth of the fish under non-feeding conditions
        double currentSpeedInit = 1.0;
        // Current profile from the Aquaexcel project

        double[] low_current = {0.041 ,  0.03 ,  0.021 ,     0.015  ,  0.0145 ,   0.013  ,   0.012};
        double[] current_depths = { 2.5 , 5 , 7.5   ,  12.5 , 17.5  , 27.5 ,  32.5};



        double modelDim = 2*(rad+8*dxy);
        int[] cageDims = new int[3];
        cageDims[0] = (int)Math.ceil(modelDim/dxy);
        cageDims[1] = cageDims[0];
        cageDims[2] = (int)Math.ceil(depth/dz)+1;

        double[][] currentProfile = new double[cageDims[2]+1][3];
        for (int i=0; i<cageDims[2]+1; i++) {

            //updating currentprofile from 0 to CurrentSpeedInit
            currentProfile[i][0] = currentSpeedInit;
        }


        double[] tempCurrentProfile = new double[cageDims[2]+1];
        interpolateVertical(tempCurrentProfile, current_depths, low_current, cageDims[2],dz);

        double [][] newCurrentProfile = new double[cageDims[2]+1][3];

        for (int i = 0; i< tempCurrentProfile.length;i++) {
            newCurrentProfile[i][0] = tempCurrentProfile[i];

        }




        hydro = SimpleTankHydraulics.getProfileHydraulicField(cageDims, newCurrentProfile);


        //System.out.print(newCurrentProfile[0]);
        for (int i=0; i<cageDims[2]+1; i++) {

            //updating currentprofile from 0 to CurrentSpeedInit
            System.out.println(newCurrentProfile[i][0]) ;
        }

        double per = 10.0;
        double paal = 2 ;
        System.out.println("lengde hydro: ");
        //System.out.println(per/paal);

        System.out.println(hydro.length);
        System.out.println(hydro[0].length);
        System.out.println(hydro[0][0].length);


        for (int i=0; i<hydro[0][3].length; i++) {

            //System.out.println(hydro[i][i][0][0]) ;
            //System.out.println(hydro[0][i][0][0]) ;
            //System.out.println(hydro[0][0][i][0]) ;//current
            //System.out.println(hydro[0][0][0][i]) ;
        }




        //** Alternating the current after 6 hours = 21 600 sec



    }

    /**







     // Current profile from the Aquaexcel project
     //Low_current =  [   4.1 ,   3 ,  2.1 ,     1.5  ,  1.45 ,   1.3  ,   1.2]
     //current_depths = [- 2.5 , -5 , -7.5   ,  -12.5 , -17.5  , -27.5 ,  -32.5]
*/
    // interpolateVertical(affinityProfile, affDepths, affProfile, cageDims[2], dz);

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





    /**
     * Setup for testing, remember to change filename
     * test 3 - submerged, all values outside grid = NaN, maximizing number of fish, no current,

    public static void main(String[] args) {
        LocalDate currentDate = LocalDate.now();
        int year = currentDate.getYear();
        Month month = currentDate.getMonth();
        DayOfWeek day = currentDate.getDayOfWeek();
        System.out.println("Årstall : " + year);
        System.out.println("måned : " + month.toString());
        System.out.println("dag : " + day.toString());

    }
    */
}


