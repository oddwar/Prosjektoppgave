/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package fishmodel;

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
public class Testbench {

    public static final double HYPOXIA_THRESHOLD = 6;

    public static LinearInterpolator interpol = new LinearInterpolator();

    /**
     * Setup for testing, remember to change filename
     * test 3 - submerged, all values outside grid = NaN, maximizing number of fish, no current,
     */
    public static void main(String[] args) {

        String inDataFile = "C:/Users/stein/Documents/01 Industriel kybernetikk/5. semester/Prosjektoppgave/Filer 03 aug/bjoroya_data.nc";
        boolean useInstantaneousAmbientVals = true; // true to use ambient value time series, false to use daily averages

        int initYear = 2022, initMonth = Calendar.JUNE, initDate = 22, initHour = 10, initMin = 0, initSec = 10;
        Calendar c = Calendar.getInstance();

        c.set(initYear, initMonth, initDate, initHour, initMin, initSec);
        //c.add(Calendar.DAY_OF_MONTH, 1);
        Date startTime = c.getTime();

        //------------------
        InputDataNetcdf inData = new InputDataNetcdf(inDataFile, useInstantaneousAmbientVals);
        inData.setStartTime(startTime);
        double number = 5;
        int multiplier = 10;

        double[] product = new double[multiplier]; // array with 10 elements

        for (int i = 0; i < product.length; i++) {
            product[i] = number;                    // replacing each element is 0
        }

        double[] obsCurrentDepths = inData.getCurrentDepths();

        double[] obsCurrentProfile = inData.getExtCurrentSpeedProfile();
        //System.out.println("current profile: ");
        System.out.println("value: \n");
        for (int i = 0; i < obsCurrentDepths.length; i++) {
            System.out.print(obsCurrentDepths[i]+".\n" );
        }



    }
}