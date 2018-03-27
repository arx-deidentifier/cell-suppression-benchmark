/*
 * Benchmark of cell-suppression with ARX
 * Copyright (C) 2018 TUM/MRI
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package org.deidentifier.arx.benchmark;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Iterator;

import org.deidentifier.arx.ARXAnonymizer;
import org.deidentifier.arx.ARXConfiguration;
import org.deidentifier.arx.ARXResult;
import org.deidentifier.arx.AttributeType;
import org.deidentifier.arx.AttributeType.Hierarchy;
import org.deidentifier.arx.AttributeType.Hierarchy.DefaultHierarchy;
import org.deidentifier.arx.Data;
import org.deidentifier.arx.DataHandle;
import org.deidentifier.arx.DataType;
import org.deidentifier.arx.criteria.AverageReidentificationRisk;
import org.deidentifier.arx.criteria.KAnonymity;
import org.deidentifier.arx.exceptions.RollbackRequiredException;
import org.deidentifier.arx.metric.Metric;


/**
 * This class encapsulates experiments measuring the influence of using different numbers of
 * quasi-identifiers on execution times and on data quality.
 * 
 * @author Helmut Spengler
 * @author Fabian Prasser
 * @author Johanna Eicher
 * @author Raffael Bild
 *
 */
public class Benchmark {

    /**
     * This class encapsulates parameters related to risk management. It can be used
     * for threshold definitions as well as for storing/retrieving results.
     * 
     * @author Fabian Prasser
     * @author Helmut Spengler
     * @author Johanna Eicher
     * @author Raffael Bild
     */
    private static class Risks {
        /** Threshold for the average risk */
        final private double averageRisk;;
        /** Threshold for the highest risk */
        final private double highestRisk;;
        /** Threshold for records at risk. */
        final private double recordsAtRisk;

        /**
         * Creates a new instance
         * @param averageRisk
         * @param highestRisk
         * @param recordsAtRisk
         * @param qis
         */
        public Risks(double averageRisk, double highestRisk, double recordsAtRisk) {
            this.averageRisk   = averageRisk;
            this.highestRisk   = highestRisk;
            this.recordsAtRisk = recordsAtRisk;
        }

        @Override
        public String toString() {
            return "Risks [averageRisk=" + averageRisk + ", highestRisk=" + highestRisk + ", recordsAtRisk=" + recordsAtRisk + "]";
        }
    }

    /** The directory containing the input files */
    private static String DATA_DIR = "./data/";

    /**
     * Main entry point
     * @param args
     * @throws IOException
     * @throws RollbackRequiredException
     */
    public static void main(String[] args) throws IOException, RollbackRequiredException {        
        // For each file
        for (String file : new String[]{"adult", "ihis"}) {

            // for each records at risk variant
            for (Risks risks : new Risks[] {new Risks(1.0d, 0.2d, 0.0d),
                                            new Risks(1.0d, 0.2d, 0.1d)}) {
                // Suppressed cells
                System.out.println(file);
                System.out.println("Suppressed cells (privacy = " + risks.toString()+")");
                for (int i = 1; i <= 9; i++ ) {
                    evaluateSuppressedCells(file, risks, i);
                }
                // Non-uniform entropy
                System.out.println("Non-uniform entropy (privacy = " + risks.toString()+")");
                for (int i = 1; i <= 9; i++ ) {
                    evaluateEntropy(file, risks, i);
                }
                // Execution time
                System.out.println("Execution time [ms] (privacy = " + risks.toString()+")");
                for (int i = 1; i <= 9; i++ ) {
                    evaluateExecutionTime(file, risks, i);
                } 
            }
        }
    }

    /**
     * Perform anonymization and measure the fraction of suppressed cells.
     * 
     * @param dataset
     * @param qis
     * @param risks
     * @throws IOException
     * @throws RollbackRequiredException
     */
    private static void evaluateSuppressedCells(String dataset, Risks risks, int numQis) throws IOException, RollbackRequiredException {
        DataHandle result = performCellSuppression(dataset, risks, numQis);
        System.out.println(numQis + " - " + String.valueOf(countSuppressedCells(result.iterator())));
    }

    /**
     * Perform anonymization and measure Non-Uniform Entropy.
     * 
     * @param dataset
     * @param qis
     * @param risks
     * @throws IOException
     * @throws RollbackRequiredException
     */
    private static void evaluateEntropy(String dataset, Risks risks, int numQis) throws IOException, RollbackRequiredException {
        DataHandle result = performCellSuppression(dataset, risks, numQis);
        System.out.println(numQis + " - " + result.getStatistics().getQualityStatistics().getNonUniformEntropy().getArithmeticMean(false));
    }

    /**
     * Perform anonymization and measure execution times.
     * 
     * @throws IOException
     * @throws RollbackRequiredException 
     */
    private static void evaluateExecutionTime(String dataset, Risks risks, int numQis) throws IOException, RollbackRequiredException {
        long start = System.currentTimeMillis();     
        performCellSuppression(dataset, risks, numQis);
        System.out.println(numQis + " - " + (System.currentTimeMillis() - start));
    }

    /**
     * Load and configure a dataset with the given number of quasi-identifiers
     * 
     * @param dataset
     * @param risks
     * @param numQis
     * @return
     * @throws IOException 
     */
    private static Data getData(String dataset, int numQis) throws IOException {        
        String[] qis = Arrays.copyOfRange(getQis(dataset), 0, numQis);

        Data data = Data.create(DATA_DIR + dataset + ".csv", StandardCharsets.UTF_8, ';');

        int numColumns = data.getHandle().getNumColumns();
        // Declare all attributes as insensitive
        for (int i = 0; i < numColumns; i++) {
            data.getDefinition().setAttributeType(data.getHandle().getAttributeName(i), AttributeType.INSENSITIVE_ATTRIBUTE);
        }

        for (int i = 0; i < qis.length; i++) {
            String attribute = data.getHandle().getAttributeName(i);
            data.getDefinition().setAttributeType(attribute, getHierarchy(data, attribute));
        }      
        return data;
    }

    /**
     * Create an ARX configuration with the parameters used in the rest of the paper.
     * 
     * @param risks
     * @return
     */
    private static ARXConfiguration createARXConfig(Risks risks) {
        ARXConfiguration config = ARXConfiguration.create();

        double maxOutliers = 1.0d - 0.01d;
        config.setSuppressionLimit(maxOutliers);
        config.setQualityModel(Metric.createLossMetric(0d));
        if (risks.recordsAtRisk == 0d) {
            if (risks.averageRisk != 1d) {
                config.addPrivacyModel(new AverageReidentificationRisk(risks.averageRisk));
            }
            int k = getSizeThreshold(risks.highestRisk);
            if (k != 1) {
                config.addPrivacyModel(new KAnonymity(k));
            }
        } else {
            config.addPrivacyModel(new AverageReidentificationRisk(risks.averageRisk, risks.highestRisk, risks.recordsAtRisk));   
        }
        config.setHeuristicSearchEnabled(false);
        return config;
    }

    /**
     * Perform cell suppression.
     * 
     * @param arxData
     * @param config
     * @return
     * @throws IOException
     * @throws RollbackRequiredException
     */
    private static DataHandle performCellSuppression(String dataset, Risks risks, int numQis) throws IOException, RollbackRequiredException {
        // Configure and set up
        Data data = getData(dataset, numQis);
        ARXConfiguration config = createARXConfig(risks);        
        ARXAnonymizer anonymizer = new ARXAnonymizer();

        // Generalize
        ARXResult result = anonymizer.anonymize(data, config);

        // Suppress cells
        DataHandle output = result.getOutput();
        if (output != null && result.isOptimizable(output)) {
            result.optimizeIterativeFast(output, 0.01d);
        }    
        return output;
    }

    /**
     * Return the fraction of the suppressed cells in a dataset.
     * 
     * @param result
     * @return
     */
    private static double countSuppressedCells(Iterator<String[]> iterator) {
        double suppressedCells = 0;
        double numCells = 0;

        while (iterator.hasNext()) {
            String[] line = iterator.next();
            for (int i = 0; i < line.length; i++) {
                numCells++;
                if (DataType.ANY_VALUE.equals(line[i])) {
                    suppressedCells++;
                }
            }
        }
        return suppressedCells/numCells;
    }

    /**
     * Return all possible QIs for a specific dataset. The order of the QIs is the same as in the rest of the paper.
     * 
     * @param dataset
     * @return
     */
    private static String[] getQis(String dataset) {
        String[] qisAdult = { "sex", "age", "race", "marital-status", "education", "native-country", "workclass", "occupation", "salary-class" };
        String[] qisIhis  = { "YEAR", "QUARTER", "REGION", "PERNUM", "AGE", "MARSTAT", "SEX", "RACEA", "EDUC" };

        String[] qis;
        if (dataset.equals("adult")) {
            qis = qisAdult;
        } else if (dataset.equals("ihis")) {
            qis = qisIhis;
        } else {
            throw new IllegalArgumentException("invalid dataset: " + dataset);
        }
        return qis;
    }

    /**
     * Returns the generalization hierarchy for the dataset and attribute
     * 
     * @param data
     * @param attribute
     * @return
     * @throws IOException
     */
    private static Hierarchy getHierarchy(Data data, String attribute) {
        DefaultHierarchy hierarchy = Hierarchy.create();
        int col = data.getHandle().getColumnIndexOf(attribute);
        String[] values = data.getHandle().getDistinctValues(col);
        for (String value : values) {
            hierarchy.add(value, DataType.ANY_VALUE);
        }
        return hierarchy;
    }

    /**
     * Returns a minimal class size for the given risk threshold.
     * 
     * @param threshold
     * @return
     */
    private static int getSizeThreshold(double riskThreshold) {
        double size = 1d / riskThreshold;
        double floor = Math.floor(size);
        if ((1d / floor) - (1d / size) >= 0.01d * riskThreshold) {
            floor += 1d;
        }
        return (int) floor;
    }
}
