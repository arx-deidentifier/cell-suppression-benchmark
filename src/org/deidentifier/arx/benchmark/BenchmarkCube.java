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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.deidentifier.arx.ARXAnonymizer;
import org.deidentifier.arx.ARXConfiguration;
import org.deidentifier.arx.ARXResult;
import org.deidentifier.arx.AttributeType;
import org.deidentifier.arx.AttributeType.Hierarchy;
import org.deidentifier.arx.AttributeType.Hierarchy.DefaultHierarchy;
import org.deidentifier.arx.Data;
import org.deidentifier.arx.DataHandle;
import org.deidentifier.arx.DataType;
import org.deidentifier.arx.aggregates.StatisticsFrequencyDistribution;
import org.deidentifier.arx.criteria.DistinctLDiversity;
import org.deidentifier.arx.criteria.KAnonymity;
import org.deidentifier.arx.exceptions.RollbackRequiredException;
import org.deidentifier.arx.metric.Metric;


/**
 * This class encapsulates experiments determining the R-U frontiers for different
 * combination of risk measures and interpretations of missing values.
 * 
 * @author Fabian Prasser
 * @author Helmut Spengler
 * @author Johanna Eicher
 * @author Raffael Bild
 *
 */
public class BenchmarkCube {

    /** The directory containing the input files */
    private static String DATA_DIR = "./data/";

    /**
     * Main entry point
     * @param args
     * @throws IOException
     * @throws RollbackRequiredException
     */
    public static void main(String[] args) throws IOException, RollbackRequiredException {
        
        // No interpolation
        // - LM: 0.09761653690984728
        // - Median relative error (point queries): 0.0
        // - Median relative error (range queries): 0.7023809523809523
        // Simple interpolation
        // - LM: 0.09761653690984728
        // - Final results
        // - Median relative error (point queries): 0.0
        // - Median relative error (range queries): 0.6510762075962129
        // Estimated using posteriori distribution with sigmoid activation function
        // - Final results
        // - Median relative error (point queries): 0.0
        // - Median relative error (range queries): 0.36184886222834983

        // Perform experiment
        Data data = getData("adult");
        DataHandle output = getOutput(data);
        DataHandle input = data.getHandle();

        // Obtain attributes
        Set<String> attributes = new HashSet<>();
        for (int i = 0; i < 8; i++) {
            attributes.add(input.getAttributeName(i));
        }

        // Statistics
        DescriptiveStatistics statsPointQuery = new DescriptiveStatistics();
        DescriptiveStatistics statsRangeQuery = new DescriptiveStatistics();
        DescriptiveStatistics statsLM = new DescriptiveStatistics();
        
        // Posteriori likelihoods
        Map<Integer, Map<String, Double>> likelihoods = getLikelihoods(output);
        
        // For each set in the power set
        for (Set<String> set : getPowerSet(attributes)) {
            
            // Skip the empty set
            if (set.isEmpty()) {
                continue;
            }

            // Information loss
            statsLM.addValue((1d - output.getStatistics().getQualityStatistics(set).getGranularity().getArithmeticMean()));
           
            // Perform experiments
            for (int i = 0; i < 1000; i++) {

                // Create
                Map<Integer, Set<String>> pointQuery = getRandomPointQuery(input, set, true);
                Map<Integer, Set<String>> rangeQuery = getRandomRangeQuery(input, set);

                // Perform and analze
                statsPointQuery.addValue(getRelativeError(getCount(pointQuery, input, null), getCount(pointQuery, output, null)));
                statsRangeQuery.addValue(getRelativeError(getCount(rangeQuery, input, null), getCount(rangeQuery, output, likelihoods)));
            }
        }

        // Print statistics
        System.out.println(" - Final results");
        System.out.println(" - LM: " + statsLM.getMean());
        System.out.println(" - Median relative error (point queries): " + statsPointQuery.getPercentile(50d));
        System.out.println(" - Median relative error (range queries): " + statsRangeQuery.getPercentile(50d));

    }
    /**
     * Executes the query on the given handle
     * @param query
     * @param handle
     * @param likelihoods
     * @return
     */
    private static double getCount(Map<Integer, Set<String>> query, DataHandle handle, Map<Integer, Map<String, Double>> likelihoods) {
        
        // Prepare
        double count = 0;
        Map<Integer, Double> preparedLikelihoods = new HashMap<>();
        if (likelihoods != null) {
            for (Entry<Integer, Set<String>> predicate : query.entrySet()) {
                double likelihood = 0d;
                for (String value : predicate.getValue()) {
                    Double frequency = likelihoods.get(predicate.getKey()).get(value);
                    likelihood += (frequency == null ? 0d : frequency);
                }
                preparedLikelihoods.put(predicate.getKey(), likelihood);
            }
        }
        
        // Search all rows
        outer: for (int row = 0; row < handle.getNumRows(); row++) {
            
            double likelihood = 1d;
            
            // For each predicate
            for (Entry<Integer, Set<String>> predicate : query.entrySet()) {
                
                // Value
                String value = handle.getValue(row, predicate.getKey());
                
                // If not found
                if (!predicate.getValue().contains(value)) {
                    
                    // Skip if no likelihoods have been defined
                    if (likelihoods == null || !value.equals("*")) {
                        continue outer;
                    } else {
                        likelihood *= preparedLikelihoods.get(predicate.getKey());
                    }
                }
            }

            // Found
            if (likelihoods == null) {
                count++;
            } else {
                count += (1 / (1 + Math.exp(-likelihood * 5)) - 0.5) * 2;
            }
        }
        
        // Done
        return count;
    }
    
    /**
     * Load and configure a dataset with the given number of quasi-identifiers
     * 
     * @param dataset
     * @return
     * @throws IOException 
     */
    private static Data getData(String dataset) throws IOException {        
       
        Data data = Data.create(DATA_DIR + dataset + ".csv", StandardCharsets.UTF_8, ';');

        // Declare all attributes as insensitive
        for (int i = 0; i < data.getHandle().getNumColumns(); i++) {
            data.getDefinition().setAttributeType(data.getHandle().getAttributeName(i), AttributeType.INSENSITIVE_ATTRIBUTE);
        }

        for (int i = 0; i < 7; i++) {
            String attribute = data.getHandle().getAttributeName(i);
            data.getDefinition().setAttributeType(attribute, getHierarchy(data, attribute));
        }      
        data.getDefinition().setAttributeType("occupation", AttributeType.SENSITIVE_ATTRIBUTE);
        return data;
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
     * Creates a map of likelihoods
     * @param output
     * @return
     */
    private static Map<Integer, Map<String, Double>> getLikelihoods(DataHandle output) {
        Map<Integer, Map<String, Double>> likelihoods = new HashMap<>();
        for (int i=0; i<output.getNumColumns(); i++) {
            likelihoods.put(i, new HashMap<String, Double>());
            StatisticsFrequencyDistribution distribution = output.getStatistics().getFrequencyDistribution(i);
            for (int j=0; j<distribution.values.length; j++) {
                likelihoods.get(i).put(distribution.values[j], distribution.frequency[j]);
            }
        }
        return likelihoods;
    }  

    /**
     * Perform cell suppression.
     * 
     * @param data
     * @param numQis
     * @return
     * @throws IOException
     * @throws RollbackRequiredException
     */
    private static DataHandle getOutput(Data data) throws IOException, RollbackRequiredException {
        
        double o_min = 0.001d;
        ARXConfiguration config = ARXConfiguration.create();

        double maxOutliers = 1.0d - o_min;
        config.setSuppressionLimit(maxOutliers);
        config.setQualityModel(Metric.createLossMetric(0d));
        config.addPrivacyModel(new KAnonymity(5));
        config.addPrivacyModel(new DistinctLDiversity("occupation", 3));
        config.setHeuristicSearchEnabled(false);
        
        // Configure and set up   
        ARXAnonymizer anonymizer = new ARXAnonymizer();
        ARXResult result = anonymizer.anonymize(data, config);
        DataHandle output = result.getOutput();
        if (output != null && result.isOptimizable(output)) {
            result.optimizeIterativeFast(output, o_min);
        }    
        return output;
    }

    /**
     * From: https://stackoverflow.com/questions/1670862/obtaining-a-powerset-of-a-set-in-java
     * @param originalSet
     * @return
     */
    private static <T> Set<Set<T>> getPowerSet(Set<T> originalSet) {
        Set<Set<T>> sets = new HashSet<Set<T>>();
        if (originalSet.isEmpty()) {
            sets.add(new HashSet<T>());
            return sets;
        }
        List<T> list = new ArrayList<T>(originalSet);
        T head = list.get(0);
        Set<T> rest = new HashSet<T>(list.subList(1, list.size())); 
        for (Set<T> set : getPowerSet(rest)) {
            Set<T> newSet = new HashSet<T>();
            newSet.add(head);
            newSet.addAll(set);
            sets.add(newSet);
            sets.add(set);
        }       
        return sets;
    }

    /**
     * Creates a random query
     * @param input
     * @param attributes
     * @param draw
     * @return
     */
    private static Map<Integer, Set<String>> getRandomPointQuery(DataHandle input, Set<String> attributes, boolean draw) {

        // Prepare
        Map<Integer, Set<String>> query = new HashMap<>();
        
        // Select record
        int record = (int)Math.round(Math.random() * input.getNumRows());
        
        // Construct query
        for (String attribute : attributes) {
            
            int column = input.getColumnIndexOf(attribute);
            
            if (draw) {
                List<String> values = Arrays.asList(input.getDistinctValues(column));
                Collections.shuffle(values);
                query.put(column, new HashSet<String>());
                query.get(column).add(values.get(0));
            } else {
                query.put(column, new HashSet<String>());
                query.get(column).add(input.getValue(record, column));
            }
        }
        
        // Return
        return query;
    }

    /**
     * Creates a random query
     * @param input
     * @param attributes
     * @return
     */
    private static Map<Integer, Set<String>> getRandomRangeQuery(DataHandle input, Set<String> attributes) {

        // Prepare
        Map<Integer, Set<String>> query = new HashMap<>();
        
        // Construct query
        for (String attribute : attributes) {
            
            // Extract
            int column = input.getColumnIndexOf(attribute);
            
            // Prepare
            String[] values = input.getDistinctValues(column);
            int lower = (int)Math.round(Math.random() * (values.length - 1));
            int upper = (int)Math.round(Math.random() * (values.length - 1));
            if (lower > upper) {
                int temp = lower;
                lower = upper;
                upper = temp;
            }
            
            // Build
            query.put(column, new HashSet<String>());
            for (int i = lower; i <= upper; i++) {
                query.get(column).add(values[i]);
            }
        }
        
        // Return
        return query;
    }

    /**
     * Calculate the relative error
     * @param countInput
     * @param countOutput
     * @return
     */
    private static double getRelativeError(double countInput, double countOutput) {
        if (countInput == 0 && countOutput == 0) {
            return 0;
        } else if (countInput == 0) {
            return Double.POSITIVE_INFINITY;
        } else {
            return Math.abs(countInput - countOutput) / countInput;
        }
    }
}
