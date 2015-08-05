/*
 * Copyright 2014-2015 ScAi, CSD, UCLA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package edu.ucla.cs.scai.clustering.syntheticgenerator;

import static edu.ucla.cs.scai.clustering.syntheticgenerator.MultidimensionalGaussianGenerator.createImage;
import static edu.ucla.cs.scai.clustering.syntheticgenerator.MultidimensionalGaussianGenerator.shuffleDataset;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Random;

/**
 *
 * @author Giuseppe M. Mazzeo <mazzeo@cs.ucla.edu>
 */
public class GenerateGaussian {

    public static void main(String args[]) throws FileNotFoundException, IOException, Exception {

        //change with your own paths
        String basePath = "/home/massimo/clubs/experiments/gaussian/data/";
        String basePath2 = "/home/massimo/clubs/experiments/gaussian/shuffled_data/";

        Random rand = new Random(100);
        int domainWidth = 2000;
        int[] dataSetSizes = {100000, 1000000, 10000000};
        int[] nClusters = {0, 1, 2, 4, 8, 16, 32};
        int[] nDimensions = {2, 4, 6, 8, 10, 12};
        double[] noiseRatios = {0.1, 0.08, 0.06, 0.04, 0.02, 0};
        double[] zeroClustersNoiseRatios = {1};
        int maxDimension = nDimensions[nDimensions.length - 1];

        HashMap<Integer, ArrayList<Double>[]> clusterPositions = new HashMap<>();
        HashMap<Integer, double[][]> clusterRadii = new HashMap<>();

        for (int nTuples : dataSetSizes) {
            for (int nOfClusters : nClusters) {
                if (!clusterPositions.containsKey(nOfClusters)) {
                    ArrayList<Double>[] positions = new ArrayList[maxDimension];
                    clusterPositions.put(nOfClusters, positions);
                    double interval = 1.0 / (nOfClusters + 1);
                    for (int dim = 0; dim < positions.length; dim++) {
                        positions[dim] = new ArrayList<>();
                        for (int clus = 0; clus < nOfClusters; clus++) {
                            positions[dim].add((clus + 1) * interval - 0.05 * interval + 0.1 * interval * rand.nextDouble());
                        }
                        Collections.shuffle(positions[dim], rand);
                    }
                    double allRadii[][] = new double[nOfClusters][maxDimension];
                    clusterRadii.put(nOfClusters, allRadii);
                    for (int clus = 0; clus < nOfClusters; clus++) {
                        for (int i = 0; i < maxDimension; i++) {
                            allRadii[clus][i] = 0.2 * interval + 0.8 * rand.nextDouble() * interval;
                        }
                    }
                }

                ArrayList<Double>[] positions = clusterPositions.get(nOfClusters);
                double allRadii[][] = clusterRadii.get(nOfClusters);

                for (double noiseRatio : nOfClusters > 0 ? noiseRatios : zeroClustersNoiseRatios) {

                    for (int dimensionality : nDimensions) {
                        String fileName = basePath + nTuples + "p_" + dimensionality + "d_" + nOfClusters + "c_" + noiseRatio + "n.data";
                        System.out.println("Generating " + fileName);
                        int[] inf = new int[dimensionality];
                        int[] sup = new int[dimensionality];
                        for (int i = 0; i < dimensionality; i++) {
                            inf[i] = 0;
                            sup[i] = domainWidth - 1;
                        }

                        Range r = new Range(inf, sup);

                        int[][] centers = new int[nOfClusters][dimensionality];
                        int[][] radii = new int[nOfClusters][dimensionality];
                        for (int i = 0; i < nOfClusters; i++) {
                            for (int j = 0; j < dimensionality; j++) {
                                centers[i][j] = ((int) (domainWidth * positions[j].get(i) + 0.5));
                                radii[i][j] = (int) (allRadii[i][j] * domainWidth);
                            }
                        }

                        MultidimensionalGaussianGenerator.generate(fileName, r, nTuples, centers, radii, noiseRatio, 100);
                        //leave the following line uncommented if you want a 2d image to be created, otherwise comment it
                        createImage(fileName, fileName + "_labels", r, true);
                        String fileNameIn = basePath + nTuples + "p_" + dimensionality + "d_" + nOfClusters + "c_" + noiseRatio + "n.data";
                        String fileNameOut = basePath2 + nTuples + "p_" + dimensionality + "d_" + nOfClusters + "c_" + noiseRatio + "n.data";
                        shuffleDataset(fileNameIn, fileNameIn + "_labels", fileNameOut, fileNameOut + "_labels", rand);
                    }
                }
            }
        }
    }
}
