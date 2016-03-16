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

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Random;
import java.util.StringTokenizer;
import javax.imageio.ImageIO;

/**
 *
 * @author Giuseppe M. Mazzeo <mazzeo@cs.ucla.edu>
 */
public class MultidimensionalGaussianGenerator {

    static Color[] colbase = {Color.red, Color.green, Color.blue, Color.magenta, Color.cyan, Color.yellow, Color.gray, Color.orange};

    static void shuffle(Random rand, int[] p) {
        for (int i = 0; i < p.length - 1; i++) {
            int j = i + rand.nextInt(p.length - i);
            int temp = p[i];
            p[i] = p[j];
            p[j] = temp;
        }
    }

    public static double getGaussian(Random rand, double mean, double devStandard) {
        return mean + rand.nextGaussian() * devStandard;
    }

    public static void generate(
            String fileName,
            Range domain,
            int numberOfPoints,
            int[][] centers,
            int[][] radii,
            double noiseRatio,
            long seed) throws FileNotFoundException, IOException {

        int nPointsPerCluster = centers.length == 0 ? 0 : (int) (0.5 + numberOfPoints * (1 - noiseRatio) / centers.length);
        int nNoisePoints = numberOfPoints - nPointsPerCluster * centers.length;
        int[] nClusterPoints = new int[centers.length];
        if (nClusterPoints.length > 0) {
            for (int i = 0; i < nClusterPoints.length; i++) {
                nClusterPoints[i] = nPointsPerCluster;
            }
        }
        generate(fileName, domain, centers, radii, nClusterPoints, nNoisePoints, seed);
    }

    public static void generate(
            String fileName,
            Range domain,
            int[][] centers,
            int[][] radii,
            int[] nClusterPoints,
            int nNoisePoints,
            long seed) throws FileNotFoundException, IOException {

        try (
                PrintWriter out = new PrintWriter(new FileOutputStream(fileName), true);
                PrintWriter out2 = new PrintWriter(new FileOutputStream(fileName + "_labels"), true);
                PrintWriter out3 = new PrintWriter(new FileOutputStream(fileName + "_info"), true);) {
            System.out.println("File opened");

            int dimensionality = domain.getDimensionality();

            Random rand = new Random(seed);

            int numberOfDenseRegions = centers.length;

            double[][] devStandard = new double[numberOfDenseRegions][dimensionality];

            for (int i = 0; i < numberOfDenseRegions; i++) {
                for (int j = 0; j < dimensionality; j++) {
                    devStandard[i][j] = radii[i][j] / 4.0; //99% of point will be within radius
                }
            }

            double[][] q = new double[numberOfDenseRegions][dimensionality];
            double[][] s = new double[numberOfDenseRegions][dimensionality];
            int[] n = new int[numberOfDenseRegions];
            double[] sTot = new double[dimensionality];
            int nTot = 0;

            for (int i = 0; i < numberOfDenseRegions; i++) {

                System.out.println("Creating region " + (i + 1) + "...");
                out3.println("Cluster " + (i + 1));
                out3.println("Points: " + nClusterPoints[i]);
                out3.print("Center: " + centers[i][0]);
                for (int j = 1; j < dimensionality; j++) {
                    out3.print("\t" + centers[i][j]);
                }
                out3.println();
                out3.print("Radius: " + radii[i][0]);
                for (int j = 1; j < dimensionality; j++) {
                    out3.print("\t" + radii[i][j]);
                }
                out3.println();
                out3.println();

                int pointsAdded = 0;

                clusterPoints:
                while (pointsAdded < nClusterPoints[i]) {

                    int[] c = new int[dimensionality];

                    for (int k = 0; k < dimensionality; k++) {
                        c[k] = (int) getGaussian(rand, centers[i][k], devStandard[i][k]);
                        if (c[k] < domain.getInfCoord(k + 1) || c[k] > domain.getSupCoord(k + 1)) {
                            continue clusterPoints;
                        }
                    }
                    if (ellipticalRelativeDistance(centers[i], radii[i], c) > 1) {
                        continue;
                    }
                    pointsAdded++;
                    for (int k = 0; k < dimensionality - 1; k++) {
                        out.print(c[k] + ",");
                        q[i][k] += c[k] * c[k];
                        s[i][k] += c[k];
                        sTot[k] += c[k];
                    }
                    out.println(c[c.length - 1]);
                    q[i][c.length - 1] += c[c.length - 1] * c[c.length - 1];
                    s[i][c.length - 1] += c[c.length - 1];
                    sTot[c.length - 1] += c[c.length - 1];
                    n[i]++;
                    nTot++;
                    out2.println(i);

                }
            }

            int nOutlier = 0;
            if (nNoisePoints > 0) {
                System.out.println("Adding noise...");

                //noise tuples are inserted		
                int nearestCluster = -1;
                for (int i = 0; i < nNoisePoints; i++) {
                    int[] p = domain.getRandomInnerCell(rand);
                    int idCluster = clusterContainer(p, centers, radii);
                    if (idCluster == -1) {
                        nOutlier++;
                        //find the nearest cluster
                        double minDist = Double.POSITIVE_INFINITY;
                        for (int r = 0; r < numberOfDenseRegions; r++) {
                            double dist = ellipticalRelativeDistance(centers[r], radii[r], p);
                            if (dist < minDist) {
                                minDist = dist;
                                nearestCluster = r;
                            }
                        }
                    } else {
                        nearestCluster = idCluster;
                    }
                    for (int k = 0; k < dimensionality - 1; k++) {
                        out.print(p[k] + ",");
                    }
                    out.println(p[p.length - 1]);
                    out2.println(idCluster);

                    if (nearestCluster >= 0) { //nearestCluster is -1 only if there are no clusters
                        for (int k = 0; k < dimensionality; k++) {
                            q[nearestCluster][k] += p[k] * p[k];
                            s[nearestCluster][k] += p[k];
                            sTot[k] += p[k];
                        }
                        n[nearestCluster]++;
                        nTot++;
                    }
                }
            }
            System.out.println(nOutlier + " outliers added while adding " + nNoisePoints + " noise points (the rest of noise was adsorbed by clusters)");
            double wc = 0;
            double bc = 0;
            for (int k=0; k<dimensionality; k++) {
                sTot[k]/=nTot;
            }
            for (int i = 0; i < numberOfDenseRegions; i++) {
                double pi=1.0*n[i]/nTot;
                double wci=0;
                for (int k=0; k<dimensionality; k++) {
                    wci+=q[i][k]-(s[i][k]/n[i])*s[i][k];
                    s[i][k]/=n[i];
                }
                double bci=pi*sqrDistance(sTot, s[i]);
                wci/=nTot;
                wc+=wci;
                bc+=bci;
            }
            if (wc!=0) {
                out3.println("VarianceRatioClusterability: "+(bc/wc));
            }
        }
    }

    private static int clusterContainer(int[] p, int[][] centers, int[][] radius) {
        double dMin = 1;
        int iMin = -1;
        for (int i = 0; i < centers.length; i++) {
            double d = ellipticalRelativeDistance(centers[i], radius[i], p);
            if (d < dMin) {
                iMin = i;
                dMin = d;
            }
        }
        return iMin;
    }

    public static double ellipticalRelativeDistance(int[] center, int[] radius, int[] p) {
        double d = 0;
        for (int i = 0; i < radius.length; i++) {
            d += Math.pow((1.0 * center[i] - p[i]) / radius[i], 2);
        }
        return d;
    }
    
    public static double sqrDistance(double[] p1, double[] p2) {
        double d = 0;
        for (int i = 0; i < p1.length; i++) {
            d += Math.pow(p1[i] - p2[i], 2);
        }
        return d;
    }
    

    public static void createImage(String fileName, String fileNameLabels, Range domain, boolean colors) throws IOException {

        BufferedReader in = new BufferedReader(new FileReader(fileName));
        BufferedReader inlabels = new BufferedReader(new FileReader(fileNameLabels));
        HashMap<Integer, Color> color = new HashMap<>();
        color.put(-1, Color.BLACK);
        int infX = domain.getInfCoord(1);
        int supX = domain.getSupCoord(1);
        int infY = domain.getInfCoord(2);
        int supY = domain.getSupCoord(2);
        int width = supX - infX + 1;
        int height = supY - infY + 1;
        Graphics2D pianoG;
        BufferedImage piano;
        int imageType = BufferedImage.TYPE_INT_RGB;
        piano = new BufferedImage(width, height, imageType);
        pianoG = piano.createGraphics();
        pianoG.setColor(Color.WHITE);
        pianoG.fillRect(0, 0, width, height);
        pianoG.setColor(Color.BLACK);
        pianoG.drawRect(0, 0, width - 1, height - 1);

        int nextColor = 0;
        String l = in.readLine();
        String ll = inlabels.readLine();
        while (l != null) {
            StringTokenizer st = new StringTokenizer(l, " ,");
            int x = Integer.parseInt(st.nextToken());
            int y = Integer.parseInt(st.nextToken());
            int label = Integer.parseInt(ll);
            Color c = color.get(label);
            if (c == null) {
                c = colbase[nextColor];
                nextColor = (nextColor + 1) % colbase.length;
                color.put(label, c);
            }
            pianoG.setColor(colors ? c : Color.BLACK);
            pianoG.drawLine(x, y, x, y);
            l = in.readLine();
            ll = inlabels.readLine();
        }
        ImageIO.write(piano, "PNG", new File(fileName + ".png"));
    }

    public static void shuffleDataset(String fileDataIn, String fileLabelsIn, String fileDataOut, String fileLabelsOut, Random rand) throws IOException {
        ArrayList<String> data = new ArrayList<>();
        ArrayList<String> labels = new ArrayList<>();
        try (
                BufferedReader in1 = new BufferedReader(new FileReader(fileDataIn));
                BufferedReader in2 = new BufferedReader(new FileReader(fileLabelsIn));) {
            String l1 = in1.readLine();
            String l2 = in2.readLine();
            while ((l1 != null && l1.length() > 0) && (l2 != null && l2.length() > 0)) {
                data.add(l1);
                labels.add(l2);
                l1 = in1.readLine();
                l2 = in2.readLine();
            }
            if (l1 != null && l1.length() > 0 || l2 != null && l2.length() > 0) {
                System.out.println("File lengths mismatch");
                System.exit(0);
            }
        }

        for (int i = 0; i < data.size() - 1; i++) {
            int j = i + rand.nextInt(data.size() - i);
            String temp = data.get(i);
            data.set(i, data.get(j));
            data.set(j, temp);
            temp = labels.get(i);
            labels.set(i, labels.get(j));
            labels.set(j, temp);
        }

        try (
                PrintWriter out1 = new PrintWriter(new FileWriter(fileDataOut), true);
                PrintWriter out2 = new PrintWriter(new FileWriter(fileLabelsOut), true);) {
            for (int i = 0; i < data.size() - 1; i++) {
                out1.println(data.get(i));
                out2.println(labels.get(i));
            }
        }
    }
}
