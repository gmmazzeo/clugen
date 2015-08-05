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

import java.util.Random;

/**
 *
 * @author Giuseppe M. Mazzeo <mazzeo@cs.ucla.edu>
 */
public class Range {

    protected int[] inf, sup;

    public Range(int[] inf, int[] sup) {
        this.inf = inf;
        this.sup = sup;
    }

    public int getInfCoord(int i) {
        return inf[i - 1];
    }

    public int getSupCoord(int i) {
        return sup[i - 1];
    }

    public int[] getInf() {
        return inf;
    }

    public int[] getSup() {
        return sup;
    }

    public void setInfCoord(int i, int v) {
        inf[i - 1] = v;
    }

    public void setSupCoord(int i, int v) {
        sup[i - 1] = v;
    }

    public double getSize() {
        double v = 1;
        for (int i = 0; i < inf.length; i++) {
            v *= sup[i] - inf[i] + 1;
        }
        return v;
    }

    public boolean contains(int[] p) {
        for (int i = 0; i < inf.length; i++) {
            if (p[i] < inf[i] || p[i] > sup[i]) {
                return false;
            }
        }
        return true;

    }

    public int[] getRandomInnerCell(Random rnd) {
        int[] coord = new int[inf.length];
        for (int k = 0; k < coord.length; k++) {
            coord[k] = inf[k] + rnd.nextInt(sup[k] - inf[k] + 1);
        }
        return coord;
    }

    public int getDimensionality() {
        return inf.length;
    }

    @Override
    public String toString() {
        String s = "";
        for (int i = 0; i < inf.length; i++) {
            s += inf[i] + ".." + sup[i] + " ";
        }
        return s;
    }
}
