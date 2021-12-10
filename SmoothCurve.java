/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package smoothcurve;

import com.sun.org.apache.bcel.internal.generic.SWAP;
import java.io.*;
import java.util.ArrayList;
import java.util.Random;
import java.util.Scanner;
import java.io.FileWriter;   // Import the FileWriter class
import java.io.IOException;  // Import the IOException class to handle errors

public class SmoothCurve {

    static ArrayList<ArrayList> ArrOfChromosomes = new ArrayList<>();

    static ArrayList<point> ArrOfpoints = new ArrayList<point>();

    static Random random = new Random();
    static int Iteration = 10;
    static int dependencyFactor = 5;
    static int totalGenerations = 5000;
    static int globalDegree = 0;
    static int globalPopSize = 0;

    static public class point {

        float x, y;

        public point(float x1, float y1) {
            x = x1;
            y = y1;
        }

    }

    //returns a unique random number within the limit, since it has all the previous random numbers
    static public int generateRandom(int limit, ArrayList<Integer> prevChoices) {

        int rand = random.nextInt(limit);
        while (prevChoices.contains(rand)) {
            rand = random.nextInt(limit);
        }
        return rand;

    }

    static public void CreateChromosomes(ArrayList<point> ArrofPoints, int degree, int populationSize) { //initialize a population

        for (int i = 0; i < populationSize; i++) {

            ArrayList<Double> Chromosome = new ArrayList<>();
            for (int j = 0; j < degree + 1; j++) {
                double Coeff = (-10) + (10 - (-10)) * random.nextDouble();
                Chromosome.add(Coeff);
            }
            ArrOfChromosomes.add(Chromosome);
        }

    }

    static public ArrayList<Double> calculateFitness(ArrayList<ArrayList> ArrOfChromosomes, int degree, int populationSize, ArrayList<point> ArrofPoints) {
        ArrayList<Double> ArrOfFitness = new ArrayList<>();
        for (int i = 0; i < populationSize; i++) {
            double TotalError = 0.0;
            for (int k = 0; k < ArrofPoints.size(); k++) //1 5
            {
                double CalError = 0.0;
                double Error = 0.0;
                for (int a = 0; a < degree + 1; a++) {
                    CalError += (Math.pow(ArrofPoints.get(k).x, a)) * (double) (ArrOfChromosomes.get(i).get(a));  //1^0 * 0.2  + //1^1 * 0.3 + //1^2 *0.5

                }

                Error = CalError - (ArrofPoints.get(k).y);

                Error *= Error;

                TotalError += Error;

            }
            TotalError /= ArrofPoints.size();
            ArrOfFitness.add(TotalError);
        }

        return ArrOfFitness;

    }

    static public ArrayList<Integer> CreateSelection(int PopulationSize, ArrayList<Double> ArrOfFitness, ArrayList<ArrayList> ArrOfChromosome, int degree) {
        ArrayList<Integer> Selected = new ArrayList<>();

        int MatingPoolSize = PopulationSize;

        ArrayList<Integer> chosenindex = new ArrayList<>();
        for (int i = 0; i < MatingPoolSize; i++) {
            int firstnum = generateRandom(PopulationSize, chosenindex);
            int secondnum = generateRandom(PopulationSize, chosenindex);

            if (ArrOfFitness.get(firstnum) <= ArrOfFitness.get(secondnum)) {
                chosenindex.add(firstnum);
                Selected.add(firstnum);

            } else {
                chosenindex.add(secondnum);
                Selected.add(secondnum);
            }

        }

        return Selected;

    }

    static public ArrayList<ArrayList> CreateCrossOver(ArrayList<Integer> Selected, int degree, ArrayList<ArrayList> ArrOfChromosome) {
        ArrayList<ArrayList> CrossedOver = new ArrayList<>();

        double pc = 0.7;
        int FirstChromosome = 0;
        int SecChromosome = 1;
        for (int i = 0; i < Selected.size(); i += 2) {
            ArrayList<Double> OffSpring1 = new ArrayList<>();
            ArrayList<Double> OffSpring2 = new ArrayList<>();
            double r2 = 0 + (1 - 0) * random.nextDouble();
            if (r2 <= pc) {

                int point1 = random.nextInt(degree) + 1;
                int point2 = random.nextInt(degree) + 1;
                while (point1 == point2) { //to get 2 different points
                    point1 = random.nextInt(degree) + 1;
                    point2 = random.nextInt(degree) + 1;
                }

                int temp;
                if (point1 > point2) { //because points should be in an ascending order to perform a 2-point crossover
                    temp = point1;
                    point1 = point2;
                    point2 = temp;
                }
                for (int j = 0; j < point1; j++) //2
                {
                    OffSpring1.add((double) ArrOfChromosome.get(Selected.get(i)).get(j));  //0.2  0.3
                    OffSpring2.add((double) ArrOfChromosome.get(Selected.get(i + 1)).get(j)); //0.6 0.9

                }
                for (int j = point1; j < point2; j++) {
                    OffSpring1.add((double) ArrOfChromosome.get(Selected.get(i + 1)).get(j));  //0.8 
                    OffSpring2.add((double) ArrOfChromosome.get(Selected.get(i)).get(j));   //0.4
                }
                for (int j = point2; j < degree + 1; j++) {
                    OffSpring1.add((double) ArrOfChromosome.get(Selected.get(i)).get(j));  //0.5
                    OffSpring2.add((double) ArrOfChromosome.get(Selected.get(i + 1)).get(j));   //0.7
                }
                CrossedOver.add(OffSpring1);
                CrossedOver.add(OffSpring2);

            } else {

                CrossedOver.add(ArrOfChromosome.get(Selected.get(i)));
                CrossedOver.add(ArrOfChromosome.get(Selected.get(i + 1)));

            }

        }

        return CrossedOver;
    }

    static public ArrayList<ArrayList> CreateMutation(ArrayList<ArrayList> CrossedOver, int degree, int currentIteration, int totalIterations) {
        ArrayList<ArrayList> Mutated = new ArrayList<>();
        double pm = 0.1;
        for (int i = 0; i < CrossedOver.size(); i++) {
            ArrayList<Double> Chromosome = new ArrayList<>();

            for (int j = 0; j < degree + 1; j++) {

                double isMutating = 0 + (1 - 0) * random.nextDouble();
                if (isMutating <= pm) {

                    double r = 0 + (1 - 0) * random.nextDouble();
                    int change;
                    double DeltaL = (double) CrossedOver.get(i).get(j) + 10;
                    double DeltaU = 10 - ((double) CrossedOver.get(i).get(j));
                    double deltaTandY;
                    double y;
                    double r2 = 0 + (1 - 0) * random.nextDouble();
                    if (r2 <= 0.5) {
                        y = DeltaL;
                        change = -1;
                    } else {

                        y = DeltaU;
                        change = 1;
                    }

                    double currentOverTotal = (double) (currentIteration) / totalIterations;
                    double powerOfR = (double) Math.pow((1 - currentOverTotal), dependencyFactor);
                    double R = (double) Math.pow(r, powerOfR);

                    deltaTandY = y * (1 - R);
                    double newx;

                    if (change == -1) {
                        newx = (double) CrossedOver.get(i).get(j) - deltaTandY;
                    } else {
                        newx = (double) CrossedOver.get(i).get(j) + deltaTandY;
                    }
                    Chromosome.add(newx);

                } else {
                    Chromosome.add((double) CrossedOver.get(i).get(j));
                }
            }
            Mutated.add(Chromosome);
        }

        return Mutated;
    }

    static public ArrayList<ArrayList> Replacement(ArrayList<ArrayList> ArrOfChromosomes, ArrayList<ArrayList> Mutated, int degree, int populationSize, ArrayList<Integer> Selected, ArrayList<Double> ArrOfFitness, ArrayList<Double> ArrOfNewFitness) {

        ArrayList<ArrayList> ArrOfNewGeneration = new ArrayList<>();

        for (int i = 0; i < Selected.size(); i++) {
            ArrayList<Double> MovedGene = new ArrayList<>();
            if (ArrOfFitness.get(Selected.get(i)) > ArrOfNewFitness.get(i)) { //if an offspring has less fitness than its parent

                for (int j = 0; j < degree + 1; j++) {

                    MovedGene.add((double) Mutated.get(i).get(j));

                }

                ArrOfNewGeneration.add(MovedGene);

            } else {
                for (int j = 0; j < degree + 1; j++) {

                    MovedGene.add((double) ArrOfChromosomes.get(Selected.get(i)).get(j));

                }
                ArrOfNewGeneration.add(MovedGene);
            }

        }

        return ArrOfNewGeneration;
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws Exception {

        int numberoftestcases = 0;
        File file = new File("input-2.txt");
        Scanner reader = new Scanner(file);
        numberoftestcases = reader.nextInt();
        for (int i = 0; i < numberoftestcases; i++) {
            ArrayList<Double> MinOfFitness = new ArrayList<>();
            ArrayList<ArrayList> ChosenChromosomes = new ArrayList<>();
            int numberofpoints = 0;
            //to fill them again for each test case
            ArrOfpoints.clear();
            ArrOfChromosomes.clear();
            numberofpoints = reader.nextInt();
            int degree = reader.nextInt();
            globalDegree = degree;
            for (int j = 0; j < numberofpoints; j++) {
                point pointt = new point(reader.nextFloat(), reader.nextFloat());
                ArrOfpoints.add(pointt);
            }
            int populationSize;
            if (numberofpoints % 2 == 0) {
                populationSize = numberofpoints + 2;
            } else {
                populationSize = numberofpoints + 3;
            }

            globalPopSize = populationSize;

            CreateChromosomes(ArrOfpoints, degree, populationSize); //init population
            for (int b = 0; b < totalGenerations; b++) {

                ArrayList<Double> ArrOfFitness = calculateFitness(ArrOfChromosomes, globalDegree, globalPopSize, ArrOfpoints);
                ArrayList<Integer> selected = CreateSelection(globalPopSize, ArrOfFitness, ArrOfChromosomes, globalDegree);
                ArrayList<ArrayList> CrossedOver = CreateCrossOver(selected, globalDegree, ArrOfChromosomes);
                ArrayList<ArrayList> Mutated = CreateMutation(CrossedOver, globalDegree, b, totalGenerations);

                ArrayList<Double> ArrOFNewFitness = calculateFitness(Mutated, globalDegree, globalPopSize, ArrOfpoints);

                ArrOfChromosomes = Replacement(ArrOfChromosomes, Mutated, globalDegree, globalPopSize, selected, ArrOfFitness, ArrOFNewFitness);
                double min = ArrOFNewFitness.get(0);
                int indexOfmin = 0;
                for (int k = 1; k < ArrOFNewFitness.size(); k++) {
                    if (min > ArrOFNewFitness.get(k)) {
                        min = ArrOFNewFitness.get(k);
                        indexOfmin = k;
                    }

                }
                
                //getting the minimum fitness in an iteration 
                MinOfFitness.add(min);
                ArrayList<Double> MinChromosome = new ArrayList<>();
                for (int ind = 0; ind < degree + 1; ind++) {
                    MinChromosome.add((Double) ArrOfChromosomes.get(selected.get(indexOfmin)).get(ind));

                }
                ChosenChromosomes.add(MinChromosome);
            }
            
            //getting the minimum from all iterations and saving its index to add it to the file accompanied with its fitness
            double FinalError = MinOfFitness.get(0);
            int indexOfFinalError = 0;
            for (int d = 1; d < MinOfFitness.size(); d++) {
                if (FinalError > MinOfFitness.get(d)) {
                    FinalError = MinOfFitness.get(d);
                    indexOfFinalError = d;
                }
            }

            try {

                FileWriter fw = new FileWriter("Out.txt", true);
                BufferedWriter bw = new BufferedWriter(fw);
                PrintWriter pw = new PrintWriter(bw);

                for (int indx = 0; indx < degree + 1; indx++) {
                    pw.println((Double) ChosenChromosomes.get(indexOfFinalError).get(indx));

                }
                pw.println("with FinalError = " + FinalError);
                pw.flush();
                pw.close();
                bw.close();
                fw.close();
            } catch (IOException e) {
                System.out.println("An error occurred.");

            }

        }

    }
}
